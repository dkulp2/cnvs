# staircase.R - modify a set of CNV predictions by removing small intermediate segments of a "staircase"
#               then adjust the boundaries using an MLE strategy
#
# Author: David Kulp, dkulp@broadinstitute.org
#
# Usage: Called using rscript:
#
# #1: filename base - e.g. sites_cnv_segs.txt
# #2: method - which set of segments to read? ('csm','ncsm','mlcsm'). File read is base.method.Rdata
# #3: db connection - user:host:port:dbname - profile data is read from this connection
# #4: the max size of intermediate regions that was skipped when joining (in bins, e.g. 5). Use 4x this value to remove short NA segments
# #5: the window size that was used for aggregating profiles (in bins) - in staircase, don't remove intermediate CNVs larger than this 
# #6: the size of the window (in bins) when scanning for boundary using MLE method (e.g. 10 bins)
# #7: data set label

library(plyr)
library(dplyr)
library(RPostgreSQL)

DEBUG <- FALSE

cmd.args <- commandArgs(trailingOnly = TRUE)
# Sys.setenv(PGHOST="localhost",PGUSER="dkulp",PGDATABASE="seq", PGOPTIONS="--search_path=data_sfari_batch1c_27apr2017")
# cmd.args <- c('/cygwin64/home/dkulp/data/SFARI.27April2017/dataC/sites_cnv_segs.txt','csm','5','12','10','data_sfari_batch1C_27Apr2017')
cnv.seg.fn <- cmd.args[1]
cnv.seg.method <- cmd.args[2]
max.join <- as.numeric(cmd.args[3])  # size (in bins) of runs to span if flanking calls are the same, e.g. 5
profile.group.size <- as.numeric(cmd.args[4])  # length of window size in bins, e.g. 12
win.size.bins <- as.numeric(cmd.args[5])  # in bins. e.g. 10 for a ~1000 nt window.
data.label <- cmd.args[6]


load(sprintf("%s.%s.Rdata",cnv.seg.fn,cnv.seg.method)) # => cn.segs.merged
csm <- as.tbl(cn.segs.merged)
message(Sys.time(), sprintf(": Loaded %s rows", nrow(csm)))
csm$len <- csm$seg <- csm$len2 <- csm$gap <- NULL
csm$len <- csm$end.map - csm$start.map 
csm$len.bin <- csm$end.bin - csm$start.bin + 1

# identify "small" NAs and remove them
csm <- filter(csm, !is.na(cn) | len.bin > 4*max.join)
message(Sys.time(), ": Removed small NAs <= ",4*max.join,", leaving ",nrow(csm)," rows")

# merge rows with same CN now that small NAs are removed
csm <- ddply(csm, .(.id, chr), function(df) {
  # the first row of a pair of adjacent rows with the same CN
  idx <- 1:(nrow(df)-1)
  adj <- df$cn[idx] == df$cn[idx+1]
  while (any(adj, na.rm=TRUE)) {
    adj.idx <- first(which(adj))
    df[adj.idx, 'end.i'] <- df[adj.idx+1, 'end.i']
    df[adj.idx, 'end.map'] <- df[adj.idx+1, 'end.map']
    df[adj.idx, 'end.bin'] <- df[adj.idx+1, 'end.bin']
    df[adj.idx, 'label'] <- paste0(df[adj.idx, 'label'], '_', df[adj.idx+1, 'label'])
    df <- df[-c(adj.idx+1),]
    idx <- 1:(nrow(df)-1)
    adj <- df$cn[idx] == df$cn[idx+1]
  } 
  df$len <- df$end.map - df$start.map 
  df$len.bin <- df$end.bin - df$start.bin + 1
  return(df)
})
message(Sys.time(),": Merged now-adjacent segments with same CN, leaving ",nrow(csm)," rows")

# identify candidate stairstep regions
m.prev.idx <- 1:(nrow(csm)-2)
m.idx <- 2:(nrow(csm)-1)
m.next.idx <- 3:(nrow(csm))

# candidates holds the indices into csm of the middle region of a up or down 1-copy difference staircase
# NOTE: assume csm is already ordered by sample and fragment position
candidates <- m.idx[csm$len.bin[m.idx]<profile.group.size &
                      (csm$cn[m.idx] - csm$cn[m.prev.idx] == csm$cn[m.next.idx] - csm$cn[m.idx]) & 
                      (abs(csm$cn[m.next.idx] - csm$cn[m.prev.idx]) == 2) &
                      csm$.id[m.prev.idx] == csm$.id[m.next.idx]]

# NAs are where one of the segments are NA
candidates <- na.omit(candidates)

# remove the intermediate candidates
message(Sys.time(),": Removing ",length(candidates)," small segments that are in a -1, 0, +1 staircase")
if (length(candidates)>0) {
  csm <- csm[-candidates,]
}
message(Sys.time(),": leaving ",nrow(csm)," rows")


###################################################################################################################

# connect to DB
message(Sys.time(),": Connecting to db.")
db <- src_postgres()

# retrieve them all into memory so it's easy to do a bin=>pos mapping
profile.segments <- tbl(db, 'profile_segment') %>% collect(n=Inf)
# load('/cygwin64/tmp/profile_segments.Rdata')

# calculate the CI by growing greadily away from max. There are no NAs in p.
conf.int <- function(p, conf=0.95) {
  best.pos <- which.max(p)
  i <- best.pos - 1
  j <- best.pos + 1
  mass <- p[best.pos]
  while (mass < conf && (j <= length(p) || i > 0)) {
    if (i > 0) {
      if (j > length(p) || p[i] >= p[j]) {
        mass <- mass + p[i]
        i <- i - 1
      }
    } 
    if (j <= length(p) && mass < conf) {
      if (i == 0 || p[j] > p[i]) {
        mass <- mass + p[j]
        j <- j + 1
      }
    }
  }
  return(list(left=i-best.pos+1, right=j-best.pos-1))
}


# for each sample, load all profile data, for each transition compute MLE boundary and confidence intervals
# as a side effect, also write the probability of the data over all possible transitions
csm.new <- ddply(csm, .(.id), function(df) {
  # df <- csm[csm$.id=='SSC00115',]
  # df <- csm[csm$.id=='SSC02692',]
  sample <- df$.id[1]
  message(Sys.time(),": Computing MLE boundaries for ",sample)
  
  pois.cn.df <- dbGetQuery(db$con, sprintf("SELECT ps.chrom, ps.start_pos, ps.end_pos, pois.* FROM pois, profile_segment ps WHERE pois.label='%s' AND pois.sample='%s' AND pois.bin=ps.bin AND pois.chr=ps.chrom ORDER BY ps.chrom, ps.start_pos", data.label, sample))
  # save(pois.cn.df, file='/cygwin64/tmp/pois_cn_df.SSC02692.Rdata')
  # load('/cygwin64/tmp/pois_cn_df.SSC02692.Rdata')
  
  stopifnot(win.size.bins %% 2 == 0) # must be even number of bins
  half.win <- win.size.bins / 2 # each window is divided into 2 equal sides for cnA and cnB

  pois.cnL <- as.list(pois.cn.df[,grep("cnL",names(pois.cn.df))])
  pois.cnR <- as.list(pois.cn.df[,grep("cnR",names(pois.cn.df))])

  
  ml.transition2 <- function(sample, bin1, bin2, cnA, cnB) {

    binToPos <- function(b, side='start_pos') {
      pos <- filter(pois.cn.df, bin==b)[[side]]
      if (length(pos)==0) { warning(sprintf("Bin=>Pos mapping failed for %s. sample=%s, bin1=%s, bin2=%s", b, sample, bin1, bin2)) }
      pos
    }
    
    if (is.na(cnA)) {
      return(data.frame(pos=binToPos(bin2), left.bound=NA, right.bound=NA, win.size=NA, left.tail=NA, right.tail=NA, best.bin=bin2, binL=NA, binR=NA, binCI.L=NA, binCI.R=NA))
    } else if (is.na(cnB)) {
      return(data.frame(pos=binToPos(bin1), left.bound=NA, right.bound=NA, win.size=NA, left.tail=NA, right.tail=NA, best.bin=bin1, binL=NA, binR=NA, binCI.L=NA, binCI.R=NA))
    }

    cnA <- ifelse(cnA > 3, 4, cnA+1)  # map CN to index into pois.cn. CN > 3 => CN:=3
    cnB <- ifelse(cnB > 3, 4, cnB+1)

    binL <- bin1 - win.size.bins
    binL <- ifelse(binL < 1, 1, binL)
    posL <- binToPos(binL)
    
    binR <- bin2 + win.size.bins
    binR <- ifelse(binR > max(pois.cn.df$bin), max(pois.cn.df$bin), binR)
    posR <- binToPos(binR,'end_pos')

    bin.count <- binR - binL + 1
    pois.idx.start <- which(pois.cn.df$bin==binL); stopifnot(length(pois.idx.start)>0)
    pA <- pois.cnL[[cnA]][pois.idx.start:(pois.idx.start+bin.count-win.size.bins)] 
    pB <- pois.cnR[[cnB]][pois.idx.start:(pois.idx.start+bin.count-win.size.bins)]
    jp <- pA * pB    # joint likelihood
    jp.norm <- jp / sum(jp)  # normalized
    best.jp.bin <- which.max(jp.norm)

    #    cat(sprintf("%s:%.0f/%.0f => %.0f (%.4f..%.4f)\n", sample, pos1, pos2, best.pos, min(-log(jp.norm)), max(-log(jp.norm))))
    if (length(best.jp.bin)>0) { # Either cnA or cnB is NA if length is zero
      best.bin <- binL+(best.jp.bin-1)+half.win # first bin passed the transition
      best.pos <- binToPos(best.bin)

      CI <- conf.int(jp.norm)  # returns CI$left and CI$right, which are bin offsets of best.pos
      binCI.L <- best.bin + CI$left
      binCI.R <- best.bin + CI$right
      
      # debug plots
      if (DEBUG) {
        wR <- (pois.idx.start-10):(pois.idx.start+bin.count+10-win.size.bins)
        LA <- pois.cnL[[cnA]][wR]
        RB <- pois.cnR[[cnB]][wR]
        AB <- LA*RB
        WB <- pois.cn.df$bin[wR]+half.win
        gp <- rbind(data.frame(bin=WB, p=LA, var=as.character(cnA-1),part='A'),
                    data.frame(bin=WB, p=RB, var=as.character(cnB-1),part='A'),
                    data.frame(bin=WB, p=AB, var='AB',part='B'))
        b <- data.frame(x=c(best.bin,binL,binR, binL+half.win, binR-half.win, bin1, bin2), 
                        var=c('best','bounds','bounds','midbounds','midbounds','starts','starts'),
                        part=c('B','A','A','A','A','A','A'))
        maxp <- max(AB)
        print(ggplot(gp, aes(x=bin, y=p, color=var)) + geom_point() + 
                facet_grid(part~., scales="free_y") + 
                geom_vline(aes(xintercept=x, color=var), data=b) + 
                geom_segment(x=binCI.L,xend=binCI.R,y=maxp,yend=maxp))
      }
      
      # bounds are generous, including the furthest base from the best.pos in the transition bins.
      left.bound <- binToPos(binCI.L)
      right.bound <- binToPos(binCI.R, 'end_pos')
      
      return(data.frame(pos=best.pos, left.bound=left.bound, right.bound=right.bound, win.size=posR-posL+1, left.tail=jp.norm[1], right.tail=jp.norm[length(jp.norm)], binL=binL, binR=binR, best.bin=best.bin, binCI.L=binCI.L, binCI.R=binCI.R))
    } else {
      message(Sys.time(),": Failed (%s %s %s %s %s)", sample, bin1, bin2, cnA, cnB)
      best.bin <- bin1 + (bin2-bin1)%/%2 + 1
      best.pos <- binToPos(best.bin)
      return(data.frame(pos=best.pos, left.bound=NA, right.bound=NA, win.size=posR-posL+1, left.tail=NA, right.tail=NA, best.bin=best.bin, binL=binL, binR=binR, binCI.L=NA, binCI.R=NA))
    }
  }
  
  # Iterate across every boundary, where the boundary between the first and second segment is idx=1, etc.
  new.bounds <- ldply(1:(nrow(df)-1), function(idx) { 
    ml.transition2(sample, df$end.bin[idx], df$start.bin[idx+1], df$cn[idx], df$cn[idx+1]) 
  })

  # replace old values of start.map, end.map, start.bin, end.bin and add interval data.
  nrm1 <- nrow(df)-1
  df$end.map[1:nrm1] <- new.bounds$pos
  df$end.map.L <- c(new.bounds$left.bound,NA)
  df$end.map.R <- c(new.bounds$right.bound,NA)
  df$end.map.win.size <- c(new.bounds$win.size,NA)
  df$end.map.L.tail <- c(new.bounds$left.tail,NA)
  df$end.map.R.tail <- c(new.bounds$right.tail,NA)
  df$end.bin.L <- c(new.bounds$binL,NA)
  df$end.bin.R <- c(new.bounds$binR,NA)
  df$end.bin[1:nrm1] <- new.bounds$best.bin
  df$end.binCI.L <- c(new.bounds$binCI.L,NA)
  df$end.binCI.R <- c(new.bounds$binCI.R,NA)
  
  df$start.map[1:nrm1+1] <- new.bounds$pos
  df$start.map.L <- c(NA,new.bounds$left.bound)
  df$start.map.R <- c(NA,new.bounds$right.bound)
  df$start.map.win.size <- c(NA,new.bounds$win.size)
  df$start.map.L.tail <- c(NA,new.bounds$left.tail)
  df$start.map.R.tail <- c(NA,new.bounds$right.tail)
  df$start.bin.L <- c(NA,new.bounds$binL)
  df$start.bin.R <- c(NA,new.bounds$binR)
  df$start.bin[1:nrm1+1] <- new.bounds$best.bin
  df$start.binCI.L <- c(NA,new.bounds$binCI.L)
  df$start.binCI.R <- c(NA,new.bounds$binCI.R)
  
  # remove unneeded columns
  df$start.i <- df$end.i <- NULL

  # add an initial and final segment
  this.chr <- first(df$chr)
  profile.segments.chr <- filter(profile.segments, chrom==this.chr)
  df <- rbind(tibble(.id=sample, cn=2, chr=first(df$chr), start.map=1, end.map=df$start.map[1],
                     start.bin=1, end.bin=df$start.bin[1], copy.number=factor(2),
                     label=sprintf("SEG_%s_%s_START",first(df$chr), sample), 
                     len=df$start.map[1], len.bin=df$start.bin[1],
                     end.map.L=NA, end.map.R=NA, end.map.win.size=NA, end.map.L.tail=NA, end.map.R.tail=NA,
                     end.bin.L=NA, end.bin.R=NA, end.binCI.L=NA, end.binCI.R=NA,
                     start.map.L=NA, start.map.R=NA, start.map.win.size=NA, 
                     start.map.L.tail=NA, start.map.R.tail=NA, start.bin.L=NA, start.bin.R=NA,
                     start.binCI.L=NA, start.binCI.R=NA),
              df,
              tibble(.id=sample, cn=2, chr=first(df$chr), start.map=df$end.map[nrow(df)], 
                     end.map=max(profile.segments.chr$end_pos),
                     start.bin=df$end.bin[nrow(df)], end.bin=max(profile.segments.chr$bin), 
                     copy.number=factor(2),
                     label=sprintf("SEG_%s_%s_END",first(df$chr), sample), 
                     len=max(profile.segments.chr$end_pos)-df$end.map[nrow(df)], 
                     len.bin=max(profile.segments.chr$bin)-df$end.bin[nrow(df)],
                     end.map.L=NA, end.map.R=NA, end.map.win.size=NA, end.map.L.tail=NA, end.map.R.tail=NA,
                     end.bin.L=NA, end.bin.R=NA, end.binCI.L=NA, end.binCI.R=NA,
                     start.map.L=NA, start.map.R=NA, start.map.win.size=NA, 
                     start.map.L.tail=NA, start.map.R.tail=NA, start.bin.L=NA, start.bin.R=NA,
                     start.binCI.L=NA, start.binCI.R=NA))

  # It is possible that small CNVs have breakpoint bounds that have reversed order because the window scans beyond the size of the CNV.
  # Remove those CNVs.
  # For example:
  # Previously: 11111111111111333322222222222
  # Now:        11111111111111111111
  #                              333
  #                              222222222222
  #
  # Change to:
  #             11111111111111111NNN222222222
  if (any(df$end.bin < df$start.bin)) {
    reversed <- which(df$end.bin < df$start.bin)
    message(Sys.time(), ": Removing ",length(reversed), " reversed CNVs and adjusting flanks.")
    df[reversed-1,'end.map'] <- df[reversed,'end.map']
    df[reversed+1,'start.map'] <- df[reversed,'start.map']
    df[reversed,'start.map'] <- df[reversed-1,'end.map']
    df[reversed,'end.map'] <- df[reversed+1,'start.map']

    df[reversed-1,'end.bin'] <- df[reversed,'end.bin']
    df[reversed+1,'start.bin'] <- df[reversed,'start.bin']
    df[reversed,'start.bin'] <- df[reversed-1,'end.bin']
    df[reversed,'end.bin'] <- df[reversed+1,'start.bin']
    
    df[reversed,'cn'] <- NA
  }

  # It is possible that multiple CNVs start at the same breakpoint or end at the same breakpoint
  # (There should be no overlaps since that should have been resolved in the previous step.)
  # Set all such segments to CN=NA since it suggests an ambiguous site.
  idx <- 1:(nrow(df)-1)

  # Check that all adjacent segments are consistent
  stopifnot(all(df$end.bin[idx] == df$start.bin[idx+1]))

  ## # row numbers of left CNV of conflicting breakpoints
  ## trouble <- which(df$end.bin[idx] > df$start.bin[idx+1])

  ## if (length(trouble)>0) {
  ##   trouble <- unique(c(trouble+1, trouble))
  ##   message(Sys.time(), sprintf(": Eliminating %s conflicting segments after MLE adjustments.", length(trouble)))

  ##   if (!all(df$start.bin[trouble]==df$start.bin[trouble+1] | df$end.bin[trouble]==df$end.bin[trouble+1])) {
  ##     message(Sys.time(), sprintf(": WARNING. Overlapping segments for %s that don't share common start or end point. Shouldn't happen."))
  ##   }

  ##   df$cn[trouble] <- NA
  ## }

  df$dCN.R <- c(df$cn[2:nrow(df)]-df$cn[1:(nrow(df)-1)],0)
  df$dCN.L <- c(0, df$dCN.R[1:nrow(df)-1])
  df$dL <- ifelse(df$dCN.L > 0, 'G', ifelse(df$dCN.L < 0, 'L', 'N'))
  df$dR <- ifelse(df$dCN.R > 0, 'G', ifelse(df$dCN.R < 0, 'L', 'N'))

  return(df)
})

cn.segs.merged <- mutate(csm.new, 
                         seg=sprintf("SEG_%s_%s_%s", chr, start.map, end.map),
                         label=sprintf("%s_%s",seg,.id))

message(Sys.time(), ": Saving smlcsm")
save(cn.segs.merged, file=sprintf("%s.smlcsm.Rdata",cnv.seg.fn))
write.table(select(cn.segs.merged, .id, seg, chr, start.map, end.map, copy.number), file=sprintf("%s.sml.tbl",cnv.seg.fn), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(select(cn.segs.merged, .id, seg, chr, start.map.L, as.integer(start.map), start.map.R, end.map.L, as.integer(end.map), end.map.R, copy.number), file=sprintf("%s.smlCI.tbl",cnv.seg.fn), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(cn.segs.merged[!is.na(cn.segs.merged$cn),c('chr','start.map','end.map','label','cn')], file=sprintf("%s.smlmrg.bed",cnv.seg.fn), sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)

rm(db)

