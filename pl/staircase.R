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
  # df <- csm[csm$.id=='SSC03794',]
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
  # all these ifelse statements simply replace rows only where cn is not NA.
  # If we worked in breakpoints instead of CNVs, then this would be simpler.
  df$end.map <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), df$end.map, c(new.bounds$pos,NA))
  df$end.map.L <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), df$end.map, c(new.bounds$left.bound,NA))
  df$end.map.R <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), df$end.map, c(new.bounds$right.bound,NA))
  df$end.map.win.size <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), NA, c(new.bounds$win.size,NA))
  df$end.map.L.tail <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), NA, c(new.bounds$left.tail,NA))
  df$end.map.R.tail <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), NA, c(new.bounds$right.tail,NA))
  df$end.bin.L <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), NA, c(new.bounds$binL,NA))
  df$end.bin.R <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), NA, c(new.bounds$binR,NA))
  df$end.bin <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), df$end.bin, c(new.bounds$best.bin,NA))
  df$end.binCI.L <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), NA, c(new.bounds$binCI.L,NA))
  df$end.binCI.R <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), NA, c(new.bounds$binCI.R,NA))
  
  df$start.map <- ifelse(c(TRUE,is.na(df$cn[2:nrow(df)])), df$start.map, c(NA, new.bounds$pos))
  df$start.map.L <- ifelse(c(TRUE,is.na(df$cn[2:nrow(df)])), df$start.map, c(NA, new.bounds$left.bound))
  df$start.map.R <- ifelse(c(TRUE,is.na(df$cn[2:nrow(df)])), df$start.map, c(NA, new.bounds$right.bound))
  df$start.map.win.size <- ifelse(c(TRUE,is.na(df$cn[2:nrow(df)])), NA, c(NA, new.bounds$win.size))
  df$start.map.L.tail <- ifelse(c(TRUE,is.na(df$cn[2:nrow(df)])), NA, c(NA, new.bounds$left.tail))
  df$start.map.R.tail <- ifelse(c(TRUE,is.na(df$cn[2:nrow(df)])), NA, c(NA, new.bounds$right.tail))
  df$start.bin.L <- ifelse(c(TRUE,is.na(df$cn[2:nrow(df)])), NA, c(NA, new.bounds$binL))
  df$start.bin.R <- ifelse(c(TRUE,is.na(df$cn[2:nrow(df)])), NA, c(NA, new.bounds$binR))
  df$start.bin <- ifelse(c(TRUE,is.na(df$cn[2:nrow(df)])), df$start.bin, c(NA, new.bounds$best.bin))
  df$start.binCI.L <- ifelse(c(TRUE,is.na(df$cn[2:nrow(df)])), NA, c(NA, new.bounds$binCI.L))
  df$start.binCI.R <- ifelse(c(TRUE,is.na(df$cn[2:nrow(df)])), NA, c(NA, new.bounds$binCI.R))

  # It is possible that small CNVs have breakpoint bounds that have reversed order because the window scans beyond the size of the CNV.
  # Remove those CNVs.
  if (any(df$end.bin < df$start.bin)) {
    message(Sys.time(), ": Removing ",length(which(df$end.bin<df$start.bin)), " CNVs where adjusted end points eliminated a small CNV.")
    df <- filter(df, end.bin >= start.bin)
  }

  # It is possible that end.bin < start.bin at this point because a small CNV is eliminated and the bounds overlap, e.g.
  # Previously: 11111111111111333322222222222
  # Now:        11111111111111111111
  #                              222222222222

  # For now, punt and resolve by just choosing the end position of the first segment and choosing the CI based on the left
  #             11111111111111111111222222222
  # Only do this for non-NAs
  # FIXME: The fix should be to better factor ml.transition2 and apply it twice: once to all segments and a second time to just the overlapping ones.
  # FIXME: We now may, again, have adjacent segments of the same CN, usually CN=2, which should be merged.

  if (nrow(df) > 1) {
    idx <- 1:(nrow(df)-1)

    # row numbers of left CNV of conflicting breakpoints
    trouble <- which(df$end.bin[idx] > df$start.bin[idx+1] & !is.na(df$cn[idx]) & !is.na(df$cn[idx+1]))

    if (length(trouble)>0) {
      message(Sys.time(), sprintf(": FIXME. %s has %s conflicting breakpoints (overlapping flanking CNVs) after MLE adjustments.", sample, length(trouble)))
      if (any(df$cn[trouble]==df$cn[trouble+1])) {
        message(Sys.time(), ": FIXME. And some of these conflicting breakpoints have the same CN and should be merged.")
      }
    }

    # set start of right CNV to end of left CNV
    df$start.map[trouble+1] <- df$end.map[trouble]
    df$start.bin[trouble+1] <- df$end.bin[trouble]
    df$start.map.L[trouble+1] <- df$end.map.L[trouble]
    df$start.map.R[trouble+1] <- df$end.map.R[trouble]
    df$start.map.win.size[trouble+1] <- df$end.map.win.size[trouble]
    df$start.map.L.tail[trouble+1] <- df$end.map.L.tail[trouble]
    df$start.map.R.tail[trouble+1] <- df$end.map.R.tail[trouble]
    df$start.bin.L[trouble+1] <- df$end.bin.L[trouble]
    df$start.bin.R[trouble+1] <- df$end.bin.R[trouble]
    df$start.bin[trouble+1] <- df$end.bin[trouble]
    df$start.binCI.L[trouble+1] <- df$end.binCI.L[trouble]
    df$start.binCI.R[trouble+1] <- df$end.binCI.R[trouble]
  }

  df$cn.2 <- ifelse(is.na(df$cn), 2, df$cn)
  df$dCN.R <- c(df$cn.2[2:nrow(df)]-df$cn.2[1:(nrow(df)-1)],2-df$cn.2[nrow(df)])
  df$dCN.L <- c(2-df$cn.2[1], df$dCN.R[1:nrow(df)-1])
  df$dL <- ifelse(df$dCN.L > 0, 'G', ifelse(df$dCN.L < 0, 'L', 'N'))
  df$dR <- ifelse(df$dCN.R > 0, 'G', ifelse(df$dCN.R < 0, 'L', 'N'))

  df$start.i <- NULL
  df$end.i <- NULL
  
  # FIXME: AND! there are still conflicting segments where the same bp is chosen from different starting points.
  #   22222222222223344888666333222222222
  #   AAAAAAAAAAAAABBCCDDDEEEFFFGGGGGGGGG
  # Results in:
  #   22222222222223344888666333222222222
  #                B
  #                C
  #                DD
  # I think it's possible for new breakpoints to conflict and that overlapping
  # is not predictable. So sort results to find overlaps and remove them.
  # The CN=NA rows were not adjusted, so they may overlap, although not a problem.
    df <- arrange(df, start.bin, end.bin)
    df$ovlp <- c(df$end.bin[idx] > df$start.bin[idx+1] & !is.na(df$cn[idx]) & !is.na(df$cn[idx+1]),FALSE)
  # Because the MLE is performed per breakpoint, it's possible that a 
  # segment is length zero because a pair of breakpoints both optimized to the same position.
  # If so, then segment is length 0. 
    df$zlen <- df$end.bin - df$start.bin == 0
    fail <- df$ovlp | df$zlen
    if (any(fail)) {
        message(Sys.time(), sprintf(": %s segments either overlapped or had zero length. Removing", sum(fail)))
        df <- df[!fail,]
    }


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
