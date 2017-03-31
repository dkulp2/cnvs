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
# #4: the max size of intermediate regions that was skipped when joining (in bins, e.g. 5)
# #5: the window size that was used for aggregating profiles (in bases) - in staircase, don't remove intermediate CNVs larger than this 
# #6: the size of the window (in bins) when scanning for boundary using MLE method (e.g. 10 bins)
# #7: data set label

library(plyr)
library(dplyr)
library(RPostgreSQL)

cmd.args <- commandArgs(trailingOnly = TRUE)
#cmd.args <- c('C:\\cygwin64\\home\\dkulp\\data\\SFARI.run2\\data_sfari_batch1B\\sites_cnv_segs.txt','csm','dkulp:localhost:5432:seq','5','1200','10','data_sfari_batch1B')
cnv.seg.fn <- cmd.args[1]
cnv.seg.method <- cmd.args[2]
db.conn.str <- cmd.args[3]
max.join <- as.numeric(cmd.args[4])  # size (in bins) of runs to span if flanking calls are the same, e.g. 5
profile.group.size <- as.numeric(cmd.args[5])  # length of window size, e.g. 1200
win.size.bins <- as.numeric(cmd.args[6])  # in bins. e.g. 10 for a ~1000 nt window.
data.label <- cmd.args[7]


load(sprintf("%s.%s.Rdata",cnv.seg.fn,cnv.seg.method)) # => cn.segs.merged
csm <- as.tbl(cn.segs.merged)
csm$len <- csm$seg <- csm$len2 <- csm$gap <- NULL
csm$len <- csm$end.map - csm$start.map 

# identify "small" NAs and remove them
csm <- subset(csm, !is.na(cn) | len > 4*max.join)

# merge rows with same CN now that small NAs are removed
# merge across small (max.join) gaps unless flanking CNs are both 2.
csm <- ddply(csm, .(.id, chr), function(df) {
  # the first row of a pair of adjacent rows with the same CN
  idx <- 1:(nrow(df)-1)
  adj <- df$cn[idx] == df$cn[idx+1]
  while (any(adj, na.rm=TRUE)) {
    adj.idx <- first(which(adj))
    df[adj.idx, 'end.i'] <- df[adj.idx+1, 'end.i']
    df[adj.idx, 'end.map'] <- df[adj.idx+1, 'end.map']
    df[adj.idx, 'label'] <- paste0(df[adj.idx, 'label'], '_', df[adj.idx+1, 'label'])
    df <- df[-c(adj.idx+1),]
    idx <- 1:(nrow(df)-1)
    adj <- df$cn[idx] == df$cn[idx+1]
  } 
  df$len <- df$end.map - df$start.map 
  return(df)
})

# connect to DB
db.conn.params <- as.list(unlist(strsplit(db.conn.str,":")))
names(db.conn.params) <- c('user','host','port','dbname')
db <- do.call(src_postgres, db.conn.params)

# identify candidate stairstep regions
m.prev.idx <- 1:(nrow(csm)-2)
m.idx <- 2:(nrow(csm)-1)
m.next.idx <- 3:(nrow(csm))

# candidates holds the indices into csm of the middle region of a up or down 1-copy difference staircase
# NOTE: assume csm is already ordered by sample and fragment position
candidates <- m.idx[csm$len[m.idx]<profile.group.size & 
                      (csm$cn[m.idx] - csm$cn[m.prev.idx] == csm$cn[m.next.idx] - csm$cn[m.idx]) & 
                      (abs(csm$cn[m.next.idx] - csm$cn[m.prev.idx]) == 2) &
                      csm$.id[m.prev.idx] == csm$.id[m.next.idx]]

# remove the intermediate candidates
csm <- csm[-na.omit(candidates),]



###################################################################################################################

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
  #csm.new <- ddply(filter(csm, .id %in% c('08C79660','09C100176')), .(.id), function(df) {
#  df <- filter(csm, .id %in% c('MH0131634'))
  # df <- csm[csm$.id=='08C79660',]
  # df <- csm[csm$.id=='09C100176',]
  sample <- df$.id[1]
  cat(sample,"\n")
  
  pois.cn.df <- dbGetQuery(db$con, sprintf("SELECT ps.chrom, ps.start_pos, ps.end_pos, pois.* FROM pois, profile_segment ps WHERE pois.label='%s' AND pois.sample='%s' AND pois.bin=ps.bin AND pois.chr=ps.chrom ORDER BY ps.chrom, ps.start_pos", data.label, sample))

  stopifnot(win.size.bins %% 2 == 0) # must be even number of bins
  half.win <- win.size.bins / 2 # each window is divided into 2 equal sides for cnA and cnB

  pois.cnL <- as.list(pois.cn.df[,grep("cnL",names(pois.cn.df))])
  pois.cnR <- as.list(pois.cn.df[,grep("cnR",names(pois.cn.df))])
  
  ml.transition2 <- function(sample, pos1, pos2, cnA, cnB) {
    cnA <- ifelse(cnA > 3, 4, cnA+1)  # map CN to index into pois.cn. CN > 3 => CN:=3
    cnB <- ifelse(cnB > 3, 4, cnB+1)

    binL <- min(which(pois.cn.df$end_pos >= pos1)) - win.size.bins
    binL <- ifelse(binL < 1, 1, binL)
    posL <- pois.cn.df$start_pos[binL]
    
    binR <- max(which(pois.cn.df$start_pos <= pos2)) + win.size.bins
    binR <- ifelse(binR > nrow(pois.cn.df), nrow(pois.cn.df), binR)
    posR <- pois.cn.df$end_pos[binR]

    bin.count <- binR - binL + 1
    pA <- pois.cnL[[cnA]][binL:(binL+bin.count-win.size.bins)] 
    pB <- pois.cnR[[cnB]][binL:(binL+bin.count-win.size.bins)]
    jp <- pA * pB    # joint likelihood
    jp.norm <- jp / sum(jp)  # normalized
    best.jp.bin <- which.max(jp.norm)
    #    cat(sprintf("%s:%.0f/%.0f => %.0f (%.4f..%.4f)\n", sample, pos1, pos2, best.pos, min(-log(jp.norm)), max(-log(jp.norm))))
    if (length(best.jp.bin)>0) { # NAs if length is zero
      best.bin <- binL+(best.jp.bin-1)+half.win # first bin passed the transition
      best.pos <- pois.cn.df[best.bin, 'start_pos']
      CI <- conf.int(jp.norm)  # returns CI$left and CI$right, which are bin offsets of best.pos
      binCI.L <- best.bin + CI$left
      binCI.R <- best.bin + CI$right
      
      # bounds are generous, including the furthest base from the best.pos in the transition bins.
      left.bound <- pois.cn.df[binCI.L,'start_pos']
      right.bound <- pois.cn.df[binCI.R,'end_pos']
      
      return(data.frame(pos=best.pos, left.bound=left.bound, right.bound=right.bound, win.size=posR-posL+1, left.tail=jp.norm[1], right.tail=jp.norm[length(jp.norm)], binL=binL, binR=binR, best.bin=best.bin, binCI.L=binCI.L, binCI.R=binCI.R))
    } else {
      return(data.frame(pos=pos1 + (pos2-pos1)/2, left.bound=NA, right.bound=NA, win.size=posR-posL+1, left.tail=NA, right.tail=NA, best.bin=NA, binL=binL, binR=binR, binCI.L=NA, binCI.R=NA))
    }
  }
  
  # Iterate across every boundary, where the boundary between the first and second segment is idx=1, etc.
  new.bounds <- ldply(1:(nrow(df)-1), function(idx) { 
    ml.transition2(sample, df$end.map[idx], df$start.map[idx+1], df$cn[idx], df$cn[idx+1]) 
  })
  df$end.map <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), df$end.map, c(new.bounds$pos,NA))
  df$end.map.L <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), df$end.map, c(new.bounds$left.bound,NA))
  df$end.map.R <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), df$end.map, c(new.bounds$right.bound,NA))
  df$end.map.win.size <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), NA, c(new.bounds$win.size,NA))
  df$end.map.L.tail <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), NA, c(new.bounds$left.tail,NA))
  df$end.map.R.tail <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), NA, c(new.bounds$right.tail,NA))
  df$end.bin.L <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), NA, c(new.bounds$binL,NA))
  df$end.bin.R <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), NA, c(new.bounds$binR,NA))
  df$end.best.bin <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), NA, c(new.bounds$best.bin,NA))
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
  df$start.best.bin <- ifelse(c(TRUE,is.na(df$cn[2:nrow(df)])), NA, c(NA, new.bounds$best.bin))
  df$start.binCI.L <- ifelse(c(TRUE,is.na(df$cn[2:nrow(df)])), NA, c(NA, new.bounds$binCI.L))
  df$start.binCI.R <- ifelse(c(TRUE,is.na(df$cn[2:nrow(df)])), NA, c(NA, new.bounds$binCI.R))

  df$cn.2 <- ifelse(is.na(df$cn), 2, df$cn)
  df$dCN.R <- c(df$cn.2[2:nrow(df)]-df$cn.2[1:(nrow(df)-1)],2-df$cn.2[nrow(df)])
  df$dCN.L <- c(2-df$cn.2[1], df$dCN.R[1:nrow(df)-1])
  df$dL <- ifelse(df$dCN.L > 0, 'G', ifelse(df$dCN.L < 0, 'L', 'N'))
  df$dR <- ifelse(df$dCN.R > 0, 'G', ifelse(df$dCN.R < 0, 'L', 'N'))
  
  return(df)
})

cn.segs.merged <- mutate(csm.new, 
                         seg=sprintf("SEG_%s_%s_%s", chr, start.map, end.map),
                         label=sprintf("%s_%s",seg,.id))


save(cn.segs.merged, file=sprintf("%s.smlcsm.Rdata",cnv.seg.fn))
write.table(select(cn.segs.merged, .id, seg, chr, start.map, end.map, copy.number), file=sprintf("%s.sml.tbl",cnv.seg.fn), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(select(cn.segs.merged, .id, seg, chr, start.map.L, as.integer(start.map), start.map.R, end.map.L, as.integer(end.map), end.map.R, copy.number), file=sprintf("%s.smlCI.tbl",cnv.seg.fn), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(cn.segs.merged[!is.na(cn.segs.merged$cn),c('chr','start.map','end.map','label','cn')], file=sprintf("%s.smlmrg.bed",cnv.seg.fn), sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)

dbDisconnect(db$con)
