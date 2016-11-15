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
# #4: the max size of intermediate regions that was skipped when joining
# #5: the window size that was used for aggregating profiles
# #6: the size of the window when scanning for boundary using MLE method

library(plyr)
library(dplyr)
library(RPostgreSQL)
library(zoo)

cmd.args <- commandArgs(trailingOnly = TRUE)
#cmd.args <- c('C:\\cygwin64\\home\\dkulp\\data\\out\\cnv_seg.B12.L500.Q13.4\\sites_cnv_segs.txt','csm','dkulp:localhost:5432:seq','500','1200','1000')
cnv.seg.fn <- cmd.args[1]
cnv.seg.method <- cmd.args[2]
db.conn.str <- cmd.args[3]
max.join <- as.numeric(cmd.args[4])  # size (in nt) of runs to span if flanking calls are the same, e.g. 500
win.size <- as.numeric(cmd.args[5])
region.offset <- as.numeric(cmd.args[6])


load(sprintf("%s.%s.Rdata",cnv.seg.fn,cnv.seg.method)) # => cn.segs.merged
csm <- as.tbl(cn.segs.merged)
csm$len <- csm$seg <- csm$len2 <- csm$gap <- NULL
csm$len <- csm$end.map - csm$start.map 

# identify "small" NAs and remove them
csm <- subset(csm, !is.na(cn) | len > 4*max.join)

# connect to DB
db.conn.params <- as.list(unlist(strsplit(db.conn.str,":")))
names(db.conn.params) <- c('user','host','port','dbname')
db <- do.call(src_postgres, db.conn.params)

dbGetQuery(db$con, "BEGIN TRANSACTION")

# identify candidate stairstep regions
m.prev.idx <- 1:(nrow(csm)-2)
m.idx <- 2:(nrow(csm)-1)
m.next.idx <- 3:(nrow(csm))

# candidates holds the indices into csm of the middle region of a up or down 1-copy difference staircase
# NOTE: assume csm is already ordered by sample and fragment position
candidates <- m.idx[csm$len[m.idx]<win.size & 
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

bin.map <- dbGetQuery(db$con, "select * from profile_segment")
rownames(bin.map) <- bin.map$bin

win.size.bins <- first(which(cumsum(bin.map$elength)>=region.offset)) # set window size (in bins) to the size of region.offset

half.win <- win.size.bins / 2 # each window is divided into 2 equal sides for cnA and cnB

# for each sample, load all profile data, for each transition compute MLE boundary and confidence intervals
# as a side effect, also write the probability of the data over all possible transitions
csm.new <- ddply(csm, .(.id), function(df) {
#csm.new <- ddply(filter(csm, .id %in% c('08C79660','09C100176')), .(.id), function(df) {
  # df <- csm[csm$.id=='08C79660',]
  # df <- csm[csm$.id=='09C100176',]
  sample <- df$.id[1]
  cat(sample,"\n")
  all.profiles <- dbGetQuery(db$con, sprintf("select * from profile_counts where sample='%s' order by bin", sample))
  exp.cn <- llply(c(0.1,1:3), function(cn) {
    rollapply(all.profiles$expected[1:(nrow(all.profiles)-half.win)]*cn, half.win, sum)
  })
  obs.cn <- rollapply(all.profiles$observed[1:(nrow(all.profiles)-half.win)], half.win, sum)
  pois.cn <- llply(1:4, function(cn) {
    dpois(obs.cn, exp.cn[[cn]])
  })
  
  # remove half.win from left/right of pois.cn's vectors to combine them in the next step
  pois.cnL <- llply(pois.cn, function(dp) dp[1:(length(dp)-half.win)])
  pois.cnR <- llply(pois.cn, function(dp) dp[(1+half.win):length(dp)])
  
  # sum over the probability of each pair of indices in idx, e.g. for gain1.idx:
  # pois.cnL[[1]]*pois.cnR[[2]]+pois.cnL[[2]]*pois.cnR[[3]]+pois.cnL[[3]]*pois.cnL[[4]]
  sumprob <- function(idx) {
    -log10(apply(
      apply(idx, 2, function(CN_AB) {
        pois.cnL[[CN_AB[1]]]*pois.cnR[[CN_AB[2]]]
      }), 1, sum))
  }

  gain.idx <- combn(1:4,2)
  loss.idx <- gain.idx[c(2,1),]
  gain1.idx <- matrix(c(1,2,2,3,3,4), nrow=2)
  loss1.idx <- gain1.idx[c(2,1),]
  change.idx <- cbind(gain.idx, loss.idx)
  nochange.idx <- sapply(1:4, function(i) rep(i,2))
  
  bkpt.gain.ll <- sumprob(gain.idx)
  bkpt.gain1.ll <- sumprob(gain1.idx)
  bkpt.loss.ll <- sumprob(loss.idx)
  bkpt.loss1.ll <- sumprob(loss1.idx)
  bkpt.any.ll <- sumprob(change.idx)
  no.bkpt.ll <- sumprob(nochange.idx)

  bkpt.bin <- bin.map$bin[1:length(no.bkpt.ll)+half.win]
  bkpt.df <- filter(data.frame(sample=sample, chr=df$chr[1], bkpt_bin=bkpt.bin, gain_ll=bkpt.gain.ll, 
                               gain1_ll=bkpt.gain1.ll, loss_ll=bkpt.loss.ll, loss1_ll=bkpt.loss1.ll,
                               any_ll=bkpt.any.ll, no_bkpt_ll=no.bkpt.ll),
                    abs(bkpt.loss.ll)<Inf & abs(bkpt.gain.ll)<Inf & abs(no.bkpt.ll)<Inf)
  
  # remove any results from a previous run and write to DB
  dbGetQuery(db$con, sprintf("DELETE FROM bkpt WHERE sample='%s'",sample))
  dbWriteTable(db$con, "bkpt", bkpt.df, append=TRUE, row.names = FALSE)

  ml.transition2 <- function(sample, pos1, pos2, cnA, cnB) {
    cnA <- ifelse(cnA > 3, 4, cnA+1)  # map CN to index into pois.cn. CN > 3 => CN:=3
    cnB <- ifelse(cnB > 3, 4, cnB+1)
    posL <- pos1-region.offset
    posR <- pos2+region.offset
    binL <- min(which(bin.map$end_pos >= posL))
    binR <- max(which(bin.map$start_pos <= posR))
    bin.count <- binR - binL + 1
    pA <- pois.cn[[cnA]][binL:(binL+bin.count-win.size.bins)] 
    pB <- pois.cn[[cnB]][(binL+half.win):(binL+bin.count-half.win)]
    jp <- pA * pB    # joint likelihood
    jp.norm <- jp / sum(jp)  # normalized
    best.jp.bin <- which.max(jp.norm)
#    cat(sprintf("%s:%.0f/%.0f => %.0f (%.4f..%.4f)\n", sample, pos1, pos2, best.pos, min(-log(jp.norm)), max(-log(jp.norm))))
    if (length(best.jp.bin)>0) { # NAs if length is zero
      best.bin <- binL+(best.jp.bin-1)+half.win # first bin passed the transition
      best.pos <- bin.map[best.bin,'start_pos'] # left side of bin is the transition point between bins
      CI <- conf.int(jp.norm)  # returns CI$left and CI$right, which are bin offsets of best.pos
      binCI.L <- best.bin + CI$left
      binCI.R <- best.bin + CI$right
      
      # bounds are generous, including the furthest base from the best.pos in the transition bins.
      left.bound <- bin.map[binCI.L-1,'start_pos']
      right.bound <- bin.map[binCI.R,'end_pos']
      
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
  
  
  return(df)
})

cn.segs.merged <- mutate(csm.new, 
                         seg=sprintf("SEG_%s_%s_%s", chr, start.map, end.map),
                         label=sprintf("%s_%s",seg,.id))

save(cn.segs.merged, file=sprintf("%s.smlcsm.Rdata",cnv.seg.fn))
write.table(select(cn.segs.merged, .id, seg, chr, start.map, end.map, copy.number), file=sprintf("%s.sml.tbl",cnv.seg.fn), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(select(cn.segs.merged, .id, seg, chr, start.map.L, as.integer(start.map), start.map.R, end.map.L, as.integer(end.map), end.map.R, copy.number), file=sprintf("%s.smlCI.tbl",cnv.seg.fn), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(cn.segs.merged[!is.na(cn.segs.merged$cn),c('chr','start.map','end.map','label','cn')], file=sprintf("%s.smlmrg.bed",cnv.seg.fn), sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)

dbCommit(db$con)
dbSendQuery(db$con, "VACUUM ANALYZE")
dbDisconnect(db$con)
