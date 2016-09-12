# staircase.R - modify a set of CNV predictions where a "staircase" model suggests a two state, not three state transition
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
cmd.args <- c('d:/mccarroll/cnv_seg.12.500/sites_cnv_segs.txt','csm','postgres:localhost:5432:seq','500','1200','1000')
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

bin.map <- dbGetQuery(db$con, "select * from profile_segment")
rownames(bin.map) <- bin.map$bin

win.size.bins <- first(which(cumsum(bin.map$elength)>=region.offset)) # set window size (in bins) to the size of region.offset

half.win <- win.size.bins / 2 # each window is divided into 2 equal sides for cnA and cnB

csm.new <- ddply(csm, .(.id), function(df) {
  # df <- csm[csm$.id=='08C79660',]
  sample <- df$.id[1]
  cat(sample,"\n")
  all.profiles <- dbGetQuery(db$con, sprintf("select * from profile_counts where sample='%s'", sample))
  exp.cn <- llply(c(0.1,1:3), function(cn) {
    rollapply(all.profiles$expected[1:(nrow(all.profiles)-half.win)]*cn, half.win, sum)
  })
  obs.cn <- rollapply(all.profiles$observed[1:(nrow(all.profiles)-half.win)], half.win, sum)
  pois.cn <- llply(1:4, function(cn) {
    dpois(obs.cn, exp.cn[[cn]])
  })
  
  ml.transition2 <- function(sample, pos1, pos2, cnA, cnB) {
    cnA <- ifelse(cnA > 3, 4, cnA+1)  # map CN to index into pois.cn. CN > 3 => CN:=3
    cnB <- ifelse(cnB > 3, 4, cnB+1)
    binL <- min(which(bin.map$end_pos >= pos1-region.offset))
    binR <- max(which(bin.map$start_pos <= pos2+region.offset))
    bin.count <- binR - binL + 1
    pA <- pois.cn[[cnA]][binL:(binL+bin.count-win.size.bins)] 
    pB <- pois.cn[[cnB]][(binL+half.win):(binL+bin.count-half.win)]
    jp <- pA * pB    # joint likelihood
    jp.norm <- jp / sum(jp)  # normalized
    best.pos <- bin.map[binL+which.max(jp.norm)+half.win-1,'start_pos']
    #  cat(sprintf("%s:%.0f => %.0f (%.0f) (%.4f..%.4f)\n", sample, pos, best.pos, best.pos-pos, min(-log(jp.norm)), max(-log(jp.norm))))
    if (length(best.pos)>0) {
      return(best.pos)
    } else {
      return(pos1 + (pos2-pos1)/2)  # all NA or NaN
    }
  }
  
  # Iterate across every boundary (skipping small NAs), where the boundary between the first and second segment is idx=1, etc.
  new.bounds <- sapply(1:(nrow(df)-1), function(idx) { 
    ml.transition2(sample, df$end.map[idx], df$start.map[idx+1], df$cn[idx], df$cn[idx+1]) 
  })
  df$end.map <- ifelse(c(is.na(df$cn[1:nrow(df)-1]),TRUE), df$end.map, c(new.bounds,NA))
  df$start.map <- ifelse(c(TRUE,is.na(df$cn[2:nrow(df)])), df$start.map, c(NA, new.bounds))
  return(df)
})

csm.new$seg=with(csm.new, sprintf("SEG_%s_%s_%s", chr, start.map, end.map))
cn.segs.merged <- csm.new
save(cn.segs.merged, file=sprintf("%s.smlcsm.Rdata",cnv.seg.fn))
write.table(select(cn.segs.merged, .id, seg, chr, start.map, end.map, copy.number), file=sprintf("%s.sml.tbl",cnv.seg.fn), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(cn.segs.merged[!is.na(cn.segs.merged$cn),c('chr','start.map','end.map','label','cn')], file=sprintf("%s.smlmrg.bed",cnv.seg.fn), sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
