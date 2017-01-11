# likelihoods.R - scan across profiles and store likelihoods of gain, loss, nc at each position in chromosomes for each sample.
# 
# Two tables are generated with rows for each bin and sample, labeled with the data.label:
#   pois - poison likelihoods for a set of bins half the size of win.size to the left and right of each position for each copy number 
#   bkpt - various combinations of pois representing different types of transitions, i.e. gain, loss, no change.
#
#Table "public.pois"
#  Column  |       Type       | Modifiers
# ---------+------------------+-----------
#  label   | text             |
#  sample  | text             |
#  bin     | integer          |
#  cnL_0_1 | double precision |
#  cnL_1   | double precision |
#  cnL_2   | double precision |
#  cnL_3   | double precision |
#  cnR_0_1 | double precision |
#  cnR_1   | double precision |
#  cnR_2   | double precision |
#  cnR_3   | double precision |
#
#
# Table "public.bkpt"
#    Column   |         Type         | Modifiers
# ------------+----------------------+-----------
#  sample     | text                 |
#  chr        | character varying(2) |
#  bkpt_bin   | integer              |
#  gain_ll    | double precision     |
#  gain1_ll   | double precision     |
#  loss_ll    | double precision     |
#  loss1_ll   | double precision     |
#  any_ll     | double precision     |
#  no_bkpt_ll | double precision     |
#
# Author: David Kulp, dkulp@broadinstitute.org
#
# Usage: Called using rscript:
#
# #1: a sites file to identify the min/max regions per chromosome to process
# #2: db connection - user:host:port:dbname - profile data is read from this connection
# #3: the size of the sliding window in which the left half is CN_a and the right is CN_b
# #4: the label for the data set, e.g. "gpc_wave2_batch1"

library(plyr)
library(dplyr)
library(RPostgreSQL)
library(zoo)

cmd.args <- commandArgs(trailingOnly = TRUE)
#cmd.args <- c('C:\\cygwin64\\home\\dkulp\\data\\out\\cnv_seg.B12.L500.Q13.4\\sites_cnv_segs.txt','dkulp:localhost:5432:seq','1000','gpc_wave2_batch1')
cnv.seg.fn <- cmd.args[1]
db.conn.str <- cmd.args[2]
win.size <- as.numeric(cmd.args[3])
data.label <- cmd.args[4]

# connect to DB
db.conn.params <- as.list(unlist(strsplit(db.conn.str,":")))
names(db.conn.params) <- c('user','host','port','dbname')
db <- do.call(src_postgres, db.conn.params)

dbGetQuery(db$con, "BEGIN TRANSACTION")

bin.map <- dbGetQuery(db$con, "select * from profile_segment")
rownames(bin.map) <- bin.map$bin

win.size.bins <- first(which(cumsum(bin.map$elength)>=win.size)) # set window size (in bins) to the size of win.size

half.win <- win.size.bins / 2 # each window is divided into 2 equal sides for cnA and cnB

cat("Counting samples in profile_counts")
samples <- dbGetQuery(db$con, "select distinct sample from profile_counts")
#samples <- data.frame(sample=c('08C79660','09C100176'), stringsAsFactors = FALSE)
cat(sprintf("Processing %s samples\n", nrow(samples)))

sapply(samples$sample, function(sample) {
  cat(sample,"\n")
  
  # FIXME: operate per chromosome
  all.profiles <- dbGetQuery(db$con, sprintf("select * from profile_counts where sample='%s' order by chrom, bin", sample))
  
  cn.vals <- c(0.1,1:3)

  exp.cn <- llply(cn.vals, function(cn) {
    rollapply(all.profiles$expected*cn, half.win, sum)
  })
  obs.cn <- rollapply(all.profiles$observed, half.win, sum)
  pois.cn <- llply(1:length(cn.vals), function(cn) {
    dpois(obs.cn, exp.cn[[cn]])
  })
  
  # remove half.win from left/right of pois.cn's vectors to combine them in the next step
  pois.cnL <- llply(pois.cn, function(dp) dp[1:(length(dp)-half.win)])
  pois.cnR <- llply(pois.cn, function(dp) dp[(1+half.win):length(dp)])
  
  # add labels 
  names(pois.cnL) <- gsub('\\.','_',sprintf("cnL_%s", cn.vals))
  names(pois.cnR) <- gsub('\\.','_',sprintf("cnR_%s", cn.vals))
  
  # flatten pois.cnL and pois.cnR into a single data.frame
  # bin is the leftmost bin of cnL
  pois.cn.df <- cbind(mutate(data.frame(label=data.label, sample=sample, 
                                        bin=all.profiles$bin[1:(nrow(all.profiles)-(half.win-1)-half.win)]),
                             chr=bin.map$chrom[bin]), 
                      do.call(data.frame, pois.cnL), do.call(data.frame, pois.cnR))
  
  # remove any results from a previous run and write to DB
  # If table doesn't exist, then it will be created on the first write
  pois.exists <- dbExistsTable(db$con, "pois")
  if (pois.exists) {
    dbGetQuery(db$con, sprintf("DELETE FROM pois WHERE sample='%s'",sample))
  }
  dbWriteTable(db$con, "pois", pois.cn.df, append=TRUE, row.names = FALSE)
  if (!pois.exists) {  # new table
    dbGetQuery(db$con, "CREATE UNIQUE INDEX ON pois(label,sample,chr,bin)")
    dbGetQuery(db$con, "ALTER TABLE pois ADD FOREIGN KEY(chr,bin) REFERENCES profile_segment(chrom,bin)")
  }
  
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
  bkpt.df <- filter(data.frame(label=data.label, sample=sample, chr=all.profiles$chrom[1], bkpt_bin=bkpt.bin, gain_ll=bkpt.gain.ll, 
                               gain1_ll=bkpt.gain1.ll, loss_ll=bkpt.loss.ll, loss1_ll=bkpt.loss1.ll,
                               any_ll=bkpt.any.ll, no_bkpt_ll=no.bkpt.ll),
                    abs(bkpt.loss.ll)<Inf & abs(bkpt.gain.ll)<Inf & abs(no.bkpt.ll)<Inf)
  
  # remove any results from a previous run and write to DB
  dbGetQuery(db$con, sprintf("DELETE FROM bkpt WHERE sample='%s'",sample))
  dbWriteTable(db$con, "bkpt", bkpt.df, append=TRUE, row.names = FALSE)
  
})

dbCommit(db$con)
dbSendQuery(db$con, "VACUUM ANALYZE bkpt")
dbSendQuery(db$con, "VACUUM ANALYZE pois")
dbDisconnect(db$con)
