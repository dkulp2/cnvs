# priors.R - for each region with a predicted CNV breakpoint, generate a prior from the likelihoods of the 
#           samples with a CNV and save
#
#
# Author: David Kulp, dkulp@broadinstitute.org
#
# Usage: Called using rscript:
#
# #1: filename base - e.g. sites_cnv_segs.txt
# #2: method - which set of flattened segments to read? ('smlxcsm'). File read is base.method.Rdata
# #3: db connection - user:host:port:dbname - profile data is read from this connection
# #4: the label for the data set to read and write, e.g. "gpc_wave2_batch1"
# #5: the number of bin upstream and downstream from the breakpoint region to consider

library(plyr)
library(dplyr)
library(RPostgreSQL)
library(caroline)
library(reshape)

cmd.args <- commandArgs(trailingOnly = TRUE)
# Sys.setenv(PGHOST="localhost",PGUSER="dkulp",PGDATABASE="seq", PGOPTIONS="--search_path=data_sfari_batch1d_27apr2017")
# cmd.args <- unlist(strsplit('/home/unix/dkulp/data/out/27Apr2017/data_sfari_batch1D_27Apr2017/B12.L5.Q13.W10.PB0.7.ML1e7/sites_cnv_segs.txt smlx2csm data_sfari_batch1D_27Apr2017 10',' '))
# cmd.args <- c('/cygwin64/home/dkulp/data/SFARI.27April2017mod/dataD/sites_cnv_segs.txt','smlx2csm','data_sfari_batch1D_27Apr2017','10')
cnv.seg.fn <- cmd.args[1]
cnv.seg.method <- cmd.args[2]
data.label <- cmd.args[3]
PAD <- as.numeric(cmd.args[4])

load(sprintf("%s.%s.Rdata",cnv.seg.fn,cnv.seg.method)) # => cn.segs.merged
cnvx <- filter(as.tbl(cn.segs.merged), x.diff < quantile(cn.segs.merged$x.diff,.95) & y.diff < quantile(cn.segs.merged$y.diff,.95))
cnvx$freq <- unlist(lapply(strsplit(cnvx$samples,';'), function(s) { length(s) }))
cnvx$samples.quoted <- unlist(lapply(strsplit(cnvx$samples,';'), function(s) { paste(s,collapse="','") }))


# connect to DB
db <- src_postgres()

# return the bin for the genomic coordinate
posToBin <- function(chr, pos) {
  bin.df <- dbGetQuery(db$con, sprintf("SELECT ps.bin, ps.start_pos, ps.end_pos FROM profile_segment ps WHERE ps.chrom='%s' AND ps.start_pos <= %s AND ps.end_pos >= %s", chr, pos, pos))
  stopifnot(nrow(bin.df)==1)
  return(bin.df$bin)
}

# get the sample count for the data.label set, so we can store the CNV frequency for each prior
total.samples <- dbGetQuery(db$con, sprintf("SELECT count(distinct sample) FROM pois WHERE label='%s'", data.label))$count

dbGetQuery(db$con, "BEGIN TRANSACTION")

mk.prior <- function(chr, x1, x2, samples, label) {
  bin1 <- posToBin(chr, x1) - PAD
  bin2 <- posToBin(chr, x2) + PAD
  bkpts <- dbGetQuery(db$con, sprintf("SELECT b.sample, bkpt_bin as bin, loss_ll, gain_ll, no_bkpt_ll FROM bkpt b WHERE b.chr='%s' AND b.sample IN ('%s') AND b.bkpt_bin BETWEEN %s AND %s AND b.label = '%s'", 
                                      chr, samples, bin1, bin2, label))
  stopifnot(nrow(bkpts)>0)

  # Retrieve pois likelihoods for each sample and position.
  bkpts <- mutate(bkpts,
                  loss = 10^-loss_ll,
                  gain = 10^-gain_ll,
                  nc = 10^-no_bkpt_ll)
  
  # assume that loss (or gain) definitely occurs somewhere in x1..x2 window for each sample.
  # assume that p(loss_i|x_j) ~ p(x_j|loss_i)p(loss_i) where p(loss_i) is uniform.
  # so normalize the likelihoods per sample to obtain p(loss_i|x_j), i.e. p(loss|x_j)=1.
  bkpts <- ddply(bkpts, .(sample), mutate,
                  lossZ = loss/sum(loss),
                  gainZ = gain/sum(gain),
                  ncZ = nc/sum(nc)
  )

  # assume that each p(loss|x_j) has equal wait. TODO: add call confidence per sample  
  bkpt.prior <- mutate(ddply(bkpts, .(bin), summarize,
                             loss.u=mean(lossZ),
                             gain.u=mean(gainZ),
                             nc.u=mean(ncZ)),
                       Lloss.u=-log(loss.u),
                       Lgain.u=-log(gain.u),
                       Lnc.u=-log(nc.u))
  
}

plot.prior <- function(bkpt.prior) {
  bkpt.prior.m <- melt(bkpt.prior, 'start_pos', c('loss.u','gain.u'))
                  
  ggplot(bkpt.prior.m, 
         aes(x=start_pos, y=value, color=variable))+geom_point()+geom_line() + theme(axis.title.x = element_blank(),plot.margin=unit(c(0,.75,0,.6),'in'),legend.position="bottom") + ylab("Log Prob")  
  
}

save.prior <- function(bkpt.prior, chr, seg, samples, label, n, total.samples, side, dCN) {
  # remove any results from a previous run and write to DB
  # If tables don't exist, then they will be created on the first write
  prior.exists <- dbExistsTable(db$con, "prior")

  if (prior.exists) {
    dbGetQuery(db$con, sprintf("DELETE FROM prior_region WHERE seg='%s' AND label='%s' AND side='%s'",seg, label, side))
  } else {
    dbGetQuery(db$con, "CREATE TABLE prior_region (id serial primary key, chr varchar(2), binL integer, binR integer, seg text, samples text, label text, n integer, total integer, side char, dCN char, FOREIGN KEY (chr, binL) REFERENCES profile_segment(chrom,bin), FOREIGN KEY (chr, binR) REFERENCES profile_segment(chrom,bin))")
  }
  
  this.prior.region <- dbGetQuery(db$con, sprintf("SELECT * from prior_region WHERE chr='%s' and binL='%s' and side='%s' and dCN='%s'", chr, min(bkpt.prior$bin), side, dCN))
  if (nrow(this.prior.region) >= 1) {
    cat(sprintf("There are already %s prior_region rows starting at %s\n", nrow(this.prior.region), min(bkpt.prior$bin)))
  }
  
  region_id <- dbWriteTable2(db$con, "prior_region", data.frame(chr=chr, binL=min(bkpt.prior$bin), binR=max(bkpt.prior$bin), 
                                                                seg=seg, samples=samples, label=label, n=n, total=total.samples,
                                                                side=side, dCN=dCN, stringsAsFactors = FALSE), 
                             append=TRUE, row.names=FALSE)
  
  dbWriteTable(db$con, "prior", mutate(bkpt.prior,region_id=as.integer(region_id)), append=TRUE, row.names = FALSE)

  if (!prior.exists) {  # new table
    dbGetQuery(db$con, "CREATE INDEX ON prior_region(label, seg)")
    dbGetQuery(db$con, "CREATE INDEX ON prior_region(chr, binL)")
    dbGetQuery(db$con, "ALTER TABLE prior ADD FOREIGN KEY(region_id) REFERENCES prior_region(id) ON DELETE CASCADE")
    dbGetQuery(db$con, "CREATE INDEX ON prior(region_id)")
  }
  
}

# drop both tables.
if (dbExistsTable(db$con, "prior_region")) invisible(dbGetQuery(db$con, "DROP TABLE prior_region cascade"))
if (dbExistsTable(db$con, "prior")) invisible(dbGetQuery(db$con, "DROP TABLE prior cascade"))

# for each start.map and stop.map 
#  retrieve the prior data for the samples in the cluster
ddply(cnvx, .(chr,seg,side,change), function(df) {
  if (nrow(df) > 1) {
    cat(sprintf("Bug with multiple samples! Skipping\n"))
  } else {
    if (df$side == 'L') { 
      pos.min <- df$x.min; pos.max <- df$x.max 
    } else { 
      pos.min <- df$y.min; pos.max <- df$y.max 
    }
    bkpt.prior <- mk.prior(df$chr, pos.min, pos.max, df$samples.quoted, data.label)
    cat(sprintf("%s/%s %s:%s-%s ", df$side, df$change, df$chr, df$start.map, df$end.map))
    cat(sprintf("(%s-%s)\n", min(bkpt.prior$bin), max(bkpt.prior$bin)))
    save.prior(bkpt.prior, df$chr, df$seg, df$samples, data.label, df$n, total.samples, df$side, df$change)
  }
})

dbCommit(db$con)
dbSendQuery(db$con, "VACUUM ANALYZE prior")
dbDisconnect(db$con)

