# posterior.R - for each predicted CNV, update its bounds using a bayesian strategy of the
#           likelihood of the test sample and a blended prior generated from the external set and the rest of the test set. 
#
#
# Author: David Kulp, dkulp@broadinstitute.org
#
# Usage: Called using rscript:
#
# #1: filename base - e.g. sites_cnv_segs.txt
# #2: method - which set of flattened segments to read? ('smlcsm'). File read is base.method.Rdata
# #3: db connection - user:host:port:dbname - profile data is read from this connection
# #4: test label e.g. "gpc_wave2_batch2"
# #5: external label e.g. "gpc_wave2_batch1"
# #6: external blend fraction, e.g. ".7" means that the external prior is weighted 70% and the test set prior is weighted 30%.
# #7: the number of bases upstream and downstream from the breakpoint region to consider, e.g. 1000

library(plyr)
library(dplyr)
library(RPostgreSQL)
library(reshape)

cmd.args <- commandArgs(trailingOnly = TRUE)
# cmd.args <- c('C:\\cygwin64\\home\\dkulp\\data\\out\\cnv_seg.B12.L500.Q13.4\\sites_cnv_segs.txt','smlcsm','dkulp:localhost:5432:seq','gpc_wave2_batch1','gpc_wave2_batch1','.7', '1000')
cnv.seg.fn <- cmd.args[1]
cnv.seg.method <- cmd.args[2]
db.conn.str <- cmd.args[3]
test.label <- cmd.args[4]
external.label <- cmd.args[5]
external.blend <- as.numeric(cmd.args[6])
PAD <- as.numeric(cmd.args[7])

# predicted CNVs
load(sprintf("%s.%s.Rdata",cnv.seg.fn,cnv.seg.method)) # => cn.segs.merged
cnvs <- as.tbl(cn.segs.merged)

# connect to DB
db.conn.params <- as.list(unlist(strsplit(db.conn.str,":")))
names(db.conn.params) <- c('user','host','port','dbname')
db <- do.call(src_postgres, db.conn.params)

dbGetQuery(db$con, "BEGIN TRANSACTION")

fetch.prior <- function(label, chr, pos, change) {
  dbGetQuery(db$con, sprintf("SELECT p.*, pr.* FROM prior p, prior_region pr, profile_segment ps 
                             WHERE pr.label='%s' AND ps.chrom='%s' AND ps.start_pos <= %s AND ps.end_pos >= %s AND 
                             pr.chr = ps.chrom AND pr.binL <= ps.bin AND pr.binR >= ps.bin AND p.region_id = pr.id AND pr.dcn ='%s'
                             ORDER BY p.start_pos",
                             label, chr, pos, pos, change))
}

# calculate the CI by growing greadily away from max. 
conf.int <- function(p, pos=seq(1,length(p)), conf=0.95) {
  p <- p / sum(p)  # make a density
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
  return(list(best=pos[best.pos], conf.L=pos[i+1], conf.R=pos[j-1]))
}

# retrieve likelihood and overlapping external & internal priors, if any.
# Blend priors.
# Merge likelihood and priors, dealing with missing data in overlap.
# Compute posterior and return CI.
mk.posterior <- function(df, pos, change) {
  
  if (change=='N') {
    return(data.frame(best=pos, conf.L=pos, conf.R=pos, pos=pos, change=change, prior.int.id=NA_integer_, prior.ext.id=NA_integer_))
  }
  
  # load likelihoods for this sample
  bkpts <- dbGetQuery(db$con, sprintf("SELECT b.*, ps.bin, ps.start_pos, ps.end_pos FROM bkpt b, profile_segment ps 
                                      WHERE b.sample='%s' AND ps.chrom='%s' AND ps.start_pos > %s AND ps.start_pos < %s AND ps.end_pos > %s AND ps.end_pos < %s AND b.chr = ps.chrom AND b.bkpt_bin = ps.bin 
                                      AND b.label='%s' ORDER BY ps.chrom, ps.start_pos", df$.id, df$chr, pos-2*PAD, pos+PAD, pos-PAD, pos+2*PAD, test.label))
  bkpts <- mutate(bkpts,
                  loss=10^-loss_ll, 
                  gain=10^-gain_ll,
                  nc=10^-no_bkpt_ll)
  # bkpts <- mutate(bkpts,
  #                 loss=10^-loss_ll, 
  #                 gain=10^-gain_ll, 
  #                 any=10^-any_ll, 
  #                 no_bkpt=10^-no_bkpt_ll,
  #                 nc=no_bkpt,
  #                 lossZ=loss/(loss+gain+no_bkpt),
  #                 gainZ=gain/(loss+gain+no_bkpt),
  #                 no_bkptZ=no_bkpt/(loss+gain+no_bkpt),
  #                 LlossZ=-log10(lossZ),
  #                 LgainZ=-log10(gainZ),
  #                 Lno_loss=-log10(1-loss),
  #                 Lno_gain=-log10(1-gain),
  #                 Lno_lossZ=-log10(1-lossZ),
  #                 Lno_gainZ=-log10(1-gainZ))
  
  
  # load prior that overlaps the initial breakpoint, if any, from external
  prior.ext <- fetch.prior(external.label, df$chr, pos, change)
  
  # load "prior" from test data 
  # TODO: compute prior on-the-fly, excluding current sample?
  prior.int <- fetch.prior(test.label, df$chr, pos, change)
  
  if (nrow(prior.ext) == 0) {
    priors <- prior.int
  } else if (nrow(prior.int) == 0) {
    priors <- prior.ext
  } else {
    # combine internal and external into a single table on position
    priors <- merge(prior.ext, prior.int, by='bin',all=TRUE)
    
    # blend the stats of interest together
    # if the value is missing for one, then just use the value from the other.
    blend.field.names <- c('loss.u','gain.u','nc.u')
    freq.ext <- first(prior.ext$n / prior.ext$total)
    freq.int <- first(prior.int$n / prior.int$total)
    
    lapply(blend.field.names, function(name) {
      name.x <- paste0(name,'.x')
      name.y <- paste0(name,'.y')
      priors[[name]] <<- external.blend * freq.ext * priors[[name.x]] + (1-external.blend) * freq.int * priors[[name.y]]
      priors[[name]] <<- ifelse(is.na(priors[[name]]), na.omit(priors[[name.x]],priors[[name.y]]), priors[[name]])
    })
  }
  
  
  bkpt.posterior <- 
    if (nrow(priors) > 0) {
      # merge bkpts and priors
      bkpt.mrg <- merge(bkpts, priors, by='bin', all.x=TRUE)
      
      # For missing data, duplicate first or last prior row
      priors.names <- setdiff(names(priors),'bin')
      priors.Ltail <- priors[1,]
      priors.Rtail <- priors[nrow(priors),]
      bkpt.mrg[bkpt.mrg$bin < priors.Ltail$bin, priors.names] <- priors.Ltail[,priors.names]
      bkpt.mrg[bkpt.mrg$bin > priors.Rtail$bin, priors.names] <- priors.Rtail[,priors.names]
      
      # compute product of likelihood * prior
      mutate(bkpt.mrg,
             bayes_loss=loss*loss.u,
             bayes_gain=gain*gain.u,
             bayes_nc=nc*nc.u)
      
    } else {
      # no prior (flat)
      mutate(bkpts, bayes_loss=loss, bayes_gain=gain, bayes_nc=nc, bin=bkpt_bin, loss.u=NA, gain.u=NA, nc.u=NA)
    }
  
  dbWriteTable(db$con, "posterior_dist", mutate(bkpt.posterior[,c('bin','label','sample','chr','start_pos','end_pos',
                                                                  'loss','gain','nc','loss.u','gain.u','nc.u',
                                                                  'bayes_loss','bayes_gain','bayes_nc')], 
                                                seg=label), append=TRUE, row.names = FALSE)
  # 

  # which metric to use? gain or loss
  if (change=='G') {
    metric <- 'bayes_gain'
  } else {
    metric <- 'bayes_loss'
  } 
  
  metric.vals <- bkpt.posterior[,metric]
  
  # find max and compute CI.
  metric.res <- as.data.frame(conf.int(metric.vals, bkpt.posterior$start_pos))
  
  # if priors used, then augment with prior IDs
  return(mutate(metric.res, pos=pos, change=change, 
                prior.int.id=ifelse(is.null(prior.int$id), NA_integer_, first(prior.int$id)),
                prior.ext.id=ifelse(is.null(prior.ext$id), NA_integer_, first(prior.ext$id))))
  
}

if (dbExistsTable(db$con, "posterior")) {
  dbGetQuery(db$con, "DELETE FROM posterior WHERE label=$1", test.label)
}
if (dbExistsTable(db$con, "posterior_dist")) {
  dbGetQuery(db$con, "DELETE FROM posterior_dist WHERE label=$1", test.label)
}

# for each predicted CNV, compute a new normalized density for each breakpoint based on the joint probability of the likelihood and prior.
#pv <- profvis({
  res <-
    ddply(filter(cnvs, cn!=2), .(label), function(df) {
      cat(df$label,"\n")
      if (nrow(df) > 1) {
        print("BUG: FIx Me. Should only be one row per label from staircase.R")
        print(df) 
      }
      else {
        posterior.L <- mutate(mk.posterior(df, df$start.map, df$dL), side='L')
        posterior.R <- mutate(mk.posterior(df, df$end.map, df$dR), side='R')
        return(cbind(rbind(posterior.L, posterior.R), data.frame(.id=df$.id, chr=df$chr, label=df$label)))
      }
    })
#})

dbWriteTable(db$con, "posterior", mutate(res[,c('best','conf.L','conf.R','pos','change','prior.int.id','prior.ext.id','side','.id','chr','label')], 
                                         seg=label, label=test.label), append=TRUE, row.names = FALSE)


dbCommit(db$con)
dbSendQuery(db$con, "VACUUM ANALYZE posterior")

# write a reduced version of smlcsm to the database
if (dbExistsTable(db$con, "cnvs")) { dbSendQuery(db$con, "DROP TABLE cnvs") }
dbWriteTable(db$con, "cnvs", as.data.frame(cnvs[,c(".id","label","cn","chr","start.map","end.map","dCN.L","dCN.R","dL","dR","start.map.L","start.map.R","start.map.win.size","start.map.L.tail","start.map.R.tail","start.bin.L","start.bin.R","start.best.bin","start.binCI.L","start.binCI.R","end.map.L","end.map.R","end.map.win.size","end.map.L.tail","end.map.R.tail","end.bin.L","end.bin.R","end.best.bin","end.binCI.L","end.binCI.R")]))

# retrieve a new prediction set, replacing the breakpoints with those estimated here.
cnvs.post <- dbGetQuery(db$con, sprintf('SELECT c.".id", c.label, c.chr, p."conf.L" as "start.CI.L", p.best as "start.map", p."conf.R" as "start.CI.R", p2."conf.L" as "end.CI.L", p2.best as "end.map", p2."conf.R" as "end.CI.R", c.cn FROM cnvs c, posterior p, posterior p2 WHERE c.label=p.seg AND c.label=p2.seg AND p.side=\'L\' AND p2.side=\'R\' AND p.label=\'%s\' AND p2.label=\'%s\'', test.label, test.label))
cn.segs.merged <- cnvs.post

# save
save(cn.segs.merged, file=sprintf("%s.bayescsm.Rdata",cnv.seg.fn))
cn.segs.merged$copy.number <- addNA(as.factor(cn.segs.merged$cn))
write.table(select(cn.segs.merged, .id, label, chr, start.CI.L, as.integer(start.map), start.CI.R, end.CI.L, as.integer(end.map), end.CI.R, copy.number), file=sprintf("%s.bayesCI.tbl",cnv.seg.fn), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


dbDisconnect(db$con)

