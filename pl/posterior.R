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
# #7: the number of bins upstream and downstream from the breakpoint region to consider, e.g. 10

library(plyr)
library(dplyr)
library(RPostgreSQL)
library(reshape)

cmd.args <- commandArgs(trailingOnly = TRUE)
# Sys.setenv(PGHOST="localhost",PGUSER="dkulp",PGDATABASE="seq", PGOPTIONS="--search_path=data_sfari_batch1c_27apr2017")
# cmd.args <- unlist(strsplit('/home/unix/dkulp/data/out/27Apr2017/data_sfari_batch1C_27Apr2017/B12.L5.Q13.W10.PB0.7.ML1e7/sites_cnv_segs.txt smlcsm data_sfari_batch1C_27Apr2017 data_sfari_batch1C_27Apr2017 0.7 10',' '))
# cmd.args <- unlist(strsplit('/cygwin64/home/dkulp/data/SFARI.27April2017mod/dataC/sites_cnv_segs.txt smlcsm data_sfari_batch1C_27Apr2017 data_sfari_batch1C_27Apr2017 0.7 10',' '))
# cmd.args <- c('/cygwin64/home/dkulp/data/SFARI.27April2017/dataC/sites_cnv_segs.txt','smlcsm','data_sfari_batch1C_27Apr2017','data_sfari_batch1C_27Apr2017','.7', '10')
cnv.seg.fn <- cmd.args[1]
cnv.seg.method <- cmd.args[2]
test.label <- cmd.args[3]
external.label <- cmd.args[4]
external.blend <- as.numeric(cmd.args[5])
PAD <- as.numeric(cmd.args[6])

# predicted CNVs
load(sprintf("%s.%s.Rdata",cnv.seg.fn,cnv.seg.method)) # => cn.segs.merged
cnvs <- as.tbl(cn.segs.merged)
cnvs$idx <- 1:nrow(cnvs)

# connect to DB
db <- src_postgres()
prior.region <- tbl(db,'prior_region')
prior <- tbl(db,'prior')

invisible(dbGetQuery(db$con, "BEGIN TRANSACTION"))


# retrieve them all into memory so it's easy to do a bin=>pos mapping
profile.segments <- tbl(db, 'profile_segment') %>% collect(n=Inf)
# save(profile.segments, file="/cygwin64/tmp/profile_segments.Rdata")
# load("/cygwin64/tmp/profile_segments.Rdata")

# return the bin for the genomic coordinate
posToBin <- function(chr, pos) {
    stopifnot(length(pos)==1)
    bin.df <- filter(profile.segments, chrom==chr & start_pos <= pos & end_pos >= pos)
    stopifnot(nrow(bin.df)==1)
    return(bin.df$bin)
}

# returns genomic position. bin can be a vector of integers.
binToPos <- function(bin) {
    df <- right_join(profile.segments, tibble(bin=bin), by='bin')
    df$start_pos + (df$end_pos - df$start_pos) %/% 2
}

# return a null-op / no change posterior
nc <- function(bin,change) {
  return(data.frame(best=bin, conf.L=bin, conf.R=bin, bin=bin, change=change, prior.int.id=NA_integer_, prior.ext.id=NA_integer_))
}

fetch.prior <- function(label, chr, bin, change) {
  # there may be multiple loss (gain) priors that overlap our bin. Choose only one that best straddles the bin
  pr <- filter(prior.region, label==label & chr==chr & binl <= bin & binr >= bin & dcn==change) %>% collect
  if (nrow(pr) == 0) {
    message(Sys.time(),sprintf(": Warning. No prior at (%s,%s,%s,%s)", label,chr,bin,change))
    return(data.frame())
  } 
  
  if (nrow(pr)>1) {
    message(Sys.time(),sprintf(": Notice. Multiple priors at (%s,%s,%s,%s) id=%s", label,chr,bin,change, pr$id))
  }
  
  pr$dist <- pmin(bin-pr$binl,pr$binr-bin)
  pr.best <- which.max(pr$dist)
  pr.id <- pr$id[pr.best]
  pr.n <- pr$n[pr.best]
  pr.total <- pr$total[pr.best]

  filter(prior, region_id==pr.id) %>% arrange(bin) %>% collect %>% mutate(n=pr.n, total=pr.total)
}

# calculate the CI by growing greadily away from max. 
conf.int <- function(p, pos=seq(1,length(p)), conf=0.95) {
  if (any(is.na(p))) { message(Sys.time(),": Warning interval contains NA. Ignoring for now. FIX ME.") }
  p <- p / sum(p, na.rm=TRUE)  # make a density
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
mk.posterior <- function(df, bin, change) {
  
  if (is.na(change) || change=='N') {
    return(nc(bin,change))
  }
  
                                        # load likelihoods for this sample
  binL <- bin - 2*PAD
  binR <- bin + 2*PAD
  bkpts <- dbGetQuery(db$con, sprintf("SELECT b.chr, b.bkpt_bin as bin, b.sample, loss_ll, gain_ll, no_bkpt_ll, b.label FROM bkpt b WHERE b.chr='%s' AND b.sample IN ('%s') AND b.bkpt_bin BETWEEN %s AND %s AND b.label = '%s' ORDER BY b.chr, b.bkpt_bin", 
                                      df$chr, df$.id, binL, binR, test.label))

  if (nrow(bkpts) == 0) {
    message(Sys.time(),": bkpts returned 0 rows between bins ",binL,"..",binR)
    return(nc(bin,change))
  }

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
  prior.ext <- fetch.prior(external.label, df$chr, bin, change)
  
  # load "prior" from test data 
  # TODO: compute prior on-the-fly, excluding current sample?
  prior.int <- fetch.prior(test.label, df$chr, bin, change)

  if (nrow(prior.ext)==0 && nrow(prior.int)==0) {
    return(nc(bin,change))
  }
  if (first(prior.int$n)==1 && first(prior.ext$n)==1) {
    message(Sys.time(),": Skipping posterior calculation. Prior sample count is 1.")
    return(nc(bin,change))
  }
  
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
      mutate(bkpts, bayes_loss=loss, bayes_gain=gain, bayes_nc=nc, loss.u=NA_real_, gain.u=NA_real_, nc.u=NA_real_)
    }

    
  dbWriteTable(db$con, "posterior_dist", mutate(bkpt.posterior[,c('bin','label','sample','chr',
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
  metric.res <- as.data.frame(conf.int(metric.vals, bkpt.posterior$bin))
  
  # mysterious R voodoo. Crashes with 'Unsupported vector type language' unless I
  # first assign these values to temporary variables.
  pi.id <- ifelse(is.null(prior.int$id), NA_integer_, first(prior.int$id))
  pe.id <- ifelse(is.null(prior.ext$id), NA_integer_, first(prior.ext$id))

  # if priors used, then augment with prior IDs
  return(mutate(metric.res, bin=bin, change=change, 
                prior.int.id=pi.id,
                prior.ext.id=pe.id))
  
}

if (dbExistsTable(db$con, "posterior")) {
#  invisible(dbGetQuery(db$con, "DELETE FROM posterior WHERE label=$1", test.label))
  invisible(dbGetQuery(db$con, "DROP TABLE posterior"))
}
if (dbExistsTable(db$con, "posterior_dist")) {
#  invisible(dbGetQuery(db$con, "DELETE FROM posterior_dist WHERE label=$1", test.label))
  invisible(dbGetQuery(db$con, "DROP TABLE posterior_dist"))
}

# for each predicted CNV, compute a new normalized density for each breakpoint based on the joint probability of the likelihood and prior.
message(Sys.time(),sprintf(": Estimating maximum aposterior breakpoint for %s cnvs", nrow(cnvs)))
res <-
  ddply(cnvs, .(label), function(df) {
    cat(df$label,"\n")
    if (nrow(df) > 1) {
      message(Sys.time(),": BUG: FIx Me. Should only be one row per label from staircase.R")
      print(df) 
    }
    else {
      posterior.L <- mutate(mk.posterior(df, posToBin(df$chr, df$start.map), df$dL), side='L')

      posterior.R <- mutate(mk.posterior(df, posToBin(df$chr, df$end.map), df$dR), side='R')
      
      return(cbind(rbind(posterior.L, posterior.R), data.frame(.id=df$.id, chr=df$chr, label=df$label, idx=df$idx)))
    }
  })

dbWriteTable(db$con, "posterior", mutate(res[,c('best','conf.L','conf.R','bin','change','prior.int.id','prior.ext.id','side','.id','chr','label','idx')], 
                                         seg=label, label=test.label), append=TRUE, row.names = FALSE)


dbCommit(db$con)
dbGetQuery(db$con, "VACUUM ANALYZE posterior")

# write a reduced version of smlcsm to the database
if (dbExistsTable(db$con, "cnv_mle")) { dbGetQuery(db$con, "DROP TABLE cnv_mle CASCADE") }
dbWriteTable(db$con, "cnv_mle", as.data.frame(cnvs[,c(".id","label","cn","chr","start.map","end.map","dCN.L","dCN.R","dL","dR","start.map.L","start.map.R","start.map.win.size","start.map.L.tail","start.map.R.tail","start.bin.L","start.bin.R","start.bin","start.binCI.L","start.binCI.R","end.map.L","end.map.R","end.map.win.size","end.map.L.tail","end.map.R.tail","end.bin.L","end.bin.R","end.bin","end.binCI.L","end.binCI.R","idx")]))

# retrieve a new prediction set, replacing the breakpoints with those estimated here.
cnvs.post <- dbGetQuery(db$con, sprintf('SELECT c.".id", c.label, c.chr, p."conf.L" as "start.CI.L", p.best as "start.map", p."conf.R" as "start.CI.R", p2."conf.L" as "end.CI.L", p2.best as "end.map", p2."conf.R" as "end.CI.R", c.cn, c.idx FROM cnv_mle c, posterior p, posterior p2 WHERE c.label=p.seg AND c.label=p2.seg AND p.side=\'L\' AND p2.side=\'R\' AND p.label=\'%s\' AND p2.label=\'%s\' ORDER by c.idx', test.label, test.label))

# convert all coordinates to genomic for compatibility with other "cn.segs.merged" tables
cn.segs.merged <- mutate(cnvs.post,
                         start.binCI.L=start.CI.L, start.bin=start.map, start.binCI.R=start.CI.R,
                         end.binCI.L=end.CI.L, end.bin=end.map, end.binCI.R=end.CI.R, 
                         start.CI.L=binToPos(start.CI.L), start.map=binToPos(start.map), start.CI.R=binToPos(start.CI.R),
                         end.CI.L=binToPos(end.CI.L), end.map=binToPos(end.map), end.CI.R=binToPos(end.CI.R))

# save
save(cn.segs.merged, file=sprintf("%s.bayescsm.Rdata",cnv.seg.fn))
cn.segs.merged$copy.number <- addNA(as.factor(cn.segs.merged$cn))
write.table(select(cn.segs.merged, .id, label, chr, start.CI.L, as.integer(start.map), start.CI.R, end.CI.L, as.integer(end.map), end.CI.R, copy.number), file=sprintf("%s.bayesCI.tbl",cnv.seg.fn), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# write new posterior version to db
if (dbExistsTable(db$con, "cnv_post")) { invisible(dbGetQuery(db$con, "DROP TABLE cnv_post CASCADE")) }
dbWriteTable(db$con, "cnv_post", cn.segs.merged)
invisible(dbGetQuery(db$con, "CREATE INDEX on cnv_post(chr, \"start.map\", \".id\")"))

rm(db)

