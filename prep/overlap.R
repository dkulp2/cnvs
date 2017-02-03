# overlap.R - pre-process "knowns" into collapsed sets
#
# Author: dkulp@broadinstitute.org
#
# There are two primary input sets: a genome strip deletion set, to be used for sensitivity analysis, and
# a full CNV pipeline set to be used for specificity analysis.
#
# The inputs to this program are derived by first running gs_del2tab.
#
# The deletion set is filtered here for a high quality subset.
# The full CNV set is intended to be maximal, so no filtering is performed.
# Both tables are then collapsed by sample ("sampleseg") and by site.

library(plyr)
library(dplyr)
library(pryr)

# overlap - takes a data.frame (df) and returns a data.frame with all overlapping regions collapsed into one with the seg
# column concatenated together.
# The return data.frame also includes summary information for the collapsed sets, i.e. min/max CN, min/max CNQ, min/max paired read count,
#   the maximal difference of left and right bounds, whether or not there is more than one CN, and the number of segments merged.
# choose.best==FALSE: collapse all overlapping samplesegs into a maximal region
# choose.best==TRUE: return the sampleseg with the highest paired.read evidence
overlap <- function(df, start.name='start.map', end.name='end.map', label='cn', cq='cq', seg='seg', paired.reads='paired.reads',
                    choose.best=FALSE) {
  if (nrow(df)==1) {
    return(data.frame(i=1,j=1,start.map=df[[start.name]],end.map=df[[end.name]],
                      x.min=df[[start.name]], x.max=df[[start.name]],
                      y.min=df[[end.name]], y.max=df[[end.name]],
                      x.diff=0, y.diff=0, n=1, cq.min=df[[cq]], one.cn=TRUE,
                      cq.max=df[[cq]], seg=df[[seg]],
                      paired.reads.min=df[[paired.reads]], paired.reads.max=df[[paired.reads]],x.n=1,y.n=1
                      ))
  }
  df2 <- df[order(df[[start.name]]),]
  i2n <- 2:nrow(df2)
  
  # perform closure. repeat until no change.
  df2$has.ovlp2 <- rep(FALSE, nrow(df2))
  df2$end.max <- df2[[end.name]]
  again <- TRUE
  while (again) {
    df2$has.ovlp <- c(FALSE, df2[[start.name]][i2n] < df2$end.max[i2n-1])
    df2$end.max <- c(df2[[end.name]][1], ifelse(df2$has.ovlp[i2n], pmax(df2$end.max[i2n],df2$end.max[i2n-1]), df2$end.max[i2n]))
    again <- !all(df2$has.ovlp==df2$has.ovlp2)
    df2$has.ovlp2 <- df2$has.ovlp
  }
  
  # starts are where there is no overlap
  # ends are the end.max in the row before the starts. skip row 1. Add end.max of row N
  i <- which(!df2$has.ovlp)
  j <- c(which(!df2$has.ovlp[i2n]), nrow(df2))
  starts <- df2[[start.name]][i]
  ends <- df2$end.max[j]
  
  df.merged <- data.frame(i=i, j=j, starts, ends)
  names(df.merged) <- c('i','j', start.name, end.name)

  df.merge <- function(x) {
    mutate(x, i, j, 
           x.min=min(df2[[eval(start.name,parent.frame(2))]][i:j]), 
           x.max=max(df2[[eval(start.name,parent.frame(2))]][i:j]),
           y.min=min(df2[[end.name]][i:j]), y.max=max(df2[[end.name]][i:j]),
           x.diff=x.max-x.min, y.diff=y.max-y.min,
           cn.min=min(df2[[label]][i:j]), cn.max=max(df2[[label]][i:j]), 
           n=j-i+1, one.cn=cn.min==cn.max,
           cq.min=min(df2[[cq]][i:j]), cq.max=max(df2[[cq]][i:j]),
           seg=paste0(unique(df2[[seg]][i:j]), collapse = ';'),
           paired.reads.min=min(df2[[paired.reads]][i:j]),
           paired.reads.max=max(df2[[paired.reads]][i:j]),
           x.n=length(unique(df2[[start.name]][i:j])),
           y.n=length(unique(df2[[end.name]][i:j])))
  }
  df.merged <- ddply(df.merged, .(i,j), df.merge)
  
  # if choose.best, then set start.map and end.map to the bounds of the sampleseg with the highest read pair evidence
  if (choose.best) {

    df.merged <- ddply(df.merged, .(i,j), mutate,
                       k=(i:j)[which.max(df2[[paired.reads]][i:j])],
                       start.map=df2[[start.name]][k],
                       end.map=df2[[end.name]][k])

  }
  return(df.merged)
  
}

# flattenFile takes a data.frame and performs sampleseg and site overlap, writing the results to the designated output files.
# restrict: an expression that is applied to the data.frame. (Could be applied before call now, but cool example on how to pass unevaluated R) 
# remove.discordant.cn: if TRUE skips any overlap sets that have different CN values
flattenFile <- function(df, sampleseg.fn, site.fn, restrict=expression(TRUE), remove.discordant.cn=TRUE, choose.best=FALSE) {
  
  # gs.dels.flt <- ddply(gdo, .(.id), overlap)
  # gs.dels.flt0 <- ddply(subset(gdo, paired.reads > 0), .(.id), overlap)
  # gs.dels.flt2 <- overlap(subset(gs.dels.flt, n>1), label='cn.min')
  # gs.dels.flt3 <- overlap(gdo)
  
  gs.dels.flt4 <- ddply(df, .(chr), function(gdo.chr) {
    a <- ddply(subset(gdo.chr, eval(restrict)), .(.id), overlap, choose.best=choose.best)
    if (remove.discordant.cn) { subset(a, one.cn) } else { a }
  })
  
  gs.dels.flt4$cn <- ifelse(gs.dels.flt4$one.cn, gs.dels.flt4$cn.min, 99) # hack: mark collapsed CN column as discordant
  gs.dels.flt4$cq <- gs.dels.flt4$cq.min
  if (choose.best) { 
    gs.dels.flt4$paired.reads <- gs.dels.flt4$paired.reads.max
  } else {
    gs.dels.flt4$paired.reads <- gs.dels.flt4$paired.reads.min
  }
  write.table(gs.dels.flt4[,c('.id','seg', 'chr', 'start.map', 'end.map', 'cn', 'cq', 'paired.reads')], 
              file=sampleseg.fn,
              col.names=c('SAMPLE','ID','CHROM','START','END','GT','GQ', 'READCOUNT'), row.names=FALSE,
              quote=FALSE, sep="\t")
  
  
  # create a flattened set of DELS with same criteria as gs.dels.flt4, but collapsing across samples
  gs.dels.flt5 <- ddply(gs.dels.flt4, .(chr, cn), overlap)
  gs.dels.flt5$cq <- gs.dels.flt5$cq.min
  gs.dels.flt5$paired.reads <- gs.dels.flt5$paired.reads.min
  
  gs.dels.flt5$.id <- 'All'
  write.table(gs.dels.flt5[,c('.id','seg', 'chr', 'start.map', 'end.map', 'cn', 'cq', 'paired.reads')], 
              file=site.fn,
              col.names=c('SAMPLE','ID','CHROM','START','END','GT','GQ', 'READCOUNT'), row.names=FALSE,
              quote=FALSE, sep="\t")
  return(list(sampleseg=as.tbl(gs.dels.flt4), site=as.tbl(gs.dels.flt5)))
}


# read the file generated by gs_del2tab
read_known_file <- function(fn) {
  gdo <- read.table(fn, header=TRUE, sep="\t", stringsAsFactors = FALSE)
  colnames(gdo) <- c('.id','seg','chr','start.map','end.map','cn','cq','paired.reads')
  gdo$evidence <- ' Known'
  gdo$copy.number <- addNA(as.factor(gdo$cn))
  return(as.tbl(gdo))  
} 

gsdel.dir <- dirname(Sys.getenv('gsdelFile'))
gsdel.base <- sub('.vcf.gz','.txt',basename(Sys.getenv('gsdelFile')))

# create a flattened set of DELs with paired.reads > 0, CQ>13, no conflict CNs
#gstrip_dels <- read_known_file("/home/dkulp/data/gpc_wave2_batch1/gs_dels.genotypes.txt")
#gstrip_dels <- read_known_file("C:/cygwin64/home/dkulp/data/gpc_wave2_batch1/gs_dels.genotypes.txt")
gstrip_dels <- read_known_file(sprintf("%s/%s",gsdel.dir,gsdel.base))

ret <- flattenFile(gstrip_dels,
                   sprintf("%s/gs_dels_flt.genotypes.txt", gsdel.dir),
                   sprintf("%s/gs_dels_xflt.genotypes.txt", gsdel.dir),
                   restrict=quote(paired.reads > 0 & cq >= 13), remove.discordant.cn = TRUE)

ret <- flattenFile(gstrip_dels,
                   sprintf("%s/gs_dels_best.genotypes.txt",gsdel.dir),
                   "/dev/null",
                   restrict=quote(paired.reads > 0 & cq >= 13), remove.discordant.cn = TRUE, choose.best=TRUE)

gscnv.dir <- dirname(Sys.getenv('gscnvFile'))
gscnv.base <- sub('.vcf.gz','.txt',basename(Sys.getenv('gscnvFile')))

gstrip_cnv_del <- rbind(gstrip_dels, read_known_file(sprintf("%s/%s", gsdel.dir, gsdel.base)))
ret <- flattenFile(gstrip_cnv_del,
                   sprintf("%s/gs_cnv_del_flt.genotypes.txt", gscnv.dir),
                   sprintf("%s/gs_cnv_del_xflt.genotypes.txt",gscnv.dir),
                   remove.discordant.cn = FALSE)


