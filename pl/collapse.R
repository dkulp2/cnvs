# collapse.R - create a flattened set of CNV predictions
#
# Author: David Kulp, dkulp@broadinstitute.org
#
# Usage: Called using rscript:
#
# #1: filename base - e.g. sites_cnv_segs.txt
# #2: method - which set of segments to read? ('smlcsm'). File read is base.method.Rdata

library(plyr)
library(dplyr)

cmd.args <- commandArgs(trailingOnly = TRUE)
#cmd.args <- c('/cygwin64/home/dkulp/data/SFARI.27April2017/dataD/sites_cnv_segs.txt','smlcsm','smlxcsm','flt','smlx2csm')
#cmd.args <- c('C:\\cygwin64\\home\\dkulp\\data\\out\\cnv_seg.B12.L500.Q13.4\\sites_cnv_segs.txt','bayescsm','bayesxcsm','bayesflt','bayesx2csm')
cnv.seg.fn <- cmd.args[1]
cnv.seg.method <- cmd.args[2]
method.rdata.out <- cmd.args[3]
method.tbl.out <- cmd.args[4]
method.rdata2.out <- cmd.args[5]

load(sprintf("%s.%s.Rdata",cnv.seg.fn,cnv.seg.method)) # => cn.segs.merged
csm <- as.tbl(cn.segs.merged)

# overlap(df) - returns a flattened data.frame across all samples
overlap <- function(df, start.name='start.map', end.name='end.map', seg='seg') {
  if (nrow(df)==1) {
    return(data.frame(i=1,j=1,start.map=df[[start.name]],end.map=df[[end.name]],
                      x.min=df[[start.name]], x.max=df[[start.name]],
                      y.min=df[[end.name]], y.max=df[[end.name]],
                      x.diff=0, y.diff=0, n=1))
  }
  df2 <- df[order(df[[start.name]]),]
  i2n <- 2:nrow(df2)
  
  # perform closure. repeat until no change.
  df2$end.max <- df2[[end.name]]
  again <- TRUE
  while (again) {
    df2$has.ovlp <- c(FALSE, df2[[start.name]][i2n] < df2$end.max[i2n-1])
    df2$end.max2 <- c(df2[[end.name]][1], ifelse(df2$has.ovlp[i2n], pmax(df2$end.max[i2n],df2$end.max[i2n-1]), df2$end.max[i2n]))
    again <- !all(df2$end.max == df2$end.max2)
    df2$end.max <- df2$end.max2
  }
  
  # starts are where there is no overlap
  # ends are the end.max in the row before the starts. skip row 1. Add end.max of row N
  i <- which(!df2$has.ovlp)
  j <- c(which(!df2$has.ovlp[i2n]), nrow(df2))
  starts <- df2[[start.name]][i]
  ends <- df2$end.max[j]
  
  df.merged <- data.frame(i=i, j=j, starts, ends)
  names(df.merged) <- c('i','j', start.name, end.name)
  
  df.merge <- function(x, i, j) {
    mutate(x, i, j, x.min=min(df2[[eval(start.name,parent.frame(2))]][i:j]), 
           x.max=max(df2[[eval(start.name,parent.frame(2))]][i:j]),
           y.min=min(df2[[end.name]][i:j]), y.max=max(df2[[end.name]][i:j]),
           x.diff=x.max-x.min, y.diff=y.max-y.min,
           n=j-i+1,
           seg=paste0(unique(df2[[seg]][i:j]), collapse = ';'),
           samples=paste0(unique(df2[['.id']][i:j]), collapse=';'))
  }
  df.merged <- ddply(df.merged, .(i,j), df.merge)
  
  return(df.merged)
  
}

# collapse by CN
cn.segs.merged <- as.tbl(mutate(ddply(csm, .(cn, chr), overlap), copy.number=addNA(as.factor(cn))))
cn.segs.merged$.id <- 'All'
save(cn.segs.merged, file=sprintf("%s.%s.Rdata",cnv.seg.fn,method.rdata.out))

write.table(cn.segs.merged[,c('.id', 'seg', 'chr', 'start.map', 'start.map', 'start.map', 'end.map', 'end.map', 'end.map', 'copy.number')], file=sprintf("%s.%s.tbl",cnv.seg.fn,method.tbl.out), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# collapse all, merging different CNs into loss, gain
cn.segs.merged <- ddply(filter(csm, cn!=2), .(chr), function(df) {
  rbind(mutate(overlap(filter(df, dL=='L')), side='L', change='L'),
        mutate(overlap(filter(df, dL=='G')), side='L', change='G'),
        mutate(overlap(filter(df, dR=='L')), side='R', change='L'),
        mutate(overlap(filter(df, dR=='G')), side='R', change='G'))
})

save(cn.segs.merged, file=sprintf("%s.%s.Rdata",cnv.seg.fn,method.rdata2.out))




