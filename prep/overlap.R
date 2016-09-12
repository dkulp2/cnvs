
overlap <- function(df, start.name='start.map', end.name='end.map', label='cn', cq='cq', seg='seg', paired.reads='paired.reads') {
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
  
  df.merge <- function(x, i, j) {
    mutate(x, i, j, x.min=min(df2[[eval(start.name,parent.frame(2))]][i:j]), 
           x.max=max(df2[[eval(start.name,parent.frame(2))]][i:j]),
           y.min=min(df2[[end.name]][i:j]), y.max=max(df2[[end.name]][i:j]),
           x.diff=x.max-x.min, y.diff=y.max-y.min,
           cn.min=min(df2[[label]][i:j]), cn.max=max(df2[[label]][i:j]), 
           n=j-i+1, one.cn=cn.min==cn.max,
           cq.min=min(df2[[cq]][i:j]), cq.max=max(df2[[cq]][i:j]),
           seg=paste0(df2[[seg]][i:j], collapse = ';'),
           paired.reads.min=min(df2[[paired.reads]][i:j]),
           paired.reads.max=max(df2[[paired.reads]][i:j]))
  }
  df.merged <- ddply(df.merged, .(i,j), df.merge)
  
  return(df.merged)
  
}

gdo.fn <- "d:/mccarroll/cnv_seg.B12.L500.Q13/gs_dels.txt"
gdo.fn <- "d:/mccarroll/gpc_wave2_batch1/gs_dels.genotypes.txt"
gdo <- read.table(gdo.fn, header=TRUE, sep="\t", stringsAsFactors = FALSE)
colnames(gdo) <- c('.id','seg','chr','start.map','end.map','cn','cq','paired.reads')
gdo$evidence <- ' Known'
gdo$copy.number <- addNA(as.factor(gdo$cn))

gs.dels.flt <- ddply(gdo, .(.id), overlap)
gs.dels.flt0 <- ddply(subset(gdo, paired.reads > 0), .(.id), overlap)
gs.dels.flt2 <- overlap(subset(gs.dels.flt, n>1), label='cn.min')
gs.dels.flt3 <- overlap(gdo)

# create a flattened set of DELs with paired.reads > 0, CQ>13, no conflict CNs
gs.dels.flt4 <- ddply(gdo, .(chr), function(gdo.chr) {
  subset(ddply(subset(gdo.chr, paired.reads > 0 & cq >= 13), .(.id), overlap), one.cn)
})
gs.dels.flt4$cn <- gs.dels.flt4$cn.min
gs.dels.flt4$cq <- gs.dels.flt4$cq.min
gs.dels.flt4$paired.reads <- gs.dels.flt4$paired.reads.min
gdo.fn2 <- "d:/mccarroll/gpc_wave2_batch1/gs_dels_flt.genotypes.txt"
write.table(gs.dels.flt4[,c('.id','seg', 'chr', 'start.map', 'end.map', 'cn', 'cq', 'paired.reads')], 
            file=gdo.fn2,
            col.names=c('SAMPLE','ID','CHROM','START','END','GT','GQ', 'READCOUNT'), row.names=FALSE,
            quote=FALSE, sep="\t")

# > ddply(gs.dels.flt, .(n), summarize, cnt=length(n))
# n  cnt
# 1 1 1140
# 2 2  316
# 3 3  214
# 4 4   59
# 5 5   38
# 6 6    8

# > ddply(gs.dels.flt, .(one.cn), summarize, m=sum(n))
# one.cn    m
# 1  FALSE  151
# 2   TRUE 2737
# > nrow(gs.dels.flt)
# [1] 1775
# > nrow(gs.dels.orig)
# [1] 2888
# > nrow(gs.dels.flt)/nrow(gs.dels.orig)
# [1] 0.6146122
# > 151/2737
# [1] 0.05516989
# 
# > summary(gs.dels.flt[gs.dels.flt$x.diff > 0,]$x.diff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.0   148.0   187.0   296.4   245.0  3957.0 

