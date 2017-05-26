# concordance.R
#
# compare the predictions among an IBD quartet

library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(sqldf)
library(grid)
library(magrittr)
library(zoo)
library(RPostgreSQL)

load("/Cygwin64/home/dkulp/data/SFARI.19May2017/segments.Rdata")
load("/Cygwin64/home/dkulp/data/SFARI.15May2017/segments.Rdata")

MIN.OVLP.LEN <- 1400
sapply(segs.group, function(df) {
  # print(ggplot(filter(df, !is.na(ibd.state) & !(cnv %in% c('WT','DEL','BI')) & len>MIN.OVLP.LEN), aes(x=start, xend=end, y=fam, yend=fam, color=concordant)) +
  #         geom_segment(size=2) + geom_point() + facet_grid(.~cnv) + ggtitle(paste0('Samples ',df$family[1],'..',df$family[nrow(df)],"\ninterval > ",MIN.OVLP.LEN)))
  print(ggplot(filter(df, !is.na(ibd.state) & !all.wt & len>MIN.OVLP.LEN & cnv!='WT'), aes(x=start, xend=end, y=fam, yend=fam, color=concordant)) +
     geom_segment(size=2) + geom_point() + facet_grid(ibd.state~cnv) + ggtitle(paste0('Samples ',df$family[1],'..',df$family[nrow(df)],"\ninterval > ",MIN.OVLP.LEN)))
})

sa1 <- filter(segs.all[[1]], !is.na(ibd.state) & cnv=='DEL' & !cn.na)
print(ggplot(filter(sa1, len>MIN.OVLP.LEN & !concordant), aes(x=old.start, xend=old.end, y=fam, yend=fam, color=len)) + scale_colour_gradient(high = "red", low = "orange") +
        geom_segment(size=2) + geom_point() + ggtitle(sprintf("Discordant Sites\nInterval > %s nt",MIN.OVLP.LEN)))

sa2 <- filter(segs.all[[1]], !is.na(ibd.state))
print(ggplot(filter(sa2, len>MIN.OVLP.LEN), aes(x=old.start, xend=old.end, y=fam, yend=fam, color=len)) + scale_colour_gradient(high = "red", low = "orange") +
        geom_segment(size=2) + geom_point() + facet_grid(concordant~.)+ ggtitle(sprintf("Sites By Concordancy (T,F,NA)\nInterval > %s nt",MIN.OVLP.LEN)))

print(ggplot(filter(sa2, len>1e5 & concordant), aes(x=old.start, xend=old.end, y=fam, yend=fam, color=len)) + scale_colour_gradient(high = "red", low = "orange") +
        geom_segment(size=2) + geom_point() + facet_grid(concordant~.)+ ggtitle(sprintf("Sites By Concordancy (T,F,NA)\nInterval > %s nt",MIN.OVLP.LEN)))

# xlm <- xlim(0.1584e7, 0.1591e7)
# xlm <- xlim(2.949e7, 2.951e7)
# xlm <- xlim(4.1243e7, 4.125e7)
# xlm <- xlim(3.3242e7, 3.3245e7)
# 
# print(ggplot(filter(sa1, len>MIN.OVLP.LEN & !concordant), aes(x=old.start, xend=old.end, y=fam, yend=fam, color=len)) + scale_colour_gradient(high = "red", low = "orange") +
#        geom_segment(size=2) + geom_point() + ggtitle(sprintf("Discordant Sites\nInterval > %s nt",MIN.OVLP.LEN)) + xlm)

# display density of CNVs
z2<-
  ldply(c(TRUE,FALSE), function(conc) {
    ldply(c(1500,2000,5000), function(nt.len.min){
      cat(conc,nt.len.min,"\n")
      data.frame(conc=conc, min.len=nt.len.min, location=unlist(dlply(filter(sa1, len>nt.len.min & concordant==conc), .(family, chr, start), function(r) seq(r$old.start,r$old.end)), use.names = FALSE))
    })
  })
#print(plot(density(z), main=sprintf("Density of CNVs >%snt chr 20",MIN.OVLP.LEN)))
#z2$min.length <- as.factor(z2$min.len)
#ggplot(z2, aes(x=location, color=min.length)) + geom_density()
z3 <- tally(group_by(z2, min.len, conc, location))

# Expensive plot:
#ggplot(filter(z3), aes(x=location, y=n, color=conc)) + geom_step() + facet_grid(min.len+conc~., scales = "free_y")

# display location of CNVs. Only the largest CNVs are visible
#print(ggplot(filter(segs, !all.eq), aes(x=start, xend=end, y=fam, yend=fam, color=all.eq)) + geom_segment(size=10) + theme_bw() + geom_blank() + theme(panel.grid.major.y=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor.y = element_blank(), legend.position="none"))

# manually select the top 8 regions
z4 <- filter(z3, min.len==2000)
z4p <- function(z, a,b, i) {
  print(ggplot(z, aes(x=location, y=n, color=conc)) +geom_point() + xlim(a,b) + ylab('n') + xlab('pos') + ggtitle(sprintf("#%s Chr 20:%s-%s (%snt)",i, a,b,b-a)))
}
bad.loci <- list(
  c(6.05195e7,6.05235e7,1),
  c(4.1243e7,4.125e7,2),
  c(2.606e7,2.6074e7,3),
  c(2.9493e7,2.951e7,4),
  c(3.709e7,3.7095e7,5),
  c(1.2341e7, 1.23445e7,6),
  c(5.4325e6,5.4355e6,7),
  c(1.585e6,1.592e6,8))

bad.loci <- list(
  c(2.949e7, 2.951e7,1),  # a high frequency all concordant locus to the left, this is a 5000nt tail region with a few more discordant than concordant
  c(0.1584e7, 0.1591e7,2),# 4000nt part of almost 100% frequency CNV where about 1/3 are discordant
  c(4.1243e7, 4.125e7,3), # more than 40 concordant and less than 5 discordant
  c(3.3242e7, 3.3245e7,4) # these are < 1400nt
)

# pdf("discordant_loci.pdf")
# lapply(bad.loci, function(ab) {
#   z4p(z4, ab[1],ab[2], ab[3])
# })
# dev.off()

# count concordance by base and cnv
len.max <- 15000
len.min.disp <- 1200
#steps <- c(seq(0,5000,100), seq(6000,len.max,1000))
steps <- c(seq(0,50,1), seq(60,len.max/100,10))

cts <- function(segs.subset) {
  ldply(steps, function(len.min) {
    z <- ddply(filter(segs.subset, !is.na(ibd.state) & end.bin-start.bin > len.min), .(cnv), function(df) {
      mutate(data.frame(len.min,tot.bases=sum(as.numeric(df$len)), 
                        conc.bases=sum(as.numeric(df$len[df$concordant])),
                        ncnv=nrow(df), 
                        conc.cnv=sum(df$concordant)), 
             conc.base.pct=conc.bases/tot.bases,
             conc.cnv.pct=conc.cnv/ncnv)
    })
    if (nrow(z)==0) { return(data.frame())} else { return(z)}
  })
}

diff.cts <- function(x) {
  idx <- 1:(length(x)-1)
  return(c(x[idx]-x[idx+1],x[length(x)]))  
}

bin.cts.df <- function(df) {
  ddply(df, .(grp, family, cnv), 
        mutate, tot.bases.bin=diff.cts(tot.bases), conc.bases.bin=diff.cts(conc.bases),
        ncnv.bin=diff.cts(ncnv), conc.cnv.bin=diff.cts(conc.cnv), disc.cnv.bin=ncnv.bin-conc.cnv.bin)
        
}

# compute a normalized area under the curve given the points (x,y)
# Assume that y is [0,1], x is arbitrary.
auc <- function(x,y) {
  return(sum(diff(x)*rollmean(y,2))/(x[length(x)]-x[1]))
}


all.cts <- mutate(cts(segs), family='ALL', grp='All')
ggplot(filter(all.cts, cnv == 'DEL'), aes(x=len.min, y=conc.cnv.pct)) + geom_point() + geom_line() + geom_vline(xintercept = 12, linetype=2) + xlim(0,50)+ylim(0.8,1) + ggtitle("Posterior, No Prior N=1")
ggplot(filter(all.cts), aes(x=len.min, y=conc.cnv.pct, color=cnv)) + geom_point() + geom_line() + geom_vline(xintercept = 12, linetype=2) + xlim(0,50) + ggtitle("Concordancy By CNV Type")

all.cts <- mutate(all.cts, is.cnv = cnv %in% c('BI','DEL','MULTI'), is.wt=='WT') # is.wt is NA for cnv==NA, but is.cnv is not

cts.wt <- ddply(all.cts, .(is.wt, len.min), summarize, tot.bases=sum(tot.bases), conc.bases=sum(conc.bases), 
                ncnv=sum(ncnv), conc.cnv=sum(conc.cnv), cnv.base.pct=conc.bases/tot.bases, conc.cnv.pct=conc.cnv/ncnv)
ggplot(cts.wt, aes(x=len.min, y=conc.cnv.pct, color=is.wt)) + geom_point() + geom_line() + geom_vline(xintercept = 12, linetype=2) + xlim(0,50) + ggtitle("Concordancy: WT or CNV")

all.cts2 <- ddply(filter(all.cts, !is.na(cnv)), .(len.min), summarize, tot.bases=sum(tot.bases), conc.bases=sum(conc.bases), 
                  ncnv=sum(ncnv), conc.cnv=sum(conc.cnv), cnv.base.pct=conc.bases/tot.bases, conc.cnv.pct=conc.cnv/ncnv)
ggplot(all.cts2, aes(x=len.min, y=conc.cnv.pct)) + geom_point() + geom_line() + geom_vline(xintercept = 12, linetype=2) + xlim(0,50) + ggtitle("Total Base Concordancy")


# remove segments from df that overlap range ab
excl.loci <- function(df,ab) {
  filter(df, !(start < ab[2] & end > ab[1]))
}

# generate cts for each bad loci independently that can be plotted or used for AUC
rocs <- llply(bad.loci, function(ab) {
  mutate(cts(excl.loci(segs,ab)),grp=sprintf("All but #%s (20:%s-%s)",ab[3],ab[1],ab[2]))
})
rocs[['All']] <- mutate(cts(segs),grp="All")
rocs.hold1 <- do.call(rbind, rocs)
rocs.hold1$all <- rocs.hold1$grp=='All'
ggplot(filter(rocs.hold1, cnv %in% c('DEL')), 
             aes(x=len.min, y=conc.cnv.pct, color=grp))+ geom_line() + facet_grid(cnv~.) + xlab("Min Segment Size") + ylab("Concordance") + ggtitle("Quartet Concordance By CNV")  + geom_vline(xintercept = 12, linetype=2) + xlim(0,30)+ylim(0.8,1) + guides(size="none")

# generate AUC for each of the "hold one out" ROCs
auc.calc <- function(df) {
  ddply(df, .(cnv,grp), summarize, auc=auc(len.min, conc.cnv.pct))
}
auc.res <- ldply(rocs, auc.calc)

# Ranking the results by their independent impact on AUC for DELs
# > auc.res %>% filter(cnv=='DEL') %>% arrange(auc)
# .id cnv                       grp       auc
# 1     DEL   All but 1585000-1592000 0.9646529
# 2 All DEL                       All 0.9726434
# 3     DEL All but 12341000-12344500 0.9728125
# 4     DEL   All but 5432500-5435500 0.9729876
# 5     DEL All but 41243000-41250000 0.9733802
# 6     DEL All but 37090000-37095000 0.9743372
# 7     DEL All but 60519500-60523500 0.9744237
# 8     DEL All but 29493000-29510000 0.9749306
# 9     DEL All but 26060000-26074000 0.9758405


excl.order <- rev(c(6,7,2,5,1,4,3))

# put #5 first
excl.order <- rev(c(6,7,2,1,4,3,5))

# > auc.res %>% filter(cnv=='DEL') %>% arrange(auc)
# .id cnv                               grp       auc
# 1     DEL   All but #2 (20:1584000-1591000) 0.9861576
# 2 All DEL                               All 0.9865080
# 3     DEL All but #4 (20:33242000-33245000) 0.9865646
# 4     DEL All but #3 (20:41243000-41250000) 0.9866056
# 5     DEL All but #1 (20:29490000-29510000) 0.9868966
excl.order <- c(1,4,3) # 2 not included because worsens ROC

segs.ex <- list(segs)
lapply(excl.order, function(i) {
  segs.ex[[length(segs.ex)+1]] <<- excl.loci(segs.ex[[length(segs.ex)]], bad.loci[[i]])  
})
excl.label <- c('All', sprintf('All but %s', excl.order[1]))
lapply(excl.order[2:length(excl.order)], function(i) {
  excl.label <<- c(excl.label, sprintf("%s&%s", excl.label[length(excl.label)],i))
})

rocs2 <- ldply(1:length(excl.label), function(i) {
  mutate(cts(segs.ex[[i]]),grp=excl.label[i])
})

p5 <- ggplot(filter(rocs2, cnv %in% c('DEL')), 
             aes(x=len.min, y=conc.cnv.pct, color=grp))+ geom_line() + facet_grid(cnv~.) + xlab("Min Segment Size") + ylab("Concordance") + ggtitle("Quartet Concordance By CNV")  + geom_vline(xintercept = 12, linetype=2)
print(p5)
print(p5 + xlim(c(5,30)) + coord_cartesian(ylim=c(0.8,1)))  # + ggtitle("Rank #5 first"))

auc.calc(rocs2) %>% filter(cnv=='DEL') %>% arrange(auc)

p6 <- ggplot(filter(rocs3, cnv %in% c('DEL') & grp=='All'), 
       aes(x=len.min, y=conc.cnv.pct, color=method))+ geom_line() + facet_grid(cnv~.) + xlab("Min Segment Size") + ylab("Concordance") + ggtitle("CNV Concordance By Method")  + geom_vline(xintercept = 12, linetype=2)
print(p6 + xlim(c(5,30)) + coord_cartesian(ylim=c(0.8,1)))

segs$len.bin <- segs$end.bin-segs$start.bin
write.table(select(segs, family, sib1, sib2, mother, father, chr, sib1.cn, sib2.cn, par1.cn, par2.cn, ibd.state, concordant, cnv, start.bin, end.bin, len.bin), file="segs.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#######################
write.table(select(segs.map.10, family, sib1, sib2, mother, father, chr, sib1.cn, sib2.cn, par1.cn, par2.cn, ibd.state, concordant, cnv, start.bin, end.bin, len.bin), file="segs-map.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# rocs.mle <- rocs2; segs.mle <- segs; segs.mle$len.bin <- segs.mle$end.bin - segs.mle$start.bin
# save(rocs.mle, segs.mle, file="/cygwin64/tmp/mle_rocs2.Rdata")
write.table(select(segs.mle.10, family, sib1, sib2, mother, father, chr, sib1.cn, sib2.cn, par1.cn, par2.cn, ibd.state, concordant, cnv, start.bin, end.bin, len.bin), file="segs-mle.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
 
#################
# distribution of non-2 lengths and gaps
cnv.lens <- ddply(segs, .(family,all.wt), summarize, len=sum(end-start))
ggplot(na.omit(cnv.lens), aes(x=len/4, fill=all.wt)) + geom_histogram() + ggtitle("Amount of WT vs non-WT") + xlab('Length of Segments')
#################

# Retrieve the discordant CNVs in SRPB1
srpb1 <- filter(as.tbl(segs), end > 1552936 & start < 1594100 & end.bin-start.bin > 12)
srpb1.disc.fam <- unique(filter(srpb1, !concordant)$family)

label = c('sib1','sib2','mother','father')
label.cn = c('sib1.cn','sib2.cn','par1.cn','par2.cn')
srpb1_m <-
  ldply(1:4, function(idx) {
    mutate(select_(srpb1, 'family', sample=label[idx], 'chr','start','end', cn=label.cn[idx], 'ibd.state', 'concordant'), who=label[idx])
  })

# load the gold standard genotypes
srpb1.regions <- c('SD1','LEFT','CENTRAL','RIGHT')
srpb1.regions.start <- list(SD1=1552936, LEFT=1556802, CENTRAL=1561002, RIGHT=1585902)
srpb1.regions.end <- list(SD1=1556761, LEFT=1561000, CENTRAL=1585900, RIGHT=1594100)
srpb1_geno <-
  ldply(srpb1.regions, function(n) {
    mutate(read.table(file=sprintf("%s/%s.cn.dat", basedir, n), header = TRUE, stringsAsFactors = FALSE), 
           region=n, start=as.integer(srpb1.regions.start[[n]]), end=as.integer(srpb1.regions.end[[n]]))
  })

# Add the genotype and region label to each srpb1.m line
srpb1 <- sqldf("SELECT sm.*, sg.start as rstart, sg.end as rend, sg.CN as gold_cn, sg.cnq, sg.region FROM srpb1_m sm, srpb1_geno sg 
          WHERE sg.SAMPLE=sm.sample AND sm.start<sg.end AND sm.end>sg.start", drv='SQLite')
srpb1 <- mutate(srpb1, a=ifelse(start>rstart, start, rstart), b=ifelse(end<rend, end, rend))
arrange(filter(srpb1, !concordant & cn!=gold_cn), family, a, sample) %>% select(family, sample, ibd.state, who, region, a, b,cn,gold_cn) %>% mutate(len=b-a) %>% filter(len>1000)

# Then, for each segment, find all 4, e.g.
arrange(filter(srpb1, family==11470), family, a, sample) %>% select(family, sample, ibd.state, who, region, a, b,cn,gold_cn) %>% mutate(len=b-a) %>% filter(len>1000)

##################

all.cts.m2 <- mutate(cts(filter(segs, !(family %in% c(11711, 11393)))), family='ALL', grp='All but 2')
all.cts2 <- mutate(cts(filter(segs, !excl)), family='ALL', grp='No Big 3')
all.cts2a <- mutate(cts(filter(segs, !excl.alpha)), family='ALL', grp='No Alpha')

all.cts3 <- rbind(all.cts, all.cts.m2, all.cts2, all.cts2a)

ddply(all.cts3, .(grp,cnv), function(df) auc(df$len.min, df$conc.cnv.pct))



p1c <- ggplot(filter(all.cts3, cnv %in% c('DEL')), aes(x=len.min, y=conc.cnv.pct, color=grp))+ geom_line() + facet_grid(cnv~.) + xlab("Min Segment Size") + ylab("Concordance") + ggtitle("Quartet Concordance By CNV")  + geom_vline(xintercept = 12, linetype=2)
print(p1c)
# print(p1c + xlim(500,5000) + coord_cartesian(ylim=c(0.7,1)))
print(p1c + xlim(5,50) + coord_cartesian(ylim=c(0.7,1)) + xlab("Min Segment Size (Bins)"))

p1b <- ggplot(filter(all.cts, cnv %in% c('DEL')), aes(x=len.min, y=conc.base.pct)) + geom_point() + geom_line() + facet_grid(cnv~.) + xlab("Min Segment Size (Bins)") + ylab("Concordance") + ggtitle("Quartet Concordance By Base")
print(p1b)

all.cts3.bin <- bin.cts.df(all.cts3)

p1f <- ggplot(filter(all.cts3.bin, cnv %in% c('DEL')), aes(x=len.min, y=disc.cnv.bin, color=grp))+ geom_line() + facet_grid(cnv~.) + xlab("Segment Size (Bin)") + ylab("# Discordant CNVs") + ggtitle("CNV Discordance By Length")  + geom_vline(xintercept = 12, linetype=2)
print(p1f + xlim(5,60) + coord_cartesian(ylim=c(0,50)))

quartet.cts <-
  ldply(quartets$family, function(f) {
    mutate(cts(filter(segs, family==f)), family=f, grp='All')
  })
q.cts2a <- ldply(quartets$family, function(f) {
  mutate(cts(filter(segs, family==f & !excl.alpha)), family=f, grp='No Alpha')
})
cts.m <- rbind(quartet.cts, all.cts)


fams <- unique(segs$family)
pdf("quartets.pdf")
sapply(seq(1,length(fams),6), function(f) {
  p1g <- ggplot(filter(quartet.cts, family %in% fams[f:(f+5)] & cnv=='DEL'), aes(x=len.min, y=conc.cnv.pct,group=family)) + geom_point() +geom_line()+ facet_grid(family+cnv~.) + xlab("Min Segment Size") + ylab("Concordance") + ggtitle("Quartet Concordance By CNV")
  print(p1g)
})
dev.off()


p2 <- ggplot(filter(all.cts, cnv=='DEL' & len.min > len.min.disp), aes(x=len.min, y=ncnv)) + geom_point() + facet_grid(cnv~.) + xlab("Min Segment Size") + ylab("Number of CNVs") + ggtitle("Quartet CNV Count")


# as smaller segments are removed, would the flanking segments be merged because they're the same CN?
# for each family, iteratively remove segments below a minimum length, and then check the flanking genotypes for remaining
segs2 <- ddply(segs, .(family), function(df) {
  # df <- filter(segs, family==11006)
  cat(df$family[1],"\n")
  ldply(seq(len.min.disp,len.max,100), function(len.min) {
    # len.min <- 500
    df2 <- filter(df, len > len.min)
    idx <- 2:nrow(df2)
    df2$sib1.adj.eq <- c(df2$sib1.cn[idx] == df2$sib1.cn[idx+1] ,NA)  
    df2$sib2.adj.eq <- c(df2$sib2.cn[idx] == df2$sib2.cn[idx+1] ,NA)  
    df2$par1.adj.eq <- c(df2$par1.cn[idx] == df2$par1.cn[idx+1] ,NA)  
    df2$par2.adj.eq <- c(df2$par2.cn[idx] == df2$par2.cn[idx+1] ,NA)  
    df2$all.adj.eq <- df2$sib1.adj.eq & df2$sib2.adj.eq & df2$par1.adj.eq & df2$par2.adj.eq
    df2$len.min <- len.min
    df2
  })
})



# define "CNV region" as the chromosome segment that's flanked by wild type genotypes, i.e. CN=2 for all samples.
# collect stats on these regions as minimum length increases. How many CNVs and CNV regions?
# Ideally, as the minimum length increases, the ratio of CNVs to CNV regions approaches one.
cnv.regs <- function(s,grp) {
  mutate(ddply(s, .(family), function(df) {
    # df <- filter(segs, family==11835)
    ldply(steps, function(len.min) {
      df2 <- filter(df, len > len.min)
      df2$all.wt.adj <- c(FALSE, df2$all.wt[2:nrow(df2)] & df2$all.wt[1:(nrow(df2)-1)])
      
      # the number of regions is the number of WT CNVs minus one (due to boundary)
      cnv.region.ct <- sum(df2$all.wt) - sum(df2$all.wt.adj) - 1
      
      # the number of CNVs is all rows that are not WTs
      cnv.ct <- sum(!df2$all.wt)
      
      # concordant CNVs
      cnv.conc.ct <- sum(df2$concordant[!df2$all.wt], na.rm=TRUE)
      
      data.frame(len.min=len.min, cnv.region.ct=cnv.region.ct, cnv.ct=cnv.ct, cnv.conc.ct=cnv.conc.ct)
    })
  }), fam=as.factor(family), grp=grp)
}

cnv.regs1 <- cnv.regs(segs,'ALL but 2')
cnv.regs2 <- cnv.regs(filter(segs, !excl.alpha),'No Alpha')
cnv.regs3 <- cnv.regs(filter(segs, !excl), 'No Big 3')

cnv.regs123 <- rbind(mutate(cnv.regs1, grp='ALL'), filter(rbind(cnv.regs1, cnv.regs2, cnv.regs3), !(family %in% c(11711, 11393))))

cnv.regs.all <- ddply(cnv.regs123, .(len.min,grp), summarize, cnv.region.ct=sum(cnv.region.ct), cnv.ct=sum(cnv.ct), cnv.conc.ct=sum(cnv.conc.ct))
cnv.regs.all$gt.1200 <- cnv.regs.all$len.min >= 1200
cnv.regs.all$conc <- with(cnv.regs.all, cnv.conc.ct/cnv.ct)
cnv.regs.all$conc.ranges <- cut(cnv.regs.all$conc, c(0,0.8,0.85,0.9,0.95,1))

p4 <- ggplot(cnv.regs.all, aes(x=len.min, y=cnv.region.ct/cnv.ct,color=grp)) + geom_line() + ggtitle("Regions / CNVs")
print(p4)
print(p4 + xlim(500,5000) + coord_cartesian(ylim=c(0.7,1)))

#print(p4 + geom_point(aes(size=conc.ranges)) + xlim(500,5000) + coord_cartesian(ylim=c(0.7,1)))
ggplot(cnv.regs.all, aes(y=cnv.region.ct, x=cnv.ct, group=grp, color=gt.1200)) + geom_point() + geom_line() + ggtitle("Regions vs CNVs")
ggplot(cnv.regs.all, aes(y=cnv.region.ct, x=cnv.ct, color=grp)) + geom_point() + geom_line() + ggtitle("Regions vs CNVs")

ggplot(cnv.regs1, aes(x=len.min, y=cnv.region.ct/cnv.ct, color=fam)) + geom_point() + geom_line() + guides(color=FALSE)  + facet_wrap(~fam)


segs2.ct <- mutate(ddply(filter(segs2, !is.na(all.adj.eq)), .(len.min), summarize, tot=length(all.adj.eq), all.adj.eq=sum(all.adj.eq), frac.eq=all.adj.eq/tot), cnv='DEL')
#ggplot(segs2.ct, aes(x=len.min, y=all.adj.eq)) + geom_point()
p3 <- ggplot(segs2.ct, aes(x=len.min, y=frac.eq)) + geom_point() + xlab("Min Segment Size") + ylab("Fraction of Segments") + ggtitle("Segments with Identical Flanking CN") + facet_grid(cnv~.)

grid.newpage()
XL <- xlim(len.min.disp, 20000)
grid.draw(rbind(ggplotGrob(p1+XL), ggplotGrob(p2+XL), ggplotGrob(p3+XL), size = "last"))


ggplot(filter(segs, len>2000 & cnv=='DEL' & !is.na(concordant)), aes(x=start, y=fam, xend=end, yend=fam)) + geom_point() + geom_segment() + geom_point(aes(x=end)) + facet_grid(concordant+cnv~.) + ggtitle("segments > 2000 nt")

# As segments increase, concordance tends to rise
# Looking across the genome, at a high level it's clear that there are a few sites (three major ones) where there
# is high frequency and much disagreement.

# Here's a zoom into the first such site.
ggplot(filter(segs, !all.wt & start < 2.65e7 & start > 2.6e7), aes(x=start, y=fam, xend=end, yend=fam, color=concordant)) + geom_point() + geom_segment() + geom_point(aes(x=end)) + facet_grid(concordant~.)

# I checked these regions and only one of the three is an alpha satellite.

# identify sites:
# All segments are collapsed and each group is called a site
# Each site has a length, membership & frequency, and consistency
all.sibs <- cbind(rbind(cnvs.sibs[[1]], cnvs.sibs[[2]]), sib=TRUE)
all.pars <- cbind(rbind(cnvs.parents[[1]], cnvs.parents[[2]]), sib=FALSE)
cnvs.all <- rbind(all.sibs, all.pars)

# overlap(df) - returns a flattened data.frame across all samples
overlap <- function(df, start.name='start.map', end.name='end.map', seg='seg') {
  if (nrow(df)==1) {
    return(data.frame(i=1,j=1,start.map=df$start.map,end.map=df$end.map,
                      x.min=df$start.map, x.max=df$start.map,
                      y.min=df$end.map, y.max=df$end.map,
                      x.diff=0, y.diff=0, n=1))
  }
  df2 <- arrange(df, start.map, end.map)
  i2n <- 2:nrow(df2)
  
  # perform closure. repeat until no change.
  df2$has.ovlp2 <- rep(FALSE, nrow(df2))
  df2$end.max2 <- rep(0, nrow(df2))
  
  df2$end.max <- df2$end.map
  again <- TRUE
  while (again) {
    df2$has.ovlp <- c(FALSE, df2$start.map[i2n] < df2$end.max[i2n-1])
    df2$end.max <- c(df2$end.map[1], ifelse(df2$has.ovlp[i2n], pmax(df2$end.max[i2n],df2$end.max[i2n-1]), df2$end.max[i2n]))
    again <- !all(df2$has.ovlp==df2$has.ovlp2) || !all(df2$end.max==df2$end.max2)
    df2$has.ovlp2 <- df2$has.ovlp
    df2$end.max2 <- df2$end.max
  }
  
  # starts are where there is no overlap
  # ends are the end.max in the row before the starts. skip row 1. Add end.max of row N
  i <- which(!df2$has.ovlp)
  j <- c(which(!df2$has.ovlp[i2n]), nrow(df2))
  starts <- df2$start.map[i]
  ends <- df2$end.max[j]
  
  df.merged <- data.frame(i=i, j=j, start.map=starts, end.map=ends)

  df.merge <- function(x, i, j) {
    mutate(x, i, j, x.min=min(df2$start.map[i:j]), 
           x.max=max(df2$start.map[i:j]),
           y.min=min(df2$end.map[i:j]), y.max=max(df2$end.map[i:j]),
           x.diff=x.max-x.min, y.diff=y.max-y.min,
           n=j-i+1,
           xs=paste0(df2$start.map[i:j], collapse=';'),
           ys=paste0(df2$end.map[i:j], collapse=';'),
           seg=paste0(df2[[seg]][i:j], collapse = ';'),
           cns=paste0(df2$cn[i:j], collapse=';'),
           samples=paste0(df2[['.id']][i:j], collapse=';'))
  }
  df.merged <- ddply(df.merged, .(i,j), df.merge)
  
  return(df.merged)
  
}

# collapse CNVs
cnvs.c <- as.tbl(mutate(ddply(filter(cnvs.all, cn!=2), .(chr), overlap, seg='label')))
cnvs.c$len <- cnvs.c$end.map - cnvs.c$start.map
cnvs.c$site.id <- 1:nrow(cnvs.c)
cnvs.c$seg <- NULL
cnvs.c <- ddply(cnvs.c, .(site.id), function(df) {
  xs <- as.numeric(unlist(strsplit(df$xs, ';')))
  ys <- as.numeric(unlist(strsplit(df$ys, ';')))

  xs.med <- median(xs)
  ys.med <- median(ys)

  offsets <- c(xs - df$start.map, df$end.map - ys)
  stopifnot(all(offsets>=0))
  
  sd1 <- sd(offsets)
  outliers <- abs(xs - xs.med) > 2*sd1 | abs(ys-ys.med) > 2*sd1
  xs2 <- xs[!outliers]
  ys2 <- ys[!outliers]
  offsets2 <- c(xs2 - df$start.map, df$end.map - ys2)
  if (length(offsets2) < 2) { sd2 <- 0 } else { sd2 <- sd(offsets2) }
  outlier.n <- length(xs)-length(xs2)
  return(mutate(df, sd1=sd1, sd2=sd2, outlier.n=outlier.n))
})

# create a map of seg (i.e. site ID to sample name of members)
site.map <- ddply(cnvs.c, .(site.id), function(df) {
  xs <- as.numeric(unlist(strsplit(df$xs, ';')))
  ys <- as.numeric(unlist(strsplit(df$ys, ';')))
  samples <- unlist(strsplit(df$samples, ';'))
  cns <- as.numeric(unlist(strsplit(df$cns, ';')))

  return(data.frame(chr=df$chr, start.map=df$start.map, end.map=df$end.map, sample=samples, cn=cns, x=xs, y=ys, stringsAsFactors = FALSE))
})

site.map <- as.tbl(merge(site.map, family.sample))

# perform concordance. for each site, for each quartet, check concordance. If a sample is missing, then cn=2
site.concordance <-
  ddply(site.map, .(site.id, family), function(site.fam) {
    # in general there are 0, 1 or more rows per sample. If there is 0 then that sample is CN=2.
    # If 1, then a single CN!=2 for that site, e.g. CN=1.
    # If >1, then there might be multiple CNs for that sample within the site (e.g. CN=3 and CN=4)
    # or due to an outstanding bug where there are adjacent segments with the same CN.
    
    # site.fam <- filter(site.map, site.id==30 & family==11348)
    # site.fam <- filter(site.map, site.id==30 & family==11420)
    # site.fam <- filter(site.map, site.id==93 & family==11218)
    site.start <- site.fam$start.map[1]
    site.end <- site.fam$end.map[1]
    
    # get the sib1 and sib2 sample names from family.sample
    samples <- filter(family.sample, family==site.fam$family[1])
    sib1.sample <- samples$sample[samples$who=='sib1']
    sib2.sample <- samples$sample[samples$who=='sib2']
    
    site.fam <- merge(site.fam, samples, all.y=TRUE)
    missing <- is.na(site.fam$site.id)
    site.fam$cn[missing] <- 2
    site.fam$x[missing] <- site.start
    site.fam$y[missing] <- site.end
    
    # HACK! Punt on the multiple CN!=2 segments per sample as this point. FIXME.
    site.fam <- ddply(site.fam, .(who), function(x) x[1,])
    
    # filter ibd on sib1, sib2, chr and position to derive the IBD state
    this.ibd <- filter(ibd, ((SIB1==sib1.sample & SIB2==sib2.sample) | (SIB1==sib2.sample & SIB2==sib1.sample)) & START < site.end & END > site.start)
    ibd.state <- ibd.code(this.ibd)
    
    parm.cn <- site.fam$cn[site.fam$who == 'mother']
    parp.cn <- site.fam$cn[site.fam$who == 'father']
    sib.sum <- sum(site.fam$cn[site.fam$who %in% c('sib1','sib2')])
    par.sum <- sum(site.fam$cn[site.fam$who %in% c('mother','father')])
    sums.eq <- sib.sum == par.sum
    sib.eq <- site.fam$cn[site.fam$who == 'sib1'] == site.fam$cn[site.fam$who == 'sib2']
    
    # IBD1P
    sib.sum.less.m <- sib.sum - parm.cn
    sib.sum.less.m.even <- sib.sum.less.m %% 2 == 0
    sib.sum.less.m.div2.ge0 <- sib.sum.less.m / 2 >= 0
    sib.sum.less.m.div2.lef <- sib.sum.less.m / 2 <= parp.cn
    ibd1p.conc <- sib.sum.less.m.even & sib.sum.less.m.div2.ge0 & sib.sum.less.m.div2.lef
    
    # IBD1M
    sib.sum.less.p <- sib.sum - parp.cn
    sib.sum.less.p.even <- sib.sum.less.p %% 2 == 0
    sib.sum.less.p.div2.ge0 <- sib.sum.less.p / 2 >= 0
    sib.sum.less.p.div2.lem <- sib.sum.less.p / 2 <= parm.cn
    ibd1m.conc <- sib.sum.less.p.even & sib.sum.less.p.div2.ge0 & sib.sum.less.p.div2.lem
    
    site.concordance <- ifelse(ibd.state=='IBD0', sums.eq, ifelse(ibd.state=='IBD2', sib.eq, ifelse(ibd.state=='IBD1P', ibd1p.conc, ifelse(ibd.state=='IBD1M', ibd1m.conc, NA))))
    if (length(site.concordance)==0 || is.na(site.concordance)) {
      return(data.frame())
    } else {
      return(data.frame(site.concordance=site.concordance, site.len=site.end-site.start))
    }
  })

site.concordance <- arrange(site.concordance, -site.len)
site.concordance$trues <- cumsum(site.concordance$site.concordance)
site.concordance$site.fams <- 1:nrow(site.concordance)
site.concordance$frac.conc <- site.concordance$trues / site.concordance$site.fams


# treat each breakpoint as an entity instead of pairs of breakpoints representing sites.
# each breakpoint has a position, a left and right copy number, a confidence interval and ???
bp.all <- ddply(cnvs.all, .(.id, chr), function(df) {
  bp <- select(df[2:nrow(df),], sample=.id, chr, CI.L=start.CI.L, bp.map=start.map, CI.R=start.CI.R, cn.R=cn, sib)
  bp$cn.L <- bp$cn.L <- c(2, bp$cn.R[1:(nrow(bp)-1)])
  return(bp)
}) %>% select(-.id)

# > head(bp.all)
# .id chr    CI.L  bp.map    CI.R cn.R  sib cn.L
# 2 SS0013018  20  700036  700236  700336    1 TRUE    2
# 3 SS0013018  20  701480  701480  701480    2 TRUE    1
# 4 SS0013018  20 1388308 1389108 1389444    1 TRUE    2
# 5 SS0013018  20 1390844 1390844 1390844    2 TRUE    1
# 6 SS0013018  20 1560729 1561029 1561247    1 TRUE    2
# 7 SS0013018  20 1594059 1594059 1594059    2 TRUE    1

cnvsAll <- rename(cnvs.all, sample=.id) # no dots for sqldf!

# So what about the samples with no breakpoint? In order to measure bp concordance We still need to know what the CN is.
# For each quartet, cluster by overlapping CIs. For each set of BPs, fill in missing by querying for the overlapping segment.
bp <-
ddply(family.sample, .(family), function(fam) {
  # fam <- filter(family.sample, family==11006)
  q.bp <- arrange(merge(bp.all, fam), CI.L, CI.R)
  if (nrow(q.bp) > 0) {
    # here we go again. Another collapse. This time by CI range
    # this time, scan through, setting a cluster ID on each row, incrementing the ID when it does not overlap the previous one
    idx <- 2:nrow(q.bp)
    max.R <- q.bp$CI.R
    ovlp2 <- rep(FALSE, nrow(q.bp)-1)
    ovlp <- max.R[idx-1] >= q.bp$CI.L[idx]
    while (!all(ovlp==ovlp2)) {
      max.R <- c(max.R[1], ifelse(ovlp, pmax(max.R[idx], max.R[idx-1]), max.R[idx]))
      ovlp2 <- ovlp
      ovlp <- max.R[idx-1] >= q.bp$CI.L[idx]
    }
    cluster.starts <- c(1,which(!ovlp)+1)
    q.bp$cluster <- 0
    q.bp$cluster[cluster.starts] <- 1
    q.bp$cluster <- cumsum(q.bp$cluster)
    
    
    # for each cluster, complete the family by adding CN.L=2, CN.R=2 for missing samples
    q.bp <- ddply(q.bp, .(cluster), function(q.bp.c) {
      # q.bp.c <- filter(q.bp, cluster==2)
      if (nrow(q.bp.c) < 4) {
        CI.L <- min(q.bp.c$CI.L)
        CI.R <- max(q.bp.c$CI.R)
        missing <- !(fam$sample %in% q.bp.c$sample)
        missing.samples <- fam$sample[missing]
        missing.sib <- fam$who[missing] %in% c('sib1','sib2')
        missingTbl <- data.frame(sample=missing.samples)
        
        # retrieve the CN for the missing samples in the current range
        missing.cn <- sqldf(sprintf("SELECT c.sample, c.cn FROM cnvsAll c, missingTbl m WHERE c.sample=m.sample AND \"start.map\" < %s AND \"end.map\" > %s", CI.R, CI.L), drv='SQLite')
        
        if (nrow(missing.cn)+nrow(q.bp.c) == 4) {
          
          # combine bp and missing bp together
          rbind(q.bp.c,
                data.frame(sample=missing.cn$sample, chr=q.bp.c$chr[1], CI.L=CI.L, bp.map=NA, CI.R=CI.R, cn.R=missing.cn$cn, 
                           sib=missing.sib, cn.L=missing.cn$cn, family=q.bp.c$family[1], who=fam$who[missing], cluster=q.bp.c$cluster[1]))
          
        } else {
          # there could be a number of different things going on here.
          # Most commonly the location is passed the "end" of the annotated chromosome for a sample.
          # Sometimes there are multiple regions for the same sample, which shouldn't happen!
          if (any(q.bp.c$bp.map < 62900000)) {
            cat("Problem entry:\n")
            print(q.bp.c)
            print(missing.cn)
          }
          if (!(all(table(missing.cn$sample) == 1))) {
            cat("Multiple entries for some samples!\n")
          }
          return(data.frame())
        }
        
      } else {
        return(q.bp.c)  
      }
    })
    
  } else { return(data.frame()) }
})

# Add family ID to ibd table
sib.pairs <- rbind(data.frame(family=quartets$family, PAIR_NAME=paste0(quartets$sib1,'-',quartets$sib2)),
                   data.frame(family=quartets$family, PAIR_NAME=paste0(quartets$sib2,'-',quartets$sib1)))
ibd <- merge(sib.pairs, ibd)

# Extract min/max range for each bp
bpbounds <- ddply(bp, .(family,cluster), summarize, start=min(CI.L), end=max(CI.R))

# join with ibd to label each bp with an IBD code
z <- sqldf("select b.*, i.STATE, i.PATERNAL_IND, i.MATERNAL_IND from bpbounds b, ibd i where b.family=i.family and b.start <= i.END and b.end >= i.START", drv='SQLite')
z$ibd.code <- ibd.code(z)

bp <- merge(bp, subset(z, select=c(family, cluster, ibd.code)))
# bp:
#    family cluster   sample chr    CI.L  bp.map    CI.R cn.R   sib cn.L    who ibd.code
# 1   11006       1 SSC00006  20 1388908 1389108 1389208    1  TRUE    2   sib2     IBD0
# 2   11006       1 SSC00005  20 1388908 1389108 1389444    0 FALSE    2 mother     IBD0
# 3   11006       1 SSC00003  20 1388308 1389108 1390144    1  TRUE    2   sib1     IBD0
# 4   11006       1 SSC00004  20 1388308      NA 1390144    2 FALSE    2 father     IBD0
# 5   11006      10 SSC00003  20 1995267      NA 1995567    2  TRUE    2   sib1     IBD0
# 6   11006      10 SSC00005  20 1995267      NA 1995567    2 FALSE    2 mother     IBD0
# 7   11006      10 SSC00004  20 1995267 1995467 1995567    1 FALSE    2 father     IBD0
# 8   11006      10 SSC00006  20 1995267 1995367 1995467    1  TRUE    2   sib2     IBD0
# 9   11006      11 SSC00005  20 1997485      NA 1997485    2 FALSE    2 mother     IBD0
# 10  11006      11 SSC00004  20 1997485      NA 1997485    2 FALSE    2 father     IBD0

concordance.stats <- function(df, cn) {
  ibd.state <- df$ibd.code[1]
  
  parm.cn <- df[[cn]][df$who == 'mother']
  parp.cn <- df[[cn]][df$who == 'father']
  sib.sum <- sum(df[[cn]][df$who %in% c('sib1','sib2')])
  par.sum <- sum(df[[cn]][df$who %in% c('mother','father')])
  sums.eq <- sib.sum == par.sum
  sib.eq <- df[[cn]][df$who == 'sib1'] == df[[cn]][df$who == 'sib2']
  
  # IBD1P
  sib.sum.less.m <- sib.sum - parm.cn
  sib.sum.less.m.even <- sib.sum.less.m %% 2 == 0
  sib.sum.less.m.div2.ge0 <- sib.sum.less.m / 2 >= 0
  sib.sum.less.m.div2.lef <- sib.sum.less.m / 2 <= parp.cn
  ibd1p.conc <- sib.sum.less.m.even & sib.sum.less.m.div2.ge0 & sib.sum.less.m.div2.lef
  
  # IBD1M
  sib.sum.less.p <- sib.sum - parp.cn
  sib.sum.less.p.even <- sib.sum.less.p %% 2 == 0
  sib.sum.less.p.div2.ge0 <- sib.sum.less.p / 2 >= 0
  sib.sum.less.p.div2.lem <- sib.sum.less.p / 2 <= parm.cn
  ibd1m.conc <- sib.sum.less.p.even & sib.sum.less.p.div2.ge0 & sib.sum.less.p.div2.lem
  
  site.concordance <- ifelse(ibd.state=='IBD0', sums.eq, ifelse(ibd.state=='IBD2', sib.eq, ifelse(ibd.state=='IBD1P', ibd1p.conc, ifelse(ibd.state=='IBD1M', ibd1m.conc, NA))))
}

# Perform concordance analysis on each cluster within each family
conc.fam.bp <- ddply(bp, .(family, cluster), function(df) {
  
  conc.L <- concordance.stats(df, 'cn.L')
  conc.R <- concordance.stats(df, 'cn.R')
  
  if (length(conc.L)==0 || is.na(conc.L)) {
    return(data.frame())
  } else {
    return(data.frame(conc.L=conc.L, conc.R=conc.R, conc=conc.L&conc.R))
  }
})

# Now what? Needs work here. How do you filter conc.fam.bp? What are the covariates?
# We don't have the segment lengths. Perhaps I can filter on drops in CN?

