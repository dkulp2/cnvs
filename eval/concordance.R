# concordance.R
#
# compare the predictions among an IBD quartet

library(plyr)
library(dplyr)
library(ggplot2)

big2 <- 200  # any region greater than big2 nts that is CN=2 for all samples is replaced
big2.replacement <- 5 # the size of the new CN=2 region
MIN.CNV.LEN <- 1200

# load data ########################
basedir <- 'C:\\cygwin64\\home\\dkulp\\data\\SFARI'

sibs <- c('data_sfari_batch1a','data_sfari_batch1b')
parents <- c('data_sfari_batch1c','data_sfari_batch1d')

ibd.bed.fn <- sprintf("%s/merged_ibd_regions.bed", basedir)
ped.fn <- sprintf("%s/sample_pedigrees.ped", basedir)
cnv.fn <- 'sites_cnv_segs.txt.bayescsm.Rdata'

ibd <- read.table(ibd.bed.fn, header=TRUE, stringsAsFactors = FALSE)
ibd <- filter(ibd, CHR=='20')  # FIXME: remove later

# returns the highest mapped base
chrom.max <- function(chr) {
  max(ibd$END[ibd$CHR==chr])
}

# fill in the gaps with STATE=UNK
ibd.gaps <- ddply(ibd, .(PAIR_NAME, SIB1, SIB2, CHR), function(df) {
  df <- arrange(df, START)
  row.range <- 2:nrow(df)
  gap.pre <- which(df$START[row.range] > df$END[row.range-1])+1
  if (length(gap.pre)>0) {
    lbound <- df$END[gap.pre-1]
    rbound <- df$START[gap.pre]
    between <- data.frame(START=lbound, END=rbound, STATE='UNK', PATERNAL_IND=NA, MATERNAL_IND=NA, SITES=0, SITE_AGREEMENTS=0)
  } else {
    between <- data.frame()
  }

  rbind(data.frame(START=1, END=df$START[1], STATE='UNK', PATERNAL_IND=NA, MATERNAL_IND=NA, SITES=0, SITE_AGREEMENTS=0),
        between,
        data.frame(START=df$END[nrow(df)], END=chrom.max(df$CHR[1]), STATE='UNK', PATERNAL_IND=NA, MATERNAL_IND=NA, SITES=0, SITE_AGREEMENTS=0))
})

ibd <- arrange(rbind(ibd, ibd.gaps), PAIR_NAME, CHR, START)

ped <- read.table(ped.fn, col.names=c('family','sample','parent1','parent2','gender','casecontrol'), stringsAsFactors = FALSE)

gender <- function(num) { as.factor(ifelse(num==1,'M','F')) }
quartets <- ddply(ped, .(family), function(df) { 
  children <- filter(df, parent1!=0)
  child1 <- children[1,]
  child2 <- children[2,]
  father <- filter(df, parent1==0 & gender==1)
  mother <- filter(df, parent1==0 & gender==2)
  
  if (nrow(children)==2 && nrow(father)==1 && nrow(mother)==1) {
    data.frame(sib1=child1$sample, sib2=child2$sample, father=father$sample, mother=mother$sample, 
               sib1.gender=gender(child1$gender), sib2.gender=gender(child2$gender), stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
})


# load the "cn.segs.merged" table and add CN=2 for all gaps
load.cnvs <- function(d) { 
  load(sprintf("%s/%s/%s", basedir, d, cnv.fn))
  if (any(cn.segs.merged$end.map <= cn.segs.merged$start.map)) {
    cat("FIXME: Removing bogus predictions where END < START.\n")
    cn.segs.merged <- filter(cn.segs.merged, end.map > start.map)
  }
  new.cnvs <- 
    ddply(cn.segs.merged, .(.id, chr), function(df) {
      df <- arrange(filter(df, end.map-start.map > MIN.CNV.LEN), start.map)
      row.range <- 2:nrow(df)
      gap.pre <- which(df$start.map[row.range] > df$end.map[row.range-1])+1
      if (length(gap.pre)>0) {
        lbound <- df$end.map[gap.pre-1]
        rbound <- df$start.map[gap.pre]
        between <- data.frame(label=paste0(df$label[gap.pre],'G'), 
                              start.CI.L=lbound, start.map=lbound, start.CI.R=lbound,
                              end.CI.L=rbound, end.map=rbound, end.CI.R=rbound, cn=2)
      } else { between <- data.frame() }
      
      ovlps <- which(df$start.map[row.range] < df$end.map[row.range-1])+1
      if (length(ovlps) > 0) {
        cat(sprintf("FIXME: Removing %s overlapping predictions in %s\n", length(ovlps), first(df$.id)))
        print(df[ovlps,])
        df <- df[-ovlps,]
      }
      
      rbind(data.frame(label=paste0(df$.id[1],'_START'),
                       start.CI.L=1, start.map=1, start.CI.R=1,
                       end.CI.L=df$start.map[1], end.map=df$start.map[1], end.CI.R=df$start.map[1], cn=2),
            between,
            subset(df,select=-c(.id, chr)),
            data.frame(label=paste0(df$.id[1],'_END'),
                       start.CI.L=df$end.map[nrow(df)], start.map=df$end.map[nrow(df)], start.CI.R=df$end.map[nrow(df)],
                       end.CI.L=chrom.max(df$chr[1]), end.map=chrom.max(df$chr[1]), end.CI.R=chrom.max(df$chr[1]), cn=2))
    })
  return(arrange(new.cnvs, .id, chr, start.map))
}

cnvs.sibs <- lapply(sibs, load.cnvs)
cnvs.parents <- lapply(parents, load.cnvs)

# generate overlap states ###########
# increment a counter and move pointers along the rows of the quartet members
# generate rows of the form:
#
# chr start end [S1|S2|SB] [P1|P2|PB] [IBD0|IBD2|IBD1P|IBD1M]

# iterate over each quartet
t1.fn <- tempfile("segs")
t1.conn <- file(t1.fn,"w")

ibd.code <- function(df) { 
  ifelse(df$STATE=='UNK', NA, ifelse(df$STATE=='ONE', paste0('IBD1', ifelse(df$PATERNAL_IND>df$MATERNAL_IND,'P','M')), ifelse(df$STATE=='ZERO', 'IBD0', 'IBD2')))
}

# fam <- filter(quartets, mother=='SSC01112')
# fam <- filter(quartets, family=='11006')
ddply(quartets, .(family), function(fam) {
  ibd.fam <- filter(ibd, PAIR_NAME==paste(fam$sib1, fam$sib2, sep='-'))
  if (nrow(ibd.fam)==0) {  # SIB1 and SIB2 could be swapped with respect to ped file?
    ibd.fam <- filter(ibd, PAIR_NAME==paste(fam$sib2, fam$sib1, sep='-'))  
  }
  stopifnot(nrow(ibd.fam)>0)
  cat(fam$family,"\n")
  
  # iterate each chromosome
  ldply(unique(ibd.fam$CHR), function(chr) {
    # move pos forward to the minimum next end position and write a row with the current state
    cat(chr," ")
    pos <- 1 # base position in chrom
    
    sib1 <- arrange(filter(cnvs.sibs[[1]], .id==fam$sib1 & chr==chr), start.map)
    sib2 <- arrange(filter(cnvs.sibs[[2]], .id==fam$sib2 & chr==chr), start.map)
    if (nrow(sib1)>0 && nrow(sib2)>0) {
      par1 <- arrange(filter(cnvs.parents[[2]], .id==fam$father & chr==chr), start.map) # By convention it appears that father is "D" data set and mother is "C"
      par2 <- arrange(filter(cnvs.parents[[1]], .id==fam$mother & chr==chr), start.map)
      stopifnot(nrow(par1)>0 && nrow(par2)>0)
      
      ibd.fam.chr <- arrange(filter(ibd.fam, CHR==chr), START)
      sib1.idx <- sib2.idx <- par1.idx <- par2.idx <- ibd.idx <- 1
      
      while (!is.na(pos) && pos < max(ibd.fam.chr$END)) {
        next.pos <- min(ibd.fam.chr$END[ibd.idx], sib1$end.map[sib1.idx], sib2$end.map[sib2.idx], 
                        par1$end.map[par1.idx], par2$end.map[par2.idx], na.rm=TRUE)
        # print(c(ibd.fam.chr$END[ibd.idx], sib1$end.map[sib1.idx], sib2$end.map[sib2.idx], par1$end.map[par1.idx], par2$end.map[par2.idx]))
        if (next.pos <= pos) { browser() }

        # print(c(ibd.idx, sib1.idx, sib2.idx, par1.idx, par2.idx))
        if (!is.na(next.pos)) {
          cat(paste0(paste(fam$family, fam$sib1, fam$sib2, fam$mother, fam$father, chr, pos, next.pos, sib1$cn[sib1.idx], sib2$cn[sib2.idx], par1$cn[par1.idx], par2$cn[par2.idx],
                           ibd.code(ibd.fam.chr[ibd.idx,]), sep=','),"\n"), file=t1.conn)
          
          if (ibd.idx <= nrow(ibd.fam.chr) && next.pos >= ibd.fam.chr$END[ibd.idx]) { ibd.idx <- ibd.idx + 1 }
          if (sib1.idx <= nrow(sib1) && next.pos >= sib1$end.map[sib1.idx]) { sib1.idx <- sib1.idx + 1 }
          if (sib2.idx <= nrow(sib2) && next.pos >= sib2$end.map[sib2.idx]) { sib2.idx <- sib2.idx + 1 }
          if (par1.idx <= nrow(par1) && next.pos >= par1$end.map[par1.idx]) { par1.idx <- par1.idx + 1 }
          if (par2.idx <= nrow(par2) && next.pos >= par2$end.map[par2.idx]) { par2.idx <- par2.idx + 1 }
        }        

        # print(c(ibd.idx, sib1.idx, sib2.idx, par1.idx, par2.idx))
        # print(c(pos,next.pos))
        pos <- next.pos
      }
    } else {
      stopifnot(nrow(sib1)==0 && nrow(sib2)==0)
    }
    

  })
  cat("\n")
})

close(t1.conn)
segs <- read.csv(t1.fn, as.is=TRUE, check.names=FALSE, header=FALSE)
colnames(segs) <- c('family','sib1','sib2','mother','father','chr','start','end','sib1.cn','sib2.cn','par1.cn','par2.cn','ibd.state')

segs$sib.eq <- (segs$sib1.cn == segs$sib2.cn)
segs$all.eq <- (segs$sib.eq & segs$sib1.cn==segs$par1.cn & segs$sib1.cn==segs$par2.cn)
segs$sib.del <- segs$sib1.cn < 2 | segs$sib2.cn < 2
segs$sib.wt <- segs$sib.eq & segs$sib1.cn == 2
segs$all.wt <- segs$all.eq & segs$sib1.cn == 2
segs$sib.sum <- segs$sib1.cn + segs$sib2.cn
segs$par.sum <- segs$par1.cn + segs$par2.cn
segs$sums.eq <- segs$sib.sum == segs$par.sum

# Haploid CN is the CN on each chromosome. Diploid CN is the unphased CN across both samples chromsomes.
# Each Diploid CN has CN+1 possible Haploid CN pairs, e.g. dCN=4 then could be hCN=4+0, 3+1, 2+2, 1+3, 0+4.
# Phasing of hCN isn't known, so there are ceil((CN+1)/2) observable hCN pairs (e.g. 3+1 is the same as 1+3).
# If IBD1M, then sibs share same maternal chromosome, but different paternal.
# That is, the siblings share the same hCNm, but different hCNp.
# For example, suppose dCNm=2 and dCNp=3.
# One of the maternal chromosomes takes values between hCNm=0..dCNm, shared by both.
# The paternal takes two values, a=0..dCNp, b=dCNp-a
# Sib1: hCNm + a 
# Sib2: hCNm + b
#
# So each cell is sibling CN [CN1,CN2]
#
#       a+b
#     0,3           1,2          2,1           3,0
# 0  [(0+0),(0+3)] [(0+1),(0+2)] [(0+2),(0+1)] [(0+3),(0+0)]
# 1  [(1+0),(1+3)] [(1+1),(1+2)] [(1+2),(1+1)] [(1+3),(1+0)]
# 2  [(2+0),(2+3)] [(2+1),(2+2)] [(2+2),(2+1)] [(2+3),(2+0)]
#
#       a+b
#     0,3           1,2          2,1           3,0
# 0  [0,3]          [1,2]        [2,1]         [3,0]
# 1  [1,4]          [2,3]        [3,2]         [4,1]
# 2  [2,5]          [3,4]        [4,3]         [5,2]
#
# So any of these pairs of [CN1,CN2] are concordant with dCNm=2 and dCNp=3
#
# The sum of CN1+CN2 must be 0+3, 2+3, or 4+3 = 3,5,7 => 2*(0..dCNm)+dCNp
# but this is not a sufficient constraint. E.g. [1,6] is not allowed.
# Increasing dCNm to 3 you can see the pattern that the pairs with large spread are not allowed.
# This is because the parental CN must be split across sibs.
#
#       a+b
#     0,3           1,2          2,1           3,0
# 0  [0,3]          [1,2]        [2,1]         [3,0]
# 1  [1,4]          [2,3]        [3,2]         [4,1]
# 2  [2,5]          [3,4]        [4,3]         [5,2]
# 3  [3,6]          [4,5]        [5,4]         [6,3]
#
# Thus, there is an additional constraint of abs(CN1-CN2) <= dCNp 

# ibd1P requires sum of sibs to be 2*(0..dCNp)+dCNm. 
# How to check in parallel when the value of dCNp is different for each row in segs?
# The constraint is CN1+CN2 = 2*(0..dCNp)+dCNm.
# So (CN1+CN2)-dCNm = 2*(0..dCNp)
#   ((CN1+CN2)-dCNm)/2 = 0..dCNp
# Then
#   ((CN1+CN2)-dCNm)%2 = 0 and ((CN1+CN2)-dCNm)/2 <= dCNp and ((CN1+CN2)-dCNm)/2 >= 0

# IBD1P
segs$sib.sum.less.m <- segs$sib.sum - segs$par2.cn
segs$sib.sum.less.m.even <- segs$sib.sum.less.m %% 2 == 0
segs$sib.sum.less.m.div2.ge0 <- segs$sib.sum.less.m / 2 >= 0
segs$sib.sum.less.m.div2.lef <- segs$sib.sum.less.m / 2 <= segs$par1.cn
segs$ibd1p.conc <- segs$sib.sum.less.m.even & segs$sib.sum.less.m.div2.ge0 & segs$sib.sum.less.m.div2.lef

# IBD1M
segs$sib.sum.less.p <- segs$sib.sum - segs$par1.cn
segs$sib.sum.less.p.even <- segs$sib.sum.less.p %% 2 == 0
segs$sib.sum.less.p.div2.ge0 <- segs$sib.sum.less.p / 2 >= 0
segs$sib.sum.less.p.div2.lem <- segs$sib.sum.less.p / 2 <= segs$par2.cn
segs$ibd1m.conc <- segs$sib.sum.less.p.even & segs$sib.sum.less.p.div2.ge0 & segs$sib.sum.less.p.div2.lem

segs$concordant <- ifelse(segs$ibd.state=='IBD0', segs$sums.eq, ifelse(segs$ibd.state=='IBD2', segs$sib.eq, ifelse(segs$ibd.state=='IBD1P', segs$ibd1p.conc, segs$ibd1m.conc)))

segs$cnv <- ifelse(segs$sib.del, 'DEL', ifelse(segs$sib.wt, 'WT', 'DUP'))
segs$fam <- as.factor(segs$family)

segs$len <- segs$end - segs$start
segs <- arrange(segs, family, chr, start)

# display density of CNVs
lapply(unique(segs$chr), function(chr) {
  z <- unlist(dlply(filter(segs, !all.wt & chr==chr), .(family, chr, start), function(r) seq(r$start,r$end)), use.names = FALSE)
  print(plot(density(z), main=sprintf("Density of CNVs along chr %s", chr)))
})
# display location of CNVs. Only the largest CNVs are visible
print(ggplot(filter(segs, !all.eq), aes(x=start, xend=end, y=fam, yend=fam, color=all.eq)) + geom_segment(size=10) + theme_bw() + geom_blank() + theme(panel.grid.major.y=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor.y = element_blank(), legend.position="none"))

# scan through all segments simultaneously
# If we are all overlapping a CN=2 and it is big, then reduce its length for all

group.sz <- 20

segs.group <-
lapply(seq(1,100,group.sz), function(x) {
  all.cn2 <- list()  # each element is a vector of row indices, one for each sample
  remove.nt <- list() # amount to remove

  segs2 <-
    ddply(filter(segs, family %in% quartets$family[x:(x+group.sz-1)]), .(chr), function(segs.chr) {
      segs.chr.fam <- dlply(segs.chr, .(family), identity)
      num.fams <- length(segs.chr.fam)
      fam.row <- rep(1, num.fams)
      pos <- 1
      while (!is.na(pos) && pos < max(segs.chr$end)) {
        rows <- ldply(1:num.fams, function(i) { return(data.frame(family=names(segs.chr.fam)[i],
                                                                  ends=segs.chr.fam[[i]]$end[fam.row[i]], 
                                                                  starts=segs.chr.fam[[i]]$start[fam.row[i]], 
                                                                  len=segs.chr.fam[[i]]$len[fam.row[i]],
                                                                  cn2=(segs.chr.fam[[i]]$all.eq[fam.row[i]] & segs.chr.fam[[i]]$sib1.cn[fam.row[i]]==2))) })
        min.end <- min(rows$ends)
        max.start <- max(rows$starts)
        match.min.end <- rows$ends==min.end
        
        if (min.end - pos > big2 & all(rows$cn2)) {
          # push the row indices for all samples
          all.cn2[[length(all.cn2)+1]] <- fam.row
          
          # push the number of bases to remove 
          remove.nt[[length(remove.nt)+1]] <- (min.end-pos)-big2.replacement
        }
        
        # advance the row index for those families that match the right edge of the current overlap region
        fam.row[match.min.end] <- fam.row[match.min.end]+1
        
        # advance pointer
        pos <- min.end
      }
      
      
      # return a new set of data.frames for this chromosome after reducing the length of rows in all.cn2
      # I subtract by the specified amount (instead of setting to a fixed length) because the same quartet segment
      # might be part of more than one consecutive all.cn2 region because of a change in ibd.state for one of the quartets.
      # So the length of the all.cn2 regions may be a multiple o big2.replacement.
      new.segs <- ldply(1:num.fams, function(f) {
        cat(names(segs.chr.fam)[f],' ')
        sapply(1:length(all.cn2), function(idx) {
          segs.chr.fam[[f]][all.cn2[[idx]][f],'len'] <<- segs.chr.fam[[f]][all.cn2[[idx]][f],'len'] - remove.nt[[idx]]
        })
        segs.chr.fam[[f]]
      })
      cat("\n")
      
      # iterate through new segments and adjust start and end according to new length
      final.segs <- ddply(new.segs, .(family), function(df) {
        df$end <- cumsum(df$len) - 1
        df$start <- c(1, df$end[1:(nrow(df)-1)])
        df
      })
      
      return(final.segs)  
    })
  
  # # adjust length of CN=2 to a small value
  # # use rolling sum to adjust all rows
  # segs <- ddply(segs, .(family, chr), function(df) {
  #   big.twos <- df$sib1.cn==df$sib2.cn & df$sib1.cn==df$par1.cn & df$sib1.cn==df$par2.cn & df$sib1.cn==2 & df$len > 5000
  #   df$len[big.twos] <- 50
  #   df$end <- 1+cumsum(df$len)
  #   df$start[2:nrow(df)] <- df$end[1:(nrow(df)-1)]
  #   return(df)
  # })
  
#  ggplot(filter(segs2, ibd.state=='IBD2' & len>1200 & family %in% quartets$family[1:10]), aes(x=start, xend=end, y=sib1.cn, yend=sib1.cn, color=sib.eq)) + geom_segment(size=2) + facet_grid(family~.)
  
#  print(ggplot(filter(segs2, !sib.wt & len>1200), aes(x=start, xend=end, y=fam, yend=fam, color=sib.eq)) + geom_segment(size=3) + facet_grid(ibd.state~cnv))
  
#  ggplot(filter(segs2, family %in% quartets$family[1:10]), aes(x=start, xend=end, y=sib1.cn, yend=sib1.cn, color=ibd.state)) + geom_point() + facet_grid(family~.)
  
})

pdf("concordance.pdf")
MIN.OVLP.LEN <- 400
sapply(segs.group, function(df) {
  print(ggplot(filter(df, !is.na(ibd.state) & !all.wt & len>MIN.OVLP.LEN), aes(x=start, xend=end, y=fam, yend=fam, color=concordant)) + 
          geom_segment(size=2) + geom_point() + facet_grid(ibd.state~cnv) + ggtitle(paste0('Samples ',df$family[1],'..',df$family[nrow(df)],"\nseg > ",MIN.OVLP.LEN)))
})
dev.off()

