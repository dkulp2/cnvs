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

psegs <- read.table("C:\\cygwin64\\home\\dkulp\\data\\tmp\\profile_segments.txt")[,1:4]
colnames(psegs) <- c('bin','chr','start','end')
psegs[1,'start'] <- 1

big2 <- 200  # any region greater than big2 nts that is CN=2 for all samples is replaced by a small region
big2.replacement <- 5 # the size of the new CN=2 region
MIN.CNV.LEN <- 12
#USE_BINS <- FALSE
USE_BINS <- TRUE

# load data ########################
# basedir <- 'C:\\cygwin64\\home\\dkulp\\data\\SFARI'
# 
# sibs <- c('data_sfari_batch1a','data_sfari_batch1b')
# parents <- c('data_sfari_batch1c','data_sfari_batch1d')

basedir <- 'C:\\cygwin64\\home\\dkulp\\data\\SFARI.12Mar2017_test'
basedir <- 'C:\\cygwin64\\home\\dkulp\\data\\SFARI.1Apr2017'
basedir <- 'C:\\cygwin64\\home\\dkulp\\data\\SFARI.11Apr2017'
basedir <- 'C:\\cygwin64\\home\\dkulp\\data\\SFARI.11Apr2017b'
basedir <- 'C:\\cygwin64\\home\\dkulp\\data\\SFARI.27April2017'

sibs <- c('dataA','dataB')
parents <- c('dataC','dataD')

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

# create a hash lookup: family.lookup[[sample]] => family
# Only good for one sample at a time.
# create a family.sample data.frame for joining
family.lookup <- new.env()
family.sample <- 
  ddply(quartets, .(family), function(f) { 
    assign(f$sib1, f$family, envir=family.lookup)
    assign(f$sib2, f$family, envir=family.lookup)
    assign(f$father, f$family, envir=family.lookup)
    assign(f$mother, f$family, envir=family.lookup)
    return(data.frame(sample=c(f$sib1, f$sib2, f$father, f$mother), who=c('sib1','sib2','father','mother')))
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
      sample.id <- first(df$.id)
      chr <- first(df$chr)
      cat(sample.id,"\n")
      if (USE_BINS) {
        df <- arrange(filter(df, end.bin-start.bin > MIN.CNV.LEN), start.map)
      } else {
        df <- arrange(filter(df, end.map-start.map > MIN.CNV.LEN*100), start.map)
      }
      if (nrow(df)<=1) {
        between <- data.frame()
      } else {
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
          row.range <- 2:nrow(df)
        }
        
        
        # sometimes the same CN has adjacent segments.
        adjacents <- which(df$start.map[row.range]==df$end.map[row.range-1] &
                             df$cn[row.range]==df$cn[row.range-1]) + 1
        #       if (length(adjacents) > 2) {
        #         adj.range <- 2:(length(adjacents)-1)
        #         tweens <- adjacents[adj.range]==adjacents[adj.range-1]+1 & adjacents[adj.range]==adjacents[adj.range+1]-1
        # #        stopifnot(all(!tweens)) # do a closure and remove tweens
        #       }
        if (length(adjacents) > 0) {
          cat(sprintf("FIXME: Adjacent calls for %s with same CN\n", first(df$.id)))
          df$end.CI.L[adjacents-1] <- df$end.CI.L[adjacents]
          df$end.CI.R[adjacents-1] <- df$end.CI.R[adjacents]
          df$end.map[adjacents-1] <- df$end.map[adjacents]
          df <- df[-adjacents,]
        }
      }
      
      if (nrow(df) > 0) {
        rbind(data.frame(label=paste0(sample.id,'_START'),
                         start.CI.L=1, start.map=1, start.CI.R=1,
                         end.CI.L=df$start.map[1], end.map=df$start.map[1], end.CI.R=df$start.map[1], cn=2),
              between,
              select(subset(df,select=-c(.id, chr)), label, start.CI.L, start.map, start.CI.R, end.CI.L, end.map, end.CI.R, cn),
              data.frame(label=paste0(sample.id,'_END'),
                         start.CI.L=df$end.map[nrow(df)], start.map=df$end.map[nrow(df)], start.CI.R=df$end.map[nrow(df)],
                         end.CI.L=chrom.max(df$chr[1]), end.map=chrom.max(df$chr[1]), end.CI.R=chrom.max(df$chr[1]), cn=2))
      } else {
        data.frame(label=paste0(sample.id,'_START_END'),
                   start.CI.L=1, start.map=1, start.CI.R=1,
                   end.CI.L=chrom.max(chr), end.map=chrom.max(chr), end.CI.R=chrom.max(chr), cn=2)
      }
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
      par2 <- arrange(filter(cnvs.parents[[2]], .id==fam$father & chr==chr), start.map) # By convention it appears that father is "D" data set and mother is "C"
      par1 <- arrange(filter(cnvs.parents[[1]], .id==fam$mother & chr==chr), start.map)
      if (nrow(par1)==0 || nrow(par2)==0) {
        warning(Sys.time(),sprintf(": No CNVs predicted for at least one parent. P=%s M=%s", fam$father, fam$mother))
      }
      
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
segs$par.del <- segs$par1.cn < 2 | segs$par2.cn < 2
segs$any.del <- segs$sib.del | segs$par.del
segs$sib.dup <- segs$sib1.cn > 2 | segs$sib2.cn > 2
segs$par.dup <- segs$par1.cn > 2 | segs$par2.cn > 2
segs$any.dup <- segs$sib.dup | segs$par.dup
segs$is.del <- segs$any.del & !segs$any.dup  # Bob's definition
segs$cn.gt.4 <- segs$sib1.cn > 4 | segs$sib2.cn > 4 | segs$par1.cn > 4 | segs$par2.cn > 4 # 
segs$bi.all <- !segs$any.del & segs$any.dup & !segs$cn.gt.4
   
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
segs$sib.sum.less.m <- segs$sib.sum - segs$par1.cn
segs$sib.sum.less.m.even <- segs$sib.sum.less.m %% 2 == 0
segs$sib.sum.less.m.div2.ge0 <- segs$sib.sum.less.m / 2 >= 0
segs$sib.sum.less.m.div2.lef <- segs$sib.sum.less.m / 2 <= segs$par2.cn
segs$ibd1p.conc <- segs$sib.sum.less.m.even & segs$sib.sum.less.m.div2.ge0 & segs$sib.sum.less.m.div2.lef

# IBD1M
segs$sib.sum.less.p <- segs$sib.sum - segs$par2.cn
segs$sib.sum.less.p.even <- segs$sib.sum.less.p %% 2 == 0
segs$sib.sum.less.p.div2.ge0 <- segs$sib.sum.less.p / 2 >= 0
segs$sib.sum.less.p.div2.lem <- segs$sib.sum.less.p / 2 <= segs$par1.cn
segs$ibd1m.conc <- segs$sib.sum.less.p.even & segs$sib.sum.less.p.div2.ge0 & segs$sib.sum.less.p.div2.lem

segs$concordant <- ifelse(segs$ibd.state=='IBD0', segs$sums.eq, ifelse(segs$ibd.state=='IBD2', segs$sib.eq, ifelse(segs$ibd.state=='IBD1P', segs$ibd1p.conc, segs$ibd1m.conc)))

segs$cnv <- ifelse(segs$is.del, 'DEL', ifelse(segs$bi.all, 'BI', ifelse(segs$all.wt, 'WT', 'MULTI')))
#segs$cnv <- ifelse(segs$sib.del, 'DEL', ifelse(segs$all.wt, 'WT', 'OTHER')) # OLD DEFN
segs$fam <- as.factor(segs$family)

segs$len <- segs$end - segs$start

# 3 regions that, by eye, have long, high frequency, discordant deletions 
excl.ranges <- list(c(29400000,29600000),c(41200000,41300000),c(60500000,60550000))

segs$excl <- with(segs,(start<29600000 & end>29400000)|(start<41300000 & end>41200000) | (start<60550000 & end>60500000))
segs$excl.alpha <- with(segs,(start<29600000 & end>29400000))


# assign bins
bin.map <- function(pos) {
  tmp.lines <- file()
  ps.idx <- 1
  idx <- 1
  
  while (idx <= length(pos)) {
    while (psegs$end[ps.idx] < pos[idx]) { ps.idx <- ps.idx + 1 }
    stopifnot(psegs$start[ps.idx] <= pos[idx])
    cat(psegs$bin[ps.idx], "\n", file=tmp.lines)
    idx <- idx + 1
  }
  
  m <- read.table(tmp.lines)
  close(tmp.lines)
  return(m$V1)
  
}

segs <- arrange(segs, chr, start)
segs$start.bin <- bin.map(segs$start)

segs <- arrange(segs, chr, end)
segs$end.bin <- bin.map(segs$end)

segs <- arrange(segs, family, chr, start)





# scan through all segments simultaneously
# If we are all overlapping a CN=2 and it is big, then reduce its length for all


mk.segs.group <- function(group.sz) {
  group.seq <- seq(1,length(unique(segs$family)),group.sz)
  lapply(group.seq, function(x) {
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
          df$old.start <- df$start
          df$old.end <- df$end
          df$end <- cumsum(df$len) - 1
          df$start <- c(1, df$end[1:(nrow(df)-1)])
          df
        })
        
        return(final.segs)  
      })
    
  })
  
} 
segs.group <- mk.segs.group(20)
segs.all <- mk.segs.group(100)

pdf("concordance.pdf")
MIN.OVLP.LEN <- 1400
sapply(segs.group, function(df) {
 print(ggplot(filter(df, !is.na(ibd.state) & !all.wt & len>MIN.OVLP.LEN & cnv!='WT'), aes(x=start, xend=end, y=fam, yend=fam, color=concordant)) +
     geom_segment(size=2) + geom_point() + facet_grid(ibd.state~cnv) + ggtitle(paste0('Samples ',df$family[1],'..',df$family[nrow(df)],"\ninterval > ",MIN.OVLP.LEN)))
})
dev.off()

sa1 <- filter(segs.all[[1]], !is.na(ibd.state) & cnv=='DEL')
print(ggplot(filter(sa1, len>MIN.OVLP.LEN & !concordant), aes(x=old.start, xend=old.end, y=fam, yend=fam, color=len)) + scale_colour_gradient(high = "red", low = "orange") +
        geom_segment(size=2) + geom_point() + ggtitle(sprintf("Discordant Sites\nInterval > %s nt",MIN.OVLP.LEN)))

# print(ggplot(filter(sa1, len>MIN.OVLP.LEN & !concordant), aes(x=old.start, xend=old.end, y=fam, yend=fam, color=len)) + scale_colour_gradient(high = "red", low = "orange") +
#         geom_segment(size=2) + geom_point() + ggtitle(sprintf("Discordant Sites\nInterval > %s nt",MIN.OVLP.LEN)) + xlm)

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
  c(.158e7,.16e7,1),
  c(.6296e7,.6297e7,2),
  c(2.575e7,2.581e7,3),
  c(2.949e7,2.951e7,4),
  c(4.124e7,4.125e7,5),
  c(6.05175e7,6.0519e7,6)
  
)

pdf("discordant_loci.pdf")
lapply(bad.loci, function(ab) {
  z4p(z4, ab[1],ab[2], ab[3])
})
dev.off()

# count concordance by base and cnv
len.max <- 15000
len.min.disp <- 1200
#steps <- c(seq(0,5000,100), seq(6000,len.max,1000))
steps <- c(seq(0,50,1), seq(60,len.max/100,10))

cts <- function(segs.subset) {
  ldply(steps, function(len.min) {
    z <- ddply(filter(segs.subset, !is.na(ibd.state) & end.bin-start.bin > len.min), .(cnv), function(df) {
#    z <- ddply(filter(segs.subset, !is.na(ibd.state) & len > len.min), .(cnv), function(df) {
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
             aes(x=len.min, y=conc.cnv.pct, color=grp, size=all))+ geom_line() + facet_grid(cnv~.) + xlab("Min Segment Size") + ylab("Concordance") + ggtitle("Quartet Concordance By CNV")  + geom_vline(xintercept = 12, linetype=2) + xlim(0,50)+ylim(0.8,1) + guides(size="none")

# generate AUC for each of the "hold one out" ROCs
auc.res <- ldply(rocs, function(df) {
  ddply(df, .(cnv,grp), summarize, auc=auc(len.min, conc.cnv.pct))
})

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

excl.order <- c(1,6,2,4,3)

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
print(p5 + xlim(c(5,50)) + coord_cartesian(ylim=c(0.8,1)))  # + ggtitle("Rank #5 first"))

write.table(select(segs, family, sib1, sib2, mother, father, chr, start, end, sib1.cn, sib2.cn, par1.cn, par2.cn, ibd.state, concordant, cnv, start.bin, end.bin), file="segs.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

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

