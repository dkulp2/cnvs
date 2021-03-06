# mk_quartet_segments.R - merge the predictions into new data structures for subsequent analysis
#
# Author: David Kulp, dkulp@broadinstitute.org
#
# Usage: Called using rscript:
#
# #1: top-level dir with subdirectories for each quartet (concatenate {A,B,C,D})
# #2: filename base (may contain leading directory after quartet letter)
# #3: pedigree file
# #4: IBD file
# #5: method - 'bayescsm'
#
# Write a 'segments.Rdata' file in the top-level dir that contains structures:
#   cnvs.sibs - a list of data.frames of sibs
#   cnvs.parents - a list of data.frames of parents (0: father, 1: mother)
#   quartets - data.frame (sib1, sib2, father, mother, sib1.gender, sib2.gender)
#   segs - a row per for each segment with consistent IBD and CN for all family members plus lots of labels
#   segs.group - one or more groups of segs with common CN=2 regions reduced to 5 bins
#   segs.all - a single data.frame with common CN=2 regions reduced to 5 bins

library(plyr)
library(dplyr)
library(RPostgreSQL)

profile.segment.cache.fn <- paste0(Sys.getenv('TMPDIR'),"/profile_segments.Rdata")

cmd.args <- commandArgs(trailingOnly = TRUE)
# Sys.setenv(PGHOST="localhost",PGUSER="dkulp",PGDATABASE="seq", PGOPTIONS="--search_path=data_sfari_batch1c_19may2017")
# cmd.args <- unlist(strsplit('/home/dkulp/data/SFARI.19May2017/data /sites_cnv_segs.txt /home/dkulp/data/sample_pedigrees.ped /home/dkulp/data/merged_ibd_regions.bed bayescsm',' '))
basedir <- cmd.args[1]
basename <- cmd.args[2]
ped.fn <- cmd.args[3]
ibd.bed.fn <- cmd.args[4]
method <- cmd.args[5]

if (file.exists(profile.segment.cache.fn)) {
    load(profile.segment.cache.fn)
    stopifnot(exists('psegs'))
} else {
    db <- src_postgres()
    message(Sys.time(),": Connected to ",db)
    psegs <- tbl(db, 'profile_segment') %>% collect(n=Inf)
    colnames(psegs) <- c('bin','chr','start','end','elength','gc','gtotal')
    save(psegs, file=profile.segment.cache.fn)
}
psegs[1,'start'] <- 1


big2 <- 200  # any region greater than big2 nts that is CN=2 for all samples is replaced by a small region
big2.replacement <- 5 # the size of the new CN=2 region
MIN.CNV.LEN <- 12


sibs <- sprintf("%s%s%s", basedir, c('A','B'), basename)
parents <- sprintf("%s%s%s", basedir, c('C','D'), basename)

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
  load(sprintf("%s.%s.Rdata", d, method))
  
  new.cnvs <- 
    ddply(cn.segs.merged, .(.id, chr), function(df) {
      sample.id <- first(df$.id)
      chr <- first(df$chr)
      # df <- filter(cn.segs.merged, .id=='SSC05534')
      if (!is.null(df$idx)) { df <- arrange(df, idx) }

      while (any(df$end.bin < df$start.bin)) {
        reversed <- which(df$end.bin < df$start.bin)
        
        message(Sys.time(),sprintf(": FIXME: Invalidating %s bogus predictions where END < START (%s,%s) and adjusting flanking segments in %s.\n",length(reversed),df$start.bin[reversed], df$end.bin[reversed], sample.id))
        
        # A reversal of B in A-B-C results in A and C overlapping.
        # Set A's end to B's end. Set C's start to B's start. Swap B and set CN=NA
        df[reversed-1,'end.map'] <- df[reversed,'end.map']
        df[reversed+1,'start.map'] <- df[reversed,'start.map']
        df[reversed,'start.map'] <- df[reversed-1,'end.map']
        df[reversed,'end.map'] <- df[reversed+1,'start.map']
        
        df[reversed-1,'end.bin'] <- df[reversed,'end.bin']
        df[reversed+1,'start.bin'] <- df[reversed,'start.bin']
        df[reversed,'start.bin'] <- df[reversed-1,'end.bin']
        df[reversed,'end.bin'] <- df[reversed+1,'start.bin']
        
        df[reversed,'cn'] <- NA
      }
      
      # removing single bin segments
      zlen <- which(df$start.bin==df$end.bin)
      if (length(zlen)>0) {
        message(Sys.time(),sprintf(": Removing %s zero length segments at %s", length(zlen), df$start.bin[zlen]))
        df <- df[-zlen,]
      }
      
      row.range <- 2:nrow(df)
      ovlps <- which(df$start.map[row.range] < df$end.map[row.range-1])+1
      if (length(ovlps) > 0) {
        message(Sys.time(),sprintf(": FIXME: Invalidating %s overlapping predictions %s (%s,%s) in %s\n", length(ovlps), ovlps, df$start.map[ovlps], df$end.map[ovlps-1], sample.id))
        df[ovlps,'cn'] <- NA
        df[ovlps-1,'cn'] <- NA
      }

      while (any(df$start.map[row.range]==df$end.map[row.range-1] &
             df$cn[row.range]==df$cn[row.range-1],na.rm=TRUE)) {  
        # sometimes the same CN has adjacent segments.
        adjacents <- which(df$start.map[row.range]==df$end.map[row.range-1] &
                             df$cn[row.range]==df$cn[row.range-1]) + 1
        
        if (length(adjacents) > 0) {
          if (any(adjacents!=nrow(df))) {
            # only bark if not the last
            message(Sys.time(), sprintf(": FIXME: Adjacent calls for %s with same CN at (%s,%s)\n", sample.id, adjacents-1,adjacents))
          }
          df$end.CI.L[adjacents-1] <- df$end.CI.L[adjacents]
          df$end.CI.R[adjacents-1] <- df$end.CI.R[adjacents]
          df$end.map[adjacents-1] <- df$end.map[adjacents]
          df$end.binCI.L[adjacents-1] <- df$end.binCI.L[adjacents]
          df$end.binCI.R[adjacents-1] <- df$end.binCI.R[adjacents]
          df$end.bin[adjacents-1] <- df$end.bin[adjacents]
          df <- df[-adjacents,]
        }
      }

      # Mark small CNVs as NA
      too.short <- which(df$end.bin-df$start.bin < MIN.CNV.LEN)
      df[too.short,'cn'] <- NA
      
      df
      
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
# fam <- filter(quartets, family=='11010')
ddply(quartets, .(family), function(fam) {
  ibd.fam <- filter(ibd, PAIR_NAME==paste(fam$sib1, fam$sib2, sep='-'))
  if (nrow(ibd.fam)==0) {  # SIB1 and SIB2 could be swapped with respect to ped file?
    ibd.fam <- filter(ibd, PAIR_NAME==paste(fam$sib2, fam$sib1, sep='-'))  
  }
  stopifnot(nrow(ibd.fam)>0)
  
  # iterate each chromosome
  ldply(unique(ibd.fam$CHR), function(chr) {
    # move pos forward to the minimum next end position and write a row with the current state
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
#        print(c(next.pos <=pos, ibd.fam.chr$END[ibd.idx], sib1$end.map[sib1.idx], sib2$end.map[sib2.idx], par1$end.map[par1.idx], par2$end.map[par2.idx]))
        if (next.pos <= pos) { 
         print(c(ibd.fam.chr$END[ibd.idx], sib1$end.map[sib1.idx], sib2$end.map[sib2.idx],
                 par1$end.map[par1.idx], par2$end.map[par2.idx]))
        }

#        print(c(ibd.idx, sib1.idx, sib2.idx, par1.idx, par2.idx))
        if (!is.na(next.pos)) {
          cat(paste0(paste(fam$family, fam$sib1, fam$sib2, fam$mother, fam$father, chr, pos, next.pos, sib1$cn[sib1.idx], sib2$cn[sib2.idx], par1$cn[par1.idx], par2$cn[par2.idx],
                           ibd.code(ibd.fam.chr[ibd.idx,]), sep=','),"\n"), file=t1.conn)
          
          if (ibd.idx <= nrow(ibd.fam.chr) && next.pos >= ibd.fam.chr$END[ibd.idx]) { ibd.idx <- ibd.idx + 1 }
          if (sib1.idx <= nrow(sib1) && next.pos >= sib1$end.map[sib1.idx]) { sib1.idx <- sib1.idx + 1 }
          if (sib2.idx <= nrow(sib2) && next.pos >= sib2$end.map[sib2.idx]) { sib2.idx <- sib2.idx + 1 }
          if (par1.idx <= nrow(par1) && next.pos >= par1$end.map[par1.idx]) { par1.idx <- par1.idx + 1 }
          if (par2.idx <= nrow(par2) && next.pos >= par2$end.map[par2.idx]) { par2.idx <- par2.idx + 1 }
        }        

#        print(c(ibd.idx, sib1.idx, sib2.idx, par1.idx, par2.idx))
#        print(c(pos,next.pos))
        pos <- next.pos
      }
    } else {
      stopifnot(nrow(sib1)==0 && nrow(sib2)==0)
    }
    

  })
})

close(t1.conn)
segs <- read.csv(t1.fn, as.is=TRUE, check.names=FALSE, header=FALSE)
colnames(segs) <- c('family','sib1','sib2','mother','father','chr','start','end','sib1.cn','sib2.cn','par1.cn','par2.cn','ibd.state')

# Label the segments with lots of state information
segs$cn.na <- is.na(segs$sib1.cn) | is.na(segs$sib2.cn) | is.na(segs$par1.cn) | is.na(segs$par2.cn)
segs$all.cn.na <- is.na(segs$sib1.cn) & is.na(segs$sib2.cn) & is.na(segs$par1.cn) & is.na(segs$par2.cn)
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



# assign bins. The segmenting is done in base space (partly historically and partly because the IBD is in base space)
# So this reassigns bins even though the input has bins. Cross fingers that it matches up!
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
                                                                    all.na=segs.chr.fam[[i]]$all.cn.na[fam.row[i]],
                                                                    cn2=(segs.chr.fam[[i]]$all.eq[fam.row[i]] & segs.chr.fam[[i]]$sib1.cn[fam.row[i]]==2))) })
          min.end <- min(rows$ends)
          max.start <- max(rows$starts)
          match.min.end <- rows$ends==min.end
          
          if (!is.na(all(rows$cn2)) & min.end - pos > big2 & (all(rows$cn2) || all(rows$all.na))) {
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

message(Sys.time(),sprintf(": Writing objects to %s/segments.Rdata", dirname(basedir)))
save(cnvs.sibs, cnvs.parents, quartets, segs, segs.group, segs.all, file=sprintf("%s/segments.Rdata", dirname(basedir)))

