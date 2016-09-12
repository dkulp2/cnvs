library(ggplot2)
library(plyr)
library(dplyr)
cmd.args <- commandArgs(trailingOnly = TRUE)
bins <- 20
cmd.args[1] <- sprintf('d://mccarroll/cnv_seg.%s.500.90/breakpoints.txt',bins)
cmd.args[2] <- sprintf('d://mccarroll/cnv_seg.%s.500/breakpoints_oe.txt',bins)

bp <- read.table(cmd.args[1],header=FALSE,as.is=TRUE)
colnames(bp) <- c('site','sample', 'dir','offset', 'cn','geno','cnq','freq')
bp$cnq <- as.integer(bp$cnq)

# use only HET DELs
bp <- subset(bp, cn==1)
bp$freq.range <- cut_number(bp$freq, n=5)
bp$freq.range <- factor(ifelse(bp$freq == 1,'[1]',ifelse(bp$freq==2,'[2]',as.character(bp$freq.range))),
                        levels=c('[1]','[2]',levels(bp$freq.range)))

ggplot(subset(bp,abs(offset)<20), aes(x=as.factor(offset), y=cnq)) + geom_boxplot() + facet_grid(freq.range~dir) + xlab("offset") + ggtitle(sprintf("HET DEL CNQ Distribution w/rt Deletion Bounds\nStratified By CNV Frequency, %s Bin Windows",bins))

bpoe <- read.table(cmd.args[2],header=FALSE,as.is=TRUE, fill=TRUE)
colnames(bpoe) <- c('site','sample', 'dir','offset', 'cn', 'obs', 'exp', 'freq')
bpoe$obs <- as.integer(bpoe$obs)
bpoe$obsexp <- bpoe$obs / bpoe$exp

# use only HET DELs
bpoe <- subset(bpoe, cn==1)
bpoe$freq.range <- cut_number(bpoe$freq, n=5)

ggplot(subset(bpoe,abs(offset)<10), aes(x=as.factor(offset), y=obsexp)) + geom_boxplot() + facet_grid(.~dir) + xlab("offset") + ggtitle("HET DEL Obs/Exp Distribution w/rt Deletion Bounds") + ylim(0,3)
