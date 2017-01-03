# bounds.R - infer site boundaries from multiple sample
#
# Author: David Kulp, dkulp@broadinstitute.org
#
# Usage: Called using rscript:
#
# #1: filename base - e.g. sites_cnv_segs.txt
# #2: method - which set of segments to read? ('smlcsm'). File read is base.method.Rdata
# #3: collapsed - which set of collapsed segments to read? ('smlxcsm'). File read is base.collapsed.Rdata
# #3: db connection - user:host:port:dbname - profile data is read from this connection

library(plyr)
library(dplyr)
library(RPostgreSQL)
library(zoo)

#this.dir <- dirname(substring(commandArgs()[grep('--file=', commandArgs())],8))
#if (length(this.dir)==0) { this.dir <- getwd() }

cmd.args <- commandArgs(trailingOnly = TRUE)
#cmd.args <- c('d:/mccarroll/cnv_seg.12.500/sites_cnv_segs.txt','smlcsm','smlxcsm','postgres:localhost:5432:seq',600)
cmd.args <- c('C:\\cygwin64\\home\\dkulp\\data\\out\\cnv_seg.B12.L500.Q13.4\\sites_cnv_segs.txt','smlcsm','smlxcsm','dkulp:localhost:5432:seq',600)
cnv.seg.fn <- cmd.args[1]
cnv.seg.method <- cmd.args[2]
cnv.seg.collapsed <- cmd.args[3]
db.conn.str <- cmd.args[4]
bp.dist <- cmd.args[5]

# Load the cnv prediction
load(sprintf("%s.%s.Rdata",cnv.seg.fn,cnv.seg.method)) # => cn.segs.merged
csm <- as.tbl(cn.segs.merged)
csm$len <- csm$seg <- csm$len2 <- csm$gap <- NULL
csm$len <- csm$end.map - csm$start.map 

# Load the collapsed cnv prediction
load(sprintf("%s.%s.Rdata",cnv.seg.fn,cnv.seg.collapsed)) # => cn.segs.merged
csmx <- as.tbl(cn.segs.merged)
csmx$len <- csmx$seg <- csmx$len2 <- csmx$gap <- NULL
csmx$len <- csmx$end.map - csmx$start.map 

# Connect to DB
db.conn.params <- as.list(unlist(strsplit(db.conn.str,":")))
names(db.conn.params) <- c('user','host','port','dbname')
db <- do.call(src_postgres, db.conn.params)

# Convert segments to transitions
# trans == [ sample, chr, pos, CN_left, CN_right ]

# Make sure all segments are contiguous within a sample and chr
stopifnot(all(ddply(csm, .(.id, chr), function(df) all(df$start.map[2:nrow(df)] == df$end.map[1:(nrow(df)-1)]))$V1))

trans <- ddply(csm, .(.id, chr), function(df) {
  data.frame(pos=df$end.map[1:nrow(df)-1], CNL=df$cn[1:nrow(df)-1], CNR=df$cn[2:nrow(df)])
})

# cluster transitions
trans <- ddply(trans, .(chr), function(df) {
  df <- arrange(df, pos)
  i1n <- 1:(nrow(df)-1)
  df$close <- c(df$pos[i1n+1]-df$pos[i1n] < bp.dist, FALSE)
  df$end.cluster <- 1:nrow(df)
  again <- TRUE
  while (again) {
    df$end.cluster2 <- c(ifelse(df$close[i1n], pmax(df$end.cluster[i1n],df$end.cluster[i1n+1]), df$end.cluster[i1n]),df$end.cluster[nrow(df)])
    again <- !all(df$end.cluster==df$end.cluster2)
    df$end.cluster <- df$end.cluster2
  }
  df$end.cluster2 <- NULL
  
  # hack: turn end.cluster into a factor, then assign its ordinal value as the cluster ID
  df$cluster.id <- as.numeric(as.factor(df$end.cluster))
  df$end.cluster <- df$close <- NULL
  return(df)
})


????
