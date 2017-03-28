library(dplyr)
library(RPostgreSQL)
db <- src_postgres()
load("sites_cnv_segs.txt.bayescsm.Rdata")
dbWriteTable(db$con, "cnv_post", cn.segs.merged)
dbGetQuery(db$con, "CREATE INDEX on cnv_post(chr, \"start.map\", \"end.map\", \".id\")");

load("sites_cnv_segs.txt.csm.Rdata")
dbWriteTable(db$con, "cnv_mrg", select(cn.segs.merged, .id, seg, chr, start.map, end.map, cn, copy.number))
dbGetQuery(db$con, "CREATE INDEX on cnv_mrg(chr, \"start.map\", \"end.map\", \".id\")")
