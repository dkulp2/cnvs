# Run once as
# Rscript setup.R schema > env.txt
#

library(RPostgreSQL)
library(dplyr)

source("../../conf/config.R")

cmd.args <- commandArgs(trailingOnly = TRUE)
schema <- cmd.args[1]

db <- src_postgres()
invisible(dbGetQuery(db$con, sprintf("set search_path=%s", schema)))

vals <- c(do.conf("../../conf/site.conf"), do.conf("../../conf/cnv.conf"), do.conf("../../conf/shiny.conf"))
write.conf(stdout(), vals)
store.conf(db, vals)

