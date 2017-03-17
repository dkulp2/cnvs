# Run once as
# Rscript setup.R > env.txt
#

source("../../conf/config.R")

vals <- c(do.conf("../../conf/site.conf"), do.conf("../../conf/cnv.conf"), do.conf("../../conf/shiny.conf"))
write.conf(stdout(), vals)

