# create a table containing the samples in each family quartet

library(RPostgreSQL)
library(plyr)
library(dplyr)
Sys.setenv(PGOPTIONS="--search_path=data_sfari_batch1")
db <- src_postgres()

ped.fn <- '/humgen/cnp04/bobh/projects/locus_SMN/cohort_sfari/sample_pedigrees.ped'
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

copy_to(db, quartets, temporary=FALSE)
