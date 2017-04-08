# cnvs

dkulp@broadinstitute.org

conf - configuration files
eval - plots and stats
lib - R libraries
pl - pipeline for predicting cnvs
prep - data prep, transformation, filtering, etc.
ui - shiny app for displaying CNVS
util - shared scripts

======

To run on new dataset:

Edit cnvs.conf
time sh prep/prep.sh ; time sh pl/cnv_seg.sh 
