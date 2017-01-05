#!/bin/bash -eu

set -eux

THISDIR=`dirname $0`
export PATH=${THISDIR}/../util:$PATH
echo $PATH

# set ROOT and primary data sources
source ${THISDIR}/../conf/site.conf

# set params for this run
source ${THISDIR}/../conf/cnv.conf

SCRIPTS=${THISDIR}/../pl
ROLLING_WINDOWS=${SCRIPTS}/choose_windows.awk

eval "gunzip -c ${profileFile}" | awk -v OFS="\t" -v NBINS=${NBINS} -v MAXLEN=${MAXLEN} -f ${ROLLING_WINDOWS} > ${ROOT}/sites.txt

zcat ${gsdelFile} | gs_del2tab -s ${ROOT}/sites.txt > ${gsdelFile%%.vcf.gz}.txt
zcat ${gscnvFile} | gs_del2tab -s ${ROOT}/sites.txt > ${gscnvFile%%.vcf.gz}.txt

# check source. Uses environment variables.
Rscript overlap.R

