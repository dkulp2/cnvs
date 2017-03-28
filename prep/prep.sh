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

mkdir -p ${workDir}
eval "gunzip -c ${profileFile}" | awk -v OFS="\t" -v NBINS=${NBINS} -v MAXLEN=${MAXLEN} -f ${ROLLING_WINDOWS} > ${workDir}/sites.txt

# generate simplified tab files from VCF
if [ -f ${gsdelFile} ]; then
    zcat ${gsdelFile} | gs_del2tab -s ${workDir}/sites.txt > ${gsdelFile%%.vcf.gz}.txt
fi

if [ -f ${gscnvFile} ]; then
    zcat ${gscnvFile} | gs_del2tab -s ${workDir}/sites.txt > ${gscnvFile%%.vcf.gz}.txt
fi

# check source. Uses environment variables.
# Generate _flt, _xflt and _best sets from samplesegs
if [ -f ${gsdelFile} -a -f ${gscnvFile} ]; then
    Rscript overlap.R
fi

# Load database with profiles
perl ${THISDIR}/profile2db.pl -i ${profileFile} -s ${profileFile%%.dat.gz}.seg.txt.gz -p ${profileFile%%.dat.gz}.pro.txt.gz
sql ${THISDIR}/profile2db.sql

rm ${profileFile%%.dat.gz}.seg.txt.gz ${profileFile%%.dat.gz}.pro.txt.gz
