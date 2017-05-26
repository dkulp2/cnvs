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
OUTLIER_BIN_FILTER=${SCRIPTS}/outlier_bin_filter.pl

mkdir -p ${workDir}

# Assume that if the profileFile is missing, then it is of the form dat.out.gz and use the original file of dat.gz to transform.
if [ ! -f ${profileFile} ]; then
    perl ${OUTLIER_BIN_FILTER} ${OUTLIER_MULTIPLE} ${profileFile%%.out.gz}.gz | bgzip -c > ${profileFile}
    tabix -S 1 -s 2 -b 3 -e 4 -f ${profileFile}
fi

# generates rolling windows of size ELENGTH*NBINS, e.g. 100*10 = 1000
if [ ! -f ${SITES} ]; then
    zcat ${profileFile} | awk -v OFS="\t" -v NBINS=${NBINS} -v MAXLEN=${MAXLEN} -f ${ROLLING_WINDOWS} > ${SITES}
fi


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
if ! tableExists profile_segment || ! tableExists profile_counts; then
    perl ${THISDIR}/profile2db.pl -i ${profileFile} -s ${profileFile%%.dat.gz}.seg.txt.gz -p ${profileFile%%.dat.gz}.pro.txt.gz
    sql ${THISDIR}/profile2db.sql
fi
rm -f ${profileFile%%.dat.gz}.seg.txt.gz ${profileFile%%.dat.gz}.pro.txt.gz

