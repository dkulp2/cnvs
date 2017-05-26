#!/bin/sh
#
# Create a schema that are views to another schema.
# Create links to files in another run.
#
# Usage: cpPrep.sh from-schema-name from-work-dir
#
# Destination schema and workDir are known from environment variables

MIGRATION_SCHEMA=$1
MIGRATION_DIR=$2

THISDIR=`dirname $0`
export PATH=${THISDIR}/../util:$PATH
source ${THISDIR}/../conf/site.conf
source ${THISDIR}/../conf/cnv.conf

mkdir -p ${workDir}
if [ -f ${gsdelFile%%.vcf.gz}.txt ]; then
    ln -s ${gsdelFile%%.vcf.gz}.txt ${MIGRATION_DIR}
fi
if [ -f ${gscnvFile%%.vcf.gz}.txt ]; then
    ln -s ${gscnvFile%%.vcf.gz}.txt ${MIGRATION_DIR}
fi

if [ -f ${gsdelFile} -a -f ${gscnvFile} ]; then
    Rscript overlap.R
fi

psql -a <<EOF
CREATE SCHEMA IF NOT EXISTS ${MIGRATION_SCHEMA}
CREATE VIEW ${MIGRATION_SCHEMA}.profile_segment AS SELECT * FROM profile_segment;
CREATE VIEW ${MIGRATION_SCHEMA}.profile_counts AS SELECT * FROM profile_counts;
EOF
