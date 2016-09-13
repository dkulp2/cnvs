#!/bin/bash
#
# Massage known data sets for SnSp analysis in SQL

# halt on error, report uninitialized and show commands
set -eux

ROOT=/cygdrive/d/mccarroll/
SCRIPTS=${ROOT}/scripts
WORK=${ROOT}/cnv_seg.12.500     # change this only if chromosome region changes. Not necessary to rerun for each new test run.

SITES=${WORK}/sites.txt		# a file of regions to analysis. Use to limit to chr sub-regions

GSDEL2TAB=${SCRIPTS}/gs_del2tab

# Convert VCF into tab-delimited file - one line per sample del
gsdelFile=${ROOT}/gpc_wave2_batch1/gs_dels.genotypes.vcf.gz
gsdelTab=${gsdelFile%.vcf.gz}.txt
zcat ${gsdelFile} | ${GSDEL2TAB} -s $SITES > ${gsdelTab}

# Convert VCF into tab-delimited file - limit to CNVs at least 1200nt and CNQ>=13
gscnvFile=${ROOT}/gpc_wave2/gs_cnv.genotypes.vcf.gz
gscnvTab=${gscnvFile%.vcf.gz}.txt
zcat ${gscnvFile} | ${GSDEL2TAB} -s $SITES -l 1200 -q 13 > ${gscnvTab}

# For Sn analysis, use just paired reads from DEL pipeline,
# but for Sp analysis, merge CNV and DEL pipeline.
# The merging must be done only with regions that share the same sample
# and copy number.
# Use gs_del2tab a third time to create separate files for each sample/CN combo,
# convert each to bed, run bedtool's merge on each and convert back to gs_del2tab format.
gscnvdelTab=`dirname ${gscnvFile}`/gs_cnv_del.genotypes.txt
TMPMRGDIR=/tmp/mrg
mkdir -p ${TMPMRGDIR}; rm -f ${TMPMRGDIR}/*.txt

# split all into separate files
zcat ${gsdelFile} | ${GSDEL2TAB} -s $SITES -d ${TMPMRGDIR}
zcat ${gscnvFile} | ${GSDEL2TAB} -s $SITES -d ${TMPMRGDIR} -l 1200 -q 13

# convert each to BED and then back to the custom tab format
# there is no quality or readcount for merged results, so just write placeholders
# the sample and CN/GT is in the filename
for f in ${TMPMRGDIR}/*.txt; do
    echo $f
    if [ -s $f ]; then
	# custom GSDEL2TAB format is SAMPLE LABEL CHR START END GT GQ READCOUNT (here w/o header)
	# BED file is of form CHR START END LABEL
	perl -nae 'print "$F[2]\t$F[3]\t$F[4]\t$F[1]\n"' $f | sort -k1,1 -k2,3n > ${f%.txt}.bed
	bedtools merge -i ${f%.txt}.bed -c 4 -o collapse > ${f%.txt}.mrg.bed

	# extract sample and CN/GT from filename
	f2=`basename $f`
	sample=${f2%%.*}
	suffix=${f2#${sample}.}
	cngt=${suffix%%.*}
    
	# convert back to custom tab format, adding placeholders
	perl -nae "print \"$sample\\t\$F[3]\\t\$F[0]\\t\$F[1]\\t\$F[2]\\t$cngt\\t0\\t0\\n\"" ${f%.txt}.mrg.bed > ${f%.txt}.mrg.txt
    fi
done
(head -1 $gscnvTab; cat ${TMPMRGDIR}/*.mrg.txt) > ${gscnvdelTab}

