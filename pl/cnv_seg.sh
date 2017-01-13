#!/bin/bash -eu
#
# First pass: profileGenotype on fixed width rolling windows.
# CNV Calls: Run merge_cnv.R and staircase.R to generate predictions
# Second pass: genotype on CNV segments and GStrip DELs - which generates VCF of predicted genotypes per sample, which can be used for IRS.
# IRS evaluation: Determine the quality of the deletion based on probe data for CNV segments per sample

set -eux

THISDIR=`dirname $0`
export PATH=${THISDIR}/../util:$PATH
echo $PATH

# set ROOT and primary data sources
source ${THISDIR}/../conf/site.conf

# set params for this run
source ${THISDIR}/../conf/cnv.conf

# TBD: generalize this for test/train. Currently this is incestuous.
INT_LABEL=${LABEL}  # label for internal (test) prior
EXT_LABEL=${LABEL}  # label for external prior

# auxiliary scripts
SCRIPTS=${THISDIR}
UTILS=${THISDIR}/../util
ROLLING_WINDOWS=${SCRIPTS}/choose_windows.awk
VCF2TAB=${UTILS}/vcf2tab
GSDEL2TAB=${UTILS}/gs_del2tab
MERGE_CNV=${SCRIPTS}/merge_cnv.R
LIKELIHOODS=${SCRIPTS}/likelihoods.R
STAIRCASE=${SCRIPTS}/staircase.R
COLLAPSE=${SCRIPTS}/collapse.R
PRIORS=${SCRIPTS}/priors.R
POSTERIOR=${SCRIPTS}/posterior.R


# input and output files
SAMPLES=${workDir}/samples${SAMPLE_COUNT}.list
SITES=${workDir}/sites${SITE_COUNT}.txt
SITES_HIRES=${workDir}/sites${SITE_COUNT}_hires.txt

CNV_SEG_SITES_FILE=${workDir}/sites_cnv_segs.txt
CNV_SEG_SITES_FILE2=${workDir}/hires_sites_cnv_segs.txt

GS_DEL_FILE=${workDir}/gs_dels.txt # one row per sample site deletion
GS_DEL_SITES_FILE=${workDir}/sites_gs_dels.txt

# 1st and 2nd pass calls to profile genotyper
OUT1_VCF=${workDir}/windows.vcf.gz
OUT1_VCF_HIRES=${workDir}/windows.hires.vcf.gz
IN2_SITES=${workDir}/sites_in2.txt
OUT2_VCF=${workDir}/cnv_segs.vcf

# output from IRS 
IRS_REPORT=${OUT2_VCF%.vcf*}.irs
IRS_VCF=${IRS_REPORT}.vcf

################################################################################

mkdir -p $workDir

if [ ! -z ${SAMPLE_COUNT} ]; then
    SAMPLE_HEAD="| head -${SAMPLE_COUNT}"
else
    SAMPLE_HEAD=""
fi
eval "gunzip -c ${profileFile} | head -1 | transpose -t | tail -n +7 ${SAMPLE_HEAD}" > ${SAMPLES}

if [ ! -z ${SITE_COUNT} ]; then
    SITE_HEAD="| head -${SITE_COUNT}"
else
    SITE_HEAD=""
fi

# generates rolling windows of size ELENGTH*NBINS, e.g. 100*10 = 1000
eval "gunzip -c ${profileFile} ${SITE_HEAD}" | awk -v OFS="\t" -v NBINS=${NBINS} -v MAXLEN=${MAXLEN} -f ${ROLLING_WINDOWS} > ${SITES}
#eval "gunzip -c ${profileFile} ${SITE_HEAD}" | awk -v OFS="\t" -v NBINS=${NBINS2} -v MAXLEN=${MAXLEN} -f ${ROLLING_WINDOWS} > ${SITES_HIRES}

# generate first pass genotypes of rolling windows
time java -cp `cygpath.shim -wp ${SV_CLASSPATH}` -Xmx4g \
    org.broadinstitute.sv.apps.ProfileGenotyper \
    -configFile `cygpath.shim -w ${SV_DIR}/conf/genstrip_parameters.txt` \
    -R `cygpath.shim -w ${referenceFile}` \
    -profile `cygpath.shim -w ${profileFile}` \
    -ploidyMapFile `cygpath.shim -w ${referenceFile} | sed 's/.fasta$/.ploidymap.txt/'` \
    -segmentFile `cygpath.shim -w ${SITES}` \
    -sample `cygpath.shim -w ${SAMPLES}` \
    -O `cygpath.shim -w ${OUT1_VCF}` 

# # # run a second time with smaller windows
# # time java -cp `cygpath.shim -wp ${SV_CLASSPATH}` -Xmx4g \
# #     org.broadinstitute.sv.apps.ProfileGenotyper \
# #     -configFile `cygpath.shim -w ${SV_DIR}/conf/genstrip_parameters.txt` \
# #     -R `cygpath.shim -w ${referenceFile}` \
# #     -profile `cygpath.shim -w ${profileFile}` \
# #     -segmentFile `cygpath.shim -w ${SITES_HIRES}` \
# #     -sample `cygpath.shim -w ${SAMPLES}` \
# #     -O `cygpath.shim -w ${OUT1_VCF_HIRES}` 

# # convert VCF to a simpler delimited file
zcat ${OUT1_VCF} | ${VCF2TAB} > ${OUT1_VCF}.txt
# #zcat ${OUT1_VCF_HIRES} | ${VCF2TAB} > ${OUT1_VCF_HIRES}.txt

# basic CNV prediction
time Rscript ${MERGE_CNV} ${OUT1_VCF}.txt ${CNV_CALL_THRESH} ${CNQ_THRESH} ${SPAN_THRESH} ${CNV_SEG_SITES_FILE} 1>&2

# generate poisson probs
time Rscript ${LIKELIHOODS} ${OUT1_VCF}.txt ${DBCONN} ${MLE_WINSIZE} ${LABEL}

# filter and adjust
time Rscript ${STAIRCASE} ${CNV_SEG_SITES_FILE} csm ${DBCONN} ${SPAN_THRESH} $((BIN_SIZE * NBINS)) ${MLE_WINSIZE} ${LABEL} 1>&2

# generate a collapse "site" set => smlx2csm
time Rscript ${COLLAPSE} ${CNV_SEG_SITES_FILE} smlcsm smlxcsm flt smlx2csm

# generate prior distros
time Rscript ${PRIORS} ${CNV_SEG_SITES_FILE} smlx2csm ${DBCONN} ${LABEL} ${MLE_WINSIZE}

time Rscript ${POSTERIOR} ${CNV_SEG_SITES_FILE} smlcsm ${DBCONN} ${INT_LABEL} ${EXT_LABEL} ${PRIOR_BLEND} ${MLE_WINSIZE}

# One of the outputs is exploded genotype calls (one row per sample). Create a tabix index for use in viz
sort -k2n,3n ${CNV_SEG_SITES_FILE}.cnvgeno.txt > ${CNV_SEG_SITES_FILE}.cnvgeno.srt
bgzip -f ${CNV_SEG_SITES_FILE}.cnvgeno.srt
tabix -b 2 -e 3 -s 1 -S 1 ${CNV_SEG_SITES_FILE}.cnvgeno.srt.gz

#sort -k2n,3n ${CNV_SEG_SITES_FILE2}.cnvgeno.txt > ${CNV_SEG_SITES_FILE2}.cnvgeno.srt
#bgzip -f ${CNV_SEG_SITES_FILE2}.cnvgeno.srt
#tabix -b 2 -e 3 -s 1 -S 1 ${CNV_SEG_SITES_FILE2}.cnvgeno.srt.gz

# create a row per sample deletion for input into R
if [ -f $gsdelFile ]; then 
  zcat ${gsdelFile} | ${GSDEL2TAB} -s ${SITES} > ${GS_DEL_FILE}

  # create a site file for input into ProfileGentoyper
  cut -f2-5 ${GS_DEL_FILE} | tail -n +1 | uniq > ${GS_DEL_SITES_FILE}

  cat ${GS_DEL_SITES_FILE} ${CNV_SEG_SITES_FILE} | sort -n -k3,4 > ${IN2_SITES}
fi


exit 0

# rerun org.broadinstitute.sv.apps.ProfileGenotyper on new CNV segment file
time java -cp `cygpath.shim -wp ${SV_CLASSPATH}` -Xmx4g \
    org.broadinstitute.sv.apps.ProfileGenotyper \
    -configFile `cygpath.shim -w ${SV_DIR}/conf/genstrip_parameters.txt` \
    -R `cygpath.shim -w ${referenceFile}` \
    -profile `cygpath.shim -w ${profileFile}` \
    -segmentFile `cygpath.shim -w ${IN2_SITES}` \
    -sample `cygpath.shim -w ${SAMPLES}` \
    -O `cygpath.shim -w ${OUT2_VCF}` 


#
(for i in `cut -f2 ${SITES} | uniq`; do egrep "^[[:alnum:]]+[[:space:]]${i}[[:space:]]" ${arrayFile} | cut -f1-4; done) > ${workDir}/probes.txt

# run IRS on the output from these segments
time java -Xmx4g -cp `cygpath.shim -wp ${SV_CLASSPATH}` \
     org.broadinstitute.sv.main.SVAnnotator \
     -A IntensityRankSum \
     -R `cygpath.shim -w ${referenceFile}` \
     -vcf `cygpath.shim -w ${OUT2_VCF}` \
     -O `cygpath.shim -w ${IRS_VCF}` \
     -arrayIntensityFile `cygpath.shim -w ${arrayFile}` \
     -irsUseGenotypes true \
     -writeReport true \
     -reportFile `cygpath.shim -w ${IRS_REPORT}` 
