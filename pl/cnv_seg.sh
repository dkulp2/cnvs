#!/bin/bash -eu
#
# Run first pass genotyping on rolling windows.
# Create CNV segments
# Run second pass genotyping on CNV segments and GStrip DELs - which generates VCF of predicted genotypes per sample for IRS.
# Run IRS on 2nd pass genotype VCF

set -eux

ROOT=/cygdrive/d/mccarroll/

# External input data
# profileFile=${ROOT}/profiles/profile_seq_20_100.dat.gz
# arrayFile=${ROOT}/arrays/intensities.chr20.dat
profileFile=${ROOT}/gpc_wave2_batch1/profile_seq_20_100.dat.gz
arrayFile=${ROOT}/gpc_wave2_batch1/irs_matrix_OMNI25.dat
gsdelFile=${ROOT}/gpc_wave2_batch1/gs_dels.genotypes.vcf.gz
referenceFile=${ROOT}/references/Homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta

SAMPLE_COUNT=	      # number of samples or blank for all samples
SITE_COUNT=	      # number of rows from profile or blank for all data in profile
BIN_SIZE=100	      # number of bases in each bin in profile
NBINS=12	      # number of bins to combine into a single region
#NBINS2=8	      # number of bins for higher resolution "break point"
#MAXLEN=10000000  # effectively unlimited     # max length of resulting window. Is typically BIN_SIZE*NBIN = .e.g. 1000, but profile may have gaps or larger bins 
MAXLEN=2400     # max length of resulting window. Is typically BIN_SIZE*NBIN = .e.g. 1000, but profile may have gaps or larger bins 
CNV_CALL_THRESH=0.80  # fraction of samples required to use call from org.broadinstitute.sv.apps.ProfileGenotyper
CNQ_THRESH=13	      # phred-style threshold for org.broadinstitute.sv.apps.ProfileGenotyper (Q=20 => 1/10^(Q/10) => 99% confidence)
SPAN_THRESH=500      # max number of bases to span across if flanking segments are the same genotype

workDir=${ROOT}/cnv_seg.B${NBINS}.L${SPAN_THRESH}.Q${CNQ_THRESH}.4


# Java stuff
export WSP_DIR=${ROOT}/SVToolkit_Kulp
export SV_DIR=`cygpath -w ${WSP_DIR}/release/svtoolkit`  # passed to GATK
export SV_CLASSPATH=${WSP_DIR}/dist/SVToolkit-private.jar:${WSP_DIR}/public/dist/SVToolkit.jar:${WSP_DIR}/public/release/svtoolkit/lib/gatk/GenomeAnalysisTK.jar

# auxiliary scripts
SCRIPTS=${ROOT}/scripts
ROLLING_WINDOWS=${SCRIPTS}/choose_windows.awk
VCF2TAB=${SCRIPTS}/vcf2tab
MERGE_CNV=${SCRIPTS}/merge_cnv2.R
GSDEL2TAB=${SCRIPTS}/gs_del2tab

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

# # generates rolling windows of size ELENGTH*NBINS, e.g. 100*10 = 1000
# eval "gunzip -c ${profileFile} ${SITE_HEAD}" | awk -v OFS="\t" -v NBINS=${NBINS} -v MAXLEN=${MAXLEN} -f ${ROLLING_WINDOWS} > ${SITES}
# #eval "gunzip -c ${profileFile} ${SITE_HEAD}" | awk -v OFS="\t" -v NBINS=${NBINS2} -v MAXLEN=${MAXLEN} -f ${ROLLING_WINDOWS} > ${SITES_HIRES}

# # generate first pass genotypes of rolling windows
# time java -cp `cygpath -wp ${SV_CLASSPATH}` -Xmx4g \
#     org.broadinstitute.sv.apps.ProfileGenotyper \
#     -configFile `cygpath -w ${SV_DIR}/conf/genstrip_parameters.txt` \
#     -R `cygpath -w ${referenceFile}` \
#     -profile `cygpath -w ${profileFile}` \
#     -ploidyMapFile `cygpath -w ${referenceFile} | sed 's/.fasta$/.ploidymap.txt/'` \
#     -segmentFile `cygpath -w ${SITES}` \
#     -sample `cygpath -w ${SAMPLES}` \
#     -O `cygpath -w ${OUT1_VCF}` 

# # # run a second time with smaller windows
# # time java -cp `cygpath -wp ${SV_CLASSPATH}` -Xmx4g \
# #     org.broadinstitute.sv.apps.ProfileGenotyper \
# #     -configFile `cygpath -w ${SV_DIR}/conf/genstrip_parameters.txt` \
# #     -R `cygpath -w ${referenceFile}` \
# #     -profile `cygpath -w ${profileFile}` \
# #     -segmentFile `cygpath -w ${SITES_HIRES}` \
# #     -sample `cygpath -w ${SAMPLES}` \
# #     -O `cygpath -w ${OUT1_VCF_HIRES}` 

# # convert VCF to a simpler delimited file
# zcat ${OUT1_VCF} | ${VCF2TAB} > ${OUT1_VCF}.txt
# #zcat ${OUT1_VCF_HIRES} | ${VCF2TAB} > ${OUT1_VCF_HIRES}.txt

# merge windows with adjacent common CN, etc.
# write site file
time Rscript ${MERGE_CNV} ${OUT1_VCF}.txt ${CNV_CALL_THRESH} ${CNQ_THRESH} ${SPAN_THRESH} ${CNV_SEG_SITES_FILE} 1>&2
#time Rscript ${MERGE_CNV} ${OUT1_VCF_HIRES}.txt ${CNV_CALL_THRESH} ${CNQ_THRESH} ${SPAN_THRESH} ${CNV_SEG_SITES_FILE2} 1>&2

# One of the outputs is exploded genotype calls (one row per sample). Create a tabix index for use in viz
sort -k2n,3n ${CNV_SEG_SITES_FILE}.cnvgeno.txt > ${CNV_SEG_SITES_FILE}.cnvgeno.srt
bgzip -f ${CNV_SEG_SITES_FILE}.cnvgeno.srt
tabix -b 2 -e 3 -s 1 -S 1 ${CNV_SEG_SITES_FILE}.cnvgeno.srt.gz

#sort -k2n,3n ${CNV_SEG_SITES_FILE2}.cnvgeno.txt > ${CNV_SEG_SITES_FILE2}.cnvgeno.srt
#bgzip -f ${CNV_SEG_SITES_FILE2}.cnvgeno.srt
#tabix -b 2 -e 3 -s 1 -S 1 ${CNV_SEG_SITES_FILE2}.cnvgeno.srt.gz

# create a row per sample deletion for input into R
zcat ${gsdelFile} | ${GSDEL2TAB} -s ${SITES} > ${GS_DEL_FILE}

# create a site file for input into ProfileGentoyper
cut -f2-5 ${GS_DEL_FILE} | tail -n +1 | uniq > ${GS_DEL_SITES_FILE}

cat ${GS_DEL_SITES_FILE} ${CNV_SEG_SITES_FILE} | sort -n -k3,4 > ${IN2_SITES}

exit 0

# rerun org.broadinstitute.sv.apps.ProfileGenotyper on new CNV segment file
time java -cp `cygpath -wp ${SV_CLASSPATH}` -Xmx4g \
    org.broadinstitute.sv.apps.ProfileGenotyper \
    -configFile `cygpath -w ${SV_DIR}/conf/genstrip_parameters.txt` \
    -R `cygpath -w ${referenceFile}` \
    -profile `cygpath -w ${profileFile}` \
    -segmentFile `cygpath -w ${IN2_SITES}` \
    -sample `cygpath -w ${SAMPLES}` \
    -O `cygpath -w ${OUT2_VCF}` 


#
(for i in `cut -f2 ${SITES} | uniq`; do egrep "^[[:alnum:]]+[[:space:]]${i}[[:space:]]" ${arrayFile} | cut -f1-4; done) > ${workDir}/probes.txt

# run IRS on the output from these segments
time java -Xmx4g -cp `cygpath -wp ${SV_CLASSPATH}` \
     org.broadinstitute.sv.main.SVAnnotator \
     -A IntensityRankSum \
     -R `cygpath -w ${referenceFile}` \
     -vcf `cygpath -w ${OUT2_VCF}` \
     -O `cygpath -w ${IRS_VCF}` \
     -arrayIntensityFile `cygpath -w ${arrayFile}` \
     -irsUseGenotypes true \
     -writeReport true \
     -reportFile `cygpath -w ${IRS_REPORT}` 
