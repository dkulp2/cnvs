#!/bin/bash

export TAG=19May2017

qsub_std_args="-w e -j y -b y -V -cwd -r y"

# or -q [long|short] (short is up to four hours)
for Q in A B C D; do
echo   TAG=${TAG} qsub ${qsub_std_args} -l h_vmem=8g -q long -o log/$Q.log ./prep_and_run.sh $Q
done

. conf/site.conf
. conf/cnv.conf
echo When jobs are completed run:
echo Rscript pl/mk_quartet_segments.R ${ROOT}/out/${TAG}/data_sfari_batch1 _${TAG}/`basename ${workDir}`/sites_cnv_segs.txt /humgen/cnp04/bobh/projects/sfari/study_data/sample_pedigrees.ped /humgen/gsa-hpprojects/dev/cwhelan/SFARI/ibd_results/merged_ibd_regions.bed bayescsm


# can also use /humgen/cnp04/svtoolkit/scripts/queue_uger_wrapper.sh or something similar.
