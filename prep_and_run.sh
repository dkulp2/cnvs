#!/bin/sh
# run a test on a quartet set and set output to use most recent git tag
#
# Usage: prep_and_run.sh A
#
# Runs on the A quartet. Uncomment prep or cpPrep depending on circumstance

echo ====================
date

echo Using tag: ${TAG}
echo Quartet: $1
LQ=`echo $1 | tr '[:upper:]' '[:lower:]'`

export QUARTET_DIR=/humgen/cnp04/bobh/projects/segmentation/sfari/data_sfari_batch1${1}
echo QUARTET_DIR=${QUARTET_DIR}

# time sh prep/prep.sh
# time sh prep/cpPrep.sh data_sfari_batch1${1}_11apr2017b ~/data/out/11Apr2017/data_sfari_batch1${1}_11Apr2017b/data_sfari_batch1${LQ}_11apr2017b
time sh pl/cnv_seg.sh

