#!/bin/sh
# run a test on a quartet set and set output to use most recent git tag
#
# Usage: prep_and_run.sh A
#
# Runs on the A quartet. Uncomment prep or cpPrep depending on circumstance

set -eux

echo ====================
date

echo Using tag: ${TAG}
echo Quartet: $1
LQ=`echo $1 | tr '[:upper:]' '[:lower:]'`

export QUARTET_DIR=/humgen/cnp04/bobh/projects/segmentation/sfari/data_sfari_batch1${1}
echo QUARTET_DIR=${QUARTET_DIR}

OLDTAG=26May2017_24bin
time sh prep/cpPrep.sh data_sfari_batch1${LQ}_${OLDTAG} ~/data/out/${OLDTAG}/data_sfari_batch1${1}_${OLDTAG}/B24.L5.Q13.W10.PB0.7.ML1e7
time sh prep/prep.sh
time sh pl/cnv_seg.sh




