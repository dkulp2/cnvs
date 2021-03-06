#!/bin/sh
# run a test on a quartet set and set output to use most recent git tag

export TAG=11Apr2017
echo Using tag: ${TAG}

for Q in A B C D; do
  export QUARTET_DIR=/humgen/cnp04/bobh/projects/segmentation/sfari/data_sfari_batch1${Q}
  time sh prep/prep.sh
  time sh pl/cnv_seg.sh
done
