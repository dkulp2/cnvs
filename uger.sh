#!/bin/bash

export TAG=27Apr2017

qsub_std_args="-w e -j y -b y -V -cwd -r y"

# or -q [long|short] (short is up to four hours)
for Q in A B C D; do
echo   TAG=${TAG} qsub ${qsub_std_args} -l h_vmem=8g -q long -o log/$Q.log ./prep_and_run.sh $Q
done


# can also use /humgen/cnp04/svtoolkit/scripts/queue_uger_wrapper.sh or something similar.
