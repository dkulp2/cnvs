Prediction:

Adjust configuration variables in ../conf/site.conf and ../conf/cnv.conf.

Run cnv_seg.sh
Output is in $workDir


merge_cnv.R - performs initial basic prediction using simple heuristics => csm
likelihoods.R - generates poisson probs along a chromome => pois and bkpt tables
staircase.R - improves each CNV by removing small gaps, dealing with staircase, and refining ends using MLE => smlcsm
collapse.R - flattens the (smlcsm) predictions into non-overlapping segments => smlx2csm
priors.R - takes collapse set (smlx2csm) and the transition likelihoods (bkpt) to generate priors => ??
posterior.R - takes prior set(s) and CNV predictions => cnv



