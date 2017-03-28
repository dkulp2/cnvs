Prediction:

Adjust configuration variables in ../conf/site.conf and ../conf/cnv.conf.

Run cnv_seg.sh
Output is in $workDir


1. merge_cnv.R - performs initial basic prediction using simple heuristics => csm
2. likelihoods.R - generates poisson probs along a chromome => pois and bkpt tables
3. staircase.R - improves each CNV by removing small gaps, dealing with staircase, and refining ends using MLE => smlcsm
4. collapse.R - flattens the (smlcsm) predictions into non-overlapping segments => smlx2csm
5. priors.R - takes collapse set (smlx2csm) and the transition likelihoods (bkpt) to generate priors => "prior" and "prior_region" tables
6. posterior.R - takes prior set(s) and CNV predictions => cnvs (db table - not labeled (sort of temporary, fix me)) and file (bayescsm and bayesCI.tbl) 


train/test is distinguished in database by label and CNV sets stored in Rdata files are segregated by directory

1-5 is run on train
1-6 is run on test with train's label passed to posterior.R


