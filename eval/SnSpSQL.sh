#!/bin/bash

set -eux

THISDIR=`dirname $0`
export PATH=${THISDIR}/../util:$PATH
echo $PATH

# set ROOT and primary data sources
source ${THISDIR}/../conf/site.conf

# set params for this run
source ${THISDIR}/../conf/cnv.conf

psql -a <<EOF
CREATE SCHEMA eval;
SET search_path TO eval;

DROP TABLE IF EXISTS KNOWN_CNVS;
CREATE TABLE IF NOT EXISTS KNOWN_CNVS (
       SAMPLE	     varchar(30),
       ID varchar(10000),	-- this is so big because merging may create a long concatenated ID
       CHROM varchar(2),
       START_POS integer,
       END_POS integer,
       GENO varchar(3),
       QUAL float,
       READCOUNT integer);

DROP TABLE IF EXISTS ALL_KNOWN_CNVS;
CREATE TABLE ALL_KNOWN_CNVS (LIKE KNOWN_CNVS);
ALTER TABLE ALL_KNOWN_CNVS ADD COLUMN SOURCE VARCHAR(20);
ALTER TABLE ALL_KNOWN_CNVS ADD COLUMN KIND VARCHAR(20);

DROP TABLE IF EXISTS PRED_CNVS;
CREATE TABLE IF NOT EXISTS PRED_CNVS (
       SAMPLE	     varchar(30),
       ID varchar(100000),
       CHROM varchar(2),
       START_CIL integer,
       START_POS integer,
       START_CIR integer,
       END_CIL integer,
       END_POS integer,
       END_CIR integer,
       GENO varchar(3)
);
DROP TABLE IF EXISTS ALL_PRED_CNVS;
CREATE TABLE ALL_PRED_CNVS (LIKE PRED_CNVS);
ALTER TABLE ALL_PRED_CNVS ADD COLUMN SOURCE VARCHAR(20);
ALTER TABLE ALL_PRED_CNVS ADD COLUMN KIND VARCHAR(20);


-- ".smlCI" are basic predictions, extended by "staircase" and MLE, with confidence intervals
DELETE FROM PRED_CNVS;
\COPY PRED_CNVS FROM '${workDir}/sites_cnv_segs.txt.smlCI.tbl' DELIMITER E'\t' CSV NULL 'NA';
INSERT INTO ALL_PRED_CNVS SELECT *, 'SML','SAMPLESEG' FROM PRED_CNVS;
INSERT INTO ALL_PRED_CNVS SELECT *, 'SML','SAMPLESEGSP' FROM PRED_CNVS;

-- ".bayesCI" are basic predictions, extended by "staircase" and MLE, and adjusted using a
-- bayesian approach to incorporate priors.
DELETE FROM PRED_CNVS;
\COPY PRED_CNVS FROM '${workDir}/sites_cnv_segs.txt.bayesCI.tbl' DELIMITER E'\t' CSV NULL 'NA';
INSERT INTO ALL_PRED_CNVS SELECT *, 'BAYES','SAMPLESEG' FROM PRED_CNVS;
INSERT INTO ALL_PRED_CNVS SELECT *, 'BAYES','SAMPLESEGSP' FROM PRED_CNVS;

-- ".flt" is a flattened ".sml" for site-level analysis
DELETE FROM PRED_CNVS;
\COPY PRED_CNVS FROM '${workDir}/sites_cnv_segs.txt.flt.tbl' DELIMITER E'\t' CSV;
INSERT INTO ALL_PRED_CNVS SELECT *, 'SML','SITE' FROM PRED_CNVS;
INSERT INTO ALL_PRED_CNVS SELECT *, 'SML','SITESP' FROM PRED_CNVS;

-- ".bayesflt" is a flattened "bayes" for site-level analysis
DELETE FROM PRED_CNVS;
\COPY PRED_CNVS FROM '${workDir}/sites_cnv_segs.txt.bayesflt.tbl' DELIMITER E'\t' CSV;
INSERT INTO ALL_PRED_CNVS SELECT *, 'BAYES','SITE' FROM PRED_CNVS;
INSERT INTO ALL_PRED_CNVS SELECT *, 'BAYES','SITESP' FROM PRED_CNVS;

-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.B12.L500.Q13.3/sites_cnv_segs.txt.ml.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.B12.L500.Q13.2/sites_cnv_segs.txt.ml.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.B12.L500.Q13/sites_cnv_segs.txt.ml.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.B12.L500.Q13/sites_cnv_segs.txt.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.12.500/sites_cnv_segs.txt.ml.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.12.500/old/sites_cnv_segs.txt.ml.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.12.500/sites_cnv_segs.txt.ext.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.12.500/sites_cnv_segs.txt.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.B12.L500.Q20/sites_cnv_segs.txt.ext.tbl' DELIMITER E'\t' CSV;

-- GSTRIP DELETION pipeline
-- DELETE FROM KNOWN_CNVS;
-- \COPY KNOWN_CNVS FROM '/cygdrive/d/mccarroll/gpc_wave2_batch1/gs_dels.genotypes.txt' DELIMITER E'\t' CSV HEADER;
-- INSERT INTO ALL_KNOWN_CNVS SELECT *, 'GS_DELS', 'SAMPLESEG' FROM KNOWN_CNVS;

-- FLATTENED GSTRIP DELETION pipeline
DELETE FROM KNOWN_CNVS;
\COPY KNOWN_CNVS FROM '${workDir}/../../gpc_wave2_batch1/gs_dels_flt.genotypes.txt' DELIMITER E'\t' CSV HEADER;
INSERT INTO ALL_KNOWN_CNVS SELECT *, 'GS_DELS', 'SAMPLESEG' FROM KNOWN_CNVS;

-- EXTRA-FLATTENED GSTRIP DELETION pipeline (OVER SAMPLES, TOO)
DELETE FROM KNOWN_CNVS;
\COPY KNOWN_CNVS FROM '${workDir}/../../gpc_wave2_batch1/gs_dels_xflt.genotypes.txt' DELIMITER E'\t' CSV HEADER;
INSERT INTO ALL_KNOWN_CNVS SELECT *, 'GS_DELS', 'SITE' FROM KNOWN_CNVS;

-- CNV+GSDEL PIPELINE, FLATTENED
DELETE FROM KNOWN_CNVS;
\COPY KNOWN_CNVS FROM '${workDir}/../../gpc_wave2/gs_cnv_del_flt.genotypes.txt' DELIMITER E'\t' CSV HEADER
INSERT INTO ALL_KNOWN_CNVS SELECT *, 'GS_CNV_DEL', 'SAMPLESEGSP' FROM KNOWN_CNVS;

-- CNV+GSDEL PIPELINE, EXTRA-FLATTENED (OVER SAMPLES, TOO) 
DELETE FROM KNOWN_CNVS;
\COPY KNOWN_CNVS FROM '${workDir}/../../gpc_wave2/gs_cnv_del_xflt.genotypes.txt' DELIMITER E'\t' CSV HEADER
INSERT INTO ALL_KNOWN_CNVS SELECT *, 'GS_CNV_DEL', 'SITESP' FROM KNOWN_CNVS;

-- Remove manually curated "bad" DELS
DELETE FROM ALL_KNOWN_CNVS WHERE kind in ('SAMPLESEG','SITE') AND
  (id LIKE '%DEL_P0561_239%' OR id LIKE '%DEL_P0561_240%' OR id LIKE '%DEL_P0561_242%' OR
   id LIKE '%DEL_P0557_48%' or id LIKE '%DEL_P0566_78%' OR id LIKE '%DEL_P0557_41%');

DROP TABLE IF EXISTS OVERLAP;
SELECT k.SAMPLE AS k_SAMPLE, k.SOURCE as k_SOURCE, p.SAMPLE AS p_SAMPLE, p.SOURCE as p_SOURCE,
       k.ID AS k_ID, p.ID AS p_ID, k.GENO as k_GENO, p.GENO AS p_GENO, k.KIND as k_KIND, p.KIND AS p_KIND,
       LEAST(k.END_POS, p.END_POS) - GREATEST(k.START_POS, p.START_POS) AS OVERLAP_LEN,
       k.END_POS - k.START_POS AS KNOWN_LEN,
       p.END_POS - p.START_POS AS PRED_LEN,
       p.START_POS - k.START_POS AS START_OFFSET,
       p.END_POS - k.END_POS AS END_OFFSET,
       k.START_POS BETWEEN p.START_CIL AND p.START_CIR AS START_IN_CIL,
       k.END_POS BETWEEN p.END_CIL AND p.END_CIR AS END_IN_CIL
  INTO OVERLAP
  FROM ALL_KNOWN_CNVS k FULL OUTER JOIN (SELECT * FROM ALL_PRED_CNVS WHERE START_POS < END_POS AND source='BAYES') AS p ON 
       (k.SAMPLE = p.SAMPLE AND
        k.CHROM = p.CHROM AND
	k.START_POS < p.END_POS AND
        k.END_POS > p.START_POS AND
        (k.GENO = p.GENO OR k.GENO = '99') AND  -- Match on CN or any predicted CN if actual is discordant
        k.KIND = p.KIND
);

-- LEAST() and GREATEST() ignore NULLS, so if there's no overlap, then OVERLAP_LEN
-- gets set to the length of the extent. Update the table and set to zero.
UPDATE OVERLAP
   SET overlap_len = 0
 WHERE known_len IS NULL OR pred_len IS NULL;

EOF


# -- -- Sensitivity
# -- DROP TABLE IF EXISTS base_sensitivity;
# -- SELECT DISTINCT known_len AS length, tot_overlap_len::float/tot_known_len AS Sn
# --   INTO base_sensitivity
# --   FROM (
# --        SELECT known_len, 
# --               sum(known_len) OVER w AS tot_known_len,
# --               sum(overlap_len) OVER w AS tot_overlap_len
# --          FROM overlap
# --         WHERE known_len IS NOT NULL
# -- 	  AND geno = '1'
# --        WINDOW w AS (ORDER BY known_len DESC)
# --        ORDER BY known_len DESC
# --        ) AS aaa1
# -- ORDER BY known_len;

# --------------------------------------------------------------------------------
# -- More complete data set includes paired read DELs + CNV pipeline

# DELETE FROM KNOWN_CNVS;
# \COPY KNOWN_CNVS FROM '/cygdrive/d/mccarroll/gpc_wave2/gs_cnv_del.genotypes.txt' DELIMITER E'\t' CSV HEADER

# -- SAME AS ABOVE EXCEPT restrict predictions by length (and use Sp knowns)
# DROP TABLE IF EXISTS OVERLAP2;
# SELECT k.SAMPLE AS k_SAMPLE, p.SAMPLE AS p_SAMPLE, 
#        k.ID AS k_ID, p.ID AS p_ID, k.GENO as k_GENO, p.GENO as p_GENO,
#        LEAST(k.END_POS, p.END_POS) - GREATEST(k.START_POS, p.START_POS) AS OVERLAP_LEN,
#        k.END_POS - k.START_POS AS KNOWN_LEN,
#        p.END_POS - p.START_POS AS PRED_LEN
#   INTO OVERLAP2
#   FROM KNOWN_CNVS k FULL OUTER JOIN (SELECT * FROM PRED_CNVS WHERE START_POS < END_POS) AS p ON 
#        (k.SAMPLE = p.SAMPLE AND
#         k.CHROM = p.CHROM AND
# 	k.START_POS < p.END_POS AND
#         k.END_POS > p.START_POS AND
#         k.GENO = p.GENO);

# UPDATE OVERLAP2
#    SET overlap_len = 0
#  WHERE known_len IS NULL OR pred_len IS NULL;

# -- -- Specificity THIS IS BROKEN, I THINK. Now computed in R.
# -- DROP TABLE IF EXISTS base_specificity;
# -- SELECT DISTINCT pred_len AS length, tot_overlap_len::float/tot_pred_len AS Sp
# --   INTO base_specificity
# --   FROM (
# --        SELECT pred_len, 
# --               sum(pred_len) OVER w AS tot_pred_len,
# --               sum(overlap_len) OVER w AS tot_overlap_len
# --          FROM overlap2
# --         WHERE pred_len IS NOT NULL
# -- 	  AND geno = '1'
# --        WINDOW w AS (ORDER BY pred_len DESC)
# --        ORDER BY pred_len DESC
# --        ) AS aaa1
# -- ORDER BY pred_len;

# --------------------------------------------------------------------------------


# \echo 'Base Sn where known length > 500'
# SELECT k_geno as geno, SUM(overlap_len)::float / SUM(known_len) AS base_sn
#   FROM overlap
#  WHERE known_len IS NOT NULL
#    AND known_len > 500
#  GROUP BY k_geno;

# \echo 'Base Sp where predicted length < 30k'
# SELECT p_geno as geno, SUM(overlap_len)::float / SUM(pred_len) AS base_sp
#   FROM overlap2
#  WHERE pred_len IS NOT NULL
#    AND pred_len < 30000
#    AND pred_len > 0
# GROUP BY p_geno;


# --------------------------------------------------------------------------------
# -- Site Sn/Sp

# \echo 'Site Sn/Sp where overlap > 50%'

# SELECT *, hit::float/total as site_sn
#   FROM (SELECT count(distinct k_id) hit FROM overlap WHERE overlap_len::float / known_len > .5) AS A,
#        (SELECT count(distinct k_id) total FROM overlap) AS B;

# SELECT *, hit::float/total as site_sp
#   FROM (SELECT count(distinct p_id) hit FROM overlap2 WHERE p_geno IN ('0','1') AND overlap_len::float / pred_len > .5) AS A,
#        (SELECT count(distinct p_id) total FROM overlap2 WHERE p_geno IN ('0','1')) AS B;
