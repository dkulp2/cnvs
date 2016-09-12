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

DROP TABLE IF EXISTS PRED_CNVS;
CREATE TABLE IF NOT EXISTS PRED_CNVS (
       SAMPLE	     varchar(30),
       ID varchar(30),
       CHROM varchar(2),
       START_POS integer,
       END_POS integer,
       GENO varchar(3)
);

\COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.B12.L500.Q13/sites_cnv_segs.txt.sml.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.B12.L500.Q13.3/sites_cnv_segs.txt.ml.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.B12.L500.Q13.2/sites_cnv_segs.txt.ml.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.B12.L500.Q13/sites_cnv_segs.txt.ml.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.B12.L500.Q13/sites_cnv_segs.txt.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.12.500/sites_cnv_segs.txt.ml.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.12.500/old/sites_cnv_segs.txt.ml.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.12.500/sites_cnv_segs.txt.ext.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.12.500/sites_cnv_segs.txt.tbl' DELIMITER E'\t' CSV;
-- \COPY PRED_CNVS FROM '/cygdrive/d/mccarroll/cnv_seg.B12.L500.Q20/sites_cnv_segs.txt.ext.tbl' DELIMITER E'\t' CSV;

--------------------------------------------------------------------------------
-- For sensitivity, use smaller set of somewhat independently identified deletions

DELETE FROM KNOWN_CNVS;
-- \COPY KNOWN_CNVS FROM '/cygdrive/d/mccarroll/gpc_wave2_batch1/gs_dels.genotypes.txt' DELIMITER E'\t' CSV HEADER;
\COPY KNOWN_CNVS FROM '/cygdrive/d/mccarroll/gpc_wave2_batch1/gs_dels_flt.genotypes.txt' DELIMITER E'\t' CSV HEADER;

-- Manual curated "bad" DELS
DELETE FROM KNOWN_CNVS WHERE id IN ('DEL_P0561_239','DEL_P0561_240','DEL_P0561_242','DEL_P0557_48','DEL_P0566_78','DEL_P0557_41');

DROP TABLE IF EXISTS OVERLAP;
SELECT k.SAMPLE AS k_SAMPLE, p.SAMPLE AS p_SAMPLE, 
       k.ID AS k_ID, p.ID AS p_ID, k.GENO as k_GENO, p.GENO AS p_GENO,
       LEAST(k.END_POS, p.END_POS) - GREATEST(k.START_POS, p.START_POS) AS OVERLAP_LEN,
       k.END_POS - k.START_POS AS KNOWN_LEN,
       p.END_POS - p.START_POS AS PRED_LEN
  INTO OVERLAP
  FROM KNOWN_CNVS k FULL OUTER JOIN (SELECT * FROM PRED_CNVS WHERE START_POS < END_POS) AS p ON 
       (k.SAMPLE = p.SAMPLE AND
        k.CHROM = p.CHROM AND
	k.START_POS < p.END_POS AND
        k.END_POS > p.START_POS AND
        k.GENO = p.GENO);

-- LEAST() and GREATEST() ignore NULLS, so if there's no overlap, then OVERLAP_LEN
-- gets set to the length of the extent. Update the table and set to zero.
UPDATE OVERLAP
   SET overlap_len = 0
 WHERE known_len IS NULL OR pred_len IS NULL;

-- -- Sensitivity
-- DROP TABLE IF EXISTS base_sensitivity;
-- SELECT DISTINCT known_len AS length, tot_overlap_len::float/tot_known_len AS Sn
--   INTO base_sensitivity
--   FROM (
--        SELECT known_len, 
--               sum(known_len) OVER w AS tot_known_len,
--               sum(overlap_len) OVER w AS tot_overlap_len
--          FROM overlap
--         WHERE known_len IS NOT NULL
-- 	  AND geno = '1'
--        WINDOW w AS (ORDER BY known_len DESC)
--        ORDER BY known_len DESC
--        ) AS aaa1
-- ORDER BY known_len;

--------------------------------------------------------------------------------
-- More complete data set includes paired read DELs + CNV pipeline

DELETE FROM KNOWN_CNVS;
\COPY KNOWN_CNVS FROM '/cygdrive/d/mccarroll/gpc_wave2/gs_cnv_del.genotypes.txt' DELIMITER E'\t' CSV HEADER

-- SAME AS ABOVE EXCEPT restrict predictions by length (and use Sp knowns)
DROP TABLE IF EXISTS OVERLAP2;
SELECT k.SAMPLE AS k_SAMPLE, p.SAMPLE AS p_SAMPLE, 
       k.ID AS k_ID, p.ID AS p_ID, k.GENO as k_GENO, p.GENO as p_GENO,
       LEAST(k.END_POS, p.END_POS) - GREATEST(k.START_POS, p.START_POS) AS OVERLAP_LEN,
       k.END_POS - k.START_POS AS KNOWN_LEN,
       p.END_POS - p.START_POS AS PRED_LEN
  INTO OVERLAP2
  FROM KNOWN_CNVS k FULL OUTER JOIN (SELECT * FROM PRED_CNVS WHERE START_POS < END_POS) AS p ON 
       (k.SAMPLE = p.SAMPLE AND
        k.CHROM = p.CHROM AND
	k.START_POS < p.END_POS AND
        k.END_POS > p.START_POS AND
        k.GENO = p.GENO);

UPDATE OVERLAP2
   SET overlap_len = 0
 WHERE known_len IS NULL OR pred_len IS NULL;

-- -- Specificity THIS IS BROKEN, I THINK. Now computed in R.
-- DROP TABLE IF EXISTS base_specificity;
-- SELECT DISTINCT pred_len AS length, tot_overlap_len::float/tot_pred_len AS Sp
--   INTO base_specificity
--   FROM (
--        SELECT pred_len, 
--               sum(pred_len) OVER w AS tot_pred_len,
--               sum(overlap_len) OVER w AS tot_overlap_len
--          FROM overlap2
--         WHERE pred_len IS NOT NULL
-- 	  AND geno = '1'
--        WINDOW w AS (ORDER BY pred_len DESC)
--        ORDER BY pred_len DESC
--        ) AS aaa1
-- ORDER BY pred_len;

--------------------------------------------------------------------------------


\echo 'Base Sn where known length > 500'
SELECT k_geno as geno, SUM(overlap_len)::float / SUM(known_len) AS base_sn
  FROM overlap
 WHERE known_len IS NOT NULL
   AND known_len > 500
 GROUP BY k_geno;

\echo 'Base Sp where predicted length < 30k'
SELECT p_geno as geno, SUM(overlap_len)::float / SUM(pred_len) AS base_sp
  FROM overlap2
 WHERE pred_len IS NOT NULL
   AND pred_len < 30000
   AND pred_len > 0
GROUP BY p_geno;


--------------------------------------------------------------------------------
-- Site Sn/Sp

\echo 'Site Sn/Sp where overlap > 50%'

SELECT *, hit::float/total as site_sn
  FROM (SELECT count(distinct k_id) hit FROM overlap WHERE overlap_len::float / known_len > .5) AS A,
       (SELECT count(distinct k_id) total FROM overlap) AS B;

SELECT *, hit::float/total as site_sp
  FROM (SELECT count(distinct p_id) hit FROM overlap2 WHERE p_geno IN ('0','1') AND overlap_len::float / pred_len > .5) AS A,
       (SELECT count(distinct p_id) total FROM overlap2 WHERE p_geno IN ('0','1')) AS B;
