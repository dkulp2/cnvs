????

DROP TABLE IF EXISTS breakpoints_1kg;
CREATE TABLE breakpoints_1kg (
       CHROM varchar(2),
       START_POS integer,
       END_POS integer,
       ID varchar(10000),
       CNV_TYPE	     varchar(30),
       INSSEQ varchar(10000),
       HOMLEN integer,
       HOMSEQ varchar(10000),
       MUTMECH varchar(10),
       ANCESTRAL varchar(20)
);

\COPY breakpoints_1kg FROM '/cygdrive/d/mccarroll/KGen/1KG_phase3_all_bkpts.v5.chr20.L1500.txt' DELIMITER E'\t' CSV;

-- Find predictions that overlap these deletions
DROP TABLE IF EXISTS overlap3;
SELECT p.*, b.start_pos as b_start, b.end_pos as b_end
  INTO overlap3
  FROM pred_cnvs p, breakpoints_1kg b
 WHERE p.START_POS < b.END_POS AND p.END_POS > b.START_POS AND p.CHROM=b.CHROM AND p.GENO IN ('1','0')
   AND (LEAST(b.END_POS, p.END_POS) - GREATEST(b.START_POS, p.START_POS) / (b.END_POS - b.START_POS)::float) > 0.75

