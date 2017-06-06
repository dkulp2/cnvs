BEGIN TRANSACTION;

DROP TABLE IF EXISTS geno CASCADE;

CREATE TABLE geno (
       bin integer,
       chr text,
       start_pos integer,
       end_pos integer,
       sample text,
       cn integer,
       cnq double precision
);

\copy geno FROM PROGRAM 'zcat ${CNV_SEG_SITES_FILE}.cnvgeno.srt.gz' NULL 'NA' DELIMITER E'\t' CSV HEADER

CREATE UNIQUE INDEX on geno(chr, bin, sample);

COMMIT;

--VACUUM ANALYZE;

