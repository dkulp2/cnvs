CREATE TABLE geno (chr char(2), start integer, finish integer, sample varchar(20), cn integer, cnq numeric);
\copy geno from 'sites_cnv_segs.txt.cnvgeno.txt' DELIMITER E'\t' NULL 'NA' HEADER CSV
ALTER TABLE geno ADD COLUMN bin INTEGER;
CREATE INDEX foo ON geno(chr,start);
UPDATE geno SET bin = ps.bin FROM profile_segment ps WHERE ps.chrom=geno.chr AND ps.start_pos=geno.start; -- (this is non-standard SQL)
DROP INDEX foo;
CREATE UNIQUE INDEX ON geno(chr,bin,sample);
