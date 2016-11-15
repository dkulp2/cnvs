DROP TABLE profile_segment CASCADE;
CREATE TABLE profile_segment (
       BIN integer PRIMARY KEY,
       CHROM varchar(2),
       START_POS integer,
       END_POS integer,
       ELENGTH integer,
       GC integer,
       GCTOTAL integer
);

DROP TABLE profile_counts;
CREATE TABLE profile_counts (
       BIN integer REFERENCES profile_segment(bin),
       SAMPLE varchar(30),
       OBSERVED integer,
       EXPECTED float
);

\copy profile_segment FROM PROGRAM 'zcat /home/dkulp/data/gpc_wave2_batch1/profile_seq_20_100.seg.txt.gz' DELIMITER E'\t' CSV
\copy profile_counts FROM PROGRAM 'zcat /home/dkulp/data/gpc_wave2_batch1/profile_seq_20_100.pro.txt.gz' DELIMITER E'\t' CSV

--    CREATE UNIQUE INDEX pc_bin_sample_idx ON profile_counts(bin,sample);
CREATE INDEX ps_chrom_start_idx ON profile_segment(chrom,start_pos);
CREATE INDEX ps_chrom_end_idx ON profile_segment(chrom,end_pos);

CREATE UNIQUE INDEX pc_sample_bin_idx ON profile_counts(sample,bin);

CREATE TABLE bkpt (
       sample text,
       chr varchar(2),
       bkpt_bin integer REFERENCES profile_segment(bin),
       gain_ll double precision,
       gain1_ll double precision,
       loss_ll double precision,
       loss1_ll double precision,
       any_ll double precision,
       no_bkpt_ll double precision
);
CREATE UNIQUE INDEX ON bkpt(sample,chr,bkpt_bin);
CREATE INDEX ON bkpt(bkpt_bin);

VACUUM ANALYZE;

