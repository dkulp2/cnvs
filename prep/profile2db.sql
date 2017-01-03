DROP TABLE profile_segment CASCADE;
CREATE TABLE profile_segment (
       BIN integer,
       CHROM varchar(2),
       START_POS integer,
       END_POS integer,
       ELENGTH integer,
       GC integer,
       GCTOTAL integer,
       PRIMARY KEY (CHROM,BIN)
);

DROP TABLE profile_counts;
CREATE TABLE profile_counts (
       BIN integer,
       CHROM varchar(2),
       SAMPLE varchar(30),
       OBSERVED integer,
       EXPECTED float,
       FOREIGN KEY (CHROM,BIN) REFERENCES profile_segment(CHROM,BIN)
);

\copy profile_segment FROM PROGRAM 'zcat /home/dkulp/data/gpc_wave2_batch1/profile_seq_20_100.seg.txt.gz' DELIMITER E'\t' CSV
\copy profile_counts FROM PROGRAM 'zcat /home/dkulp/data/gpc_wave2_batch1/profile_seq_20_100.pro.txt.gz' DELIMITER E'\t' CSV

CREATE INDEX ON profile_segment(chrom,start_pos);
CREATE INDEX ON profile_segment(chrom,end_pos);

CREATE UNIQUE INDEX ON profile_counts(sample,chrom,bin);

CREATE TABLE bkpt (
       label text,
       sample text,
       chr varchar(2),
       bkpt_bin integer,
       gain_ll double precision,
       gain1_ll double precision,
       loss_ll double precision,
       loss1_ll double precision,
       any_ll double precision,
       no_bkpt_ll double precision,
       FOREIGN KEY(chr,bkpt_bin) REFERENCES profile_segment(chrom,bin)
);
CREATE UNIQUE INDEX ON bkpt(sample,chr,bkpt_bin);
CREATE INDEX ON bkpt(bkpt_bin);

VACUUM ANALYZE;

