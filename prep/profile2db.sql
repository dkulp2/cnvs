CREATE SCHEMA IF NOT EXISTS ${LABEL};

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

\copy profile_segment FROM PROGRAM 'zcat ${profileFile%%.dat.gz}.seg.txt.gz' DELIMITER E'\t' CSV
\copy profile_counts FROM PROGRAM 'zcat ${profileFile%%.dat.gz}.pro.txt.gz' DELIMITER E'\t' CSV

CREATE INDEX ON profile_segment(chrom,start_pos);
CREATE INDEX ON profile_segment(chrom,end_pos);

CREATE UNIQUE INDEX ON profile_counts(sample,chrom,bin);

VACUUM ANALYZE;

