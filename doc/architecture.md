# CNV Pipeline Architecture #

## Core Data ##

### Block ###

The primary piece of data is a *block*. A block is conceptually a
subset of basic genomic data, which is a two-dimensional matrix of
chromosome *regions X samples*. Thus, a block corresponds to a
contiguous region on a chromosome and a set of samples of the
region. Each cell in this matrix conceptually contain a sample's DNA.

### Region ###

A chromosome *region* is identified by a genome version, the
chromosome name, a start and end position along the chromosome
(without "chr" or similar prefix), and a bin identifier.  To avoid
conflicts with reserved words, use "startp" and "endp" to refer to
positions in the chromosome.

The interval is inclusive of the start, but not the end position. Or
you can think of the indexing as "interbase".

Thus, the region represented by the tuple (hg17,20,100,200)
corresponds to the 100th through the 199th nucleotides on chromosome
20 of build hg17.  (hg17,20,100,200) and (hg17,20,200,300) are not
overlapping and their lengths are easily computed by subtracting.

### Bin ###

A *bin* is a region identifier. Chromosomes are partitioned into bins
of the same "effective" length although the genomic coordinate lengths
may vary.  Effective length refers to masking of low complexity
regions where alignments of sequenced reads tends to be
ambiguous. Each bin has approximately the same entropy.

For convenience finding adjacent features, bins are ordered
integers. (*Is this acceptable?*)

### Segment ###

A *segment* is a generic extent along the genome representing one or
more bins.  A segment is identified by a start and end bin,
i.e. (bin1, bin2). Bins are *inclusive*.  Thus, (bin1, bin1) refers to
a segment that spans a single bin.

### CNV ###

A *CNV* is a copy number variant segment. The CNV has the same copy
number (*CN*) along the segment. For simplicity, sometimes the CN
value is the wildtype, usually 2, even though it is not a variant
value. A CNV is designated as a tuple of a segment and a CN, (bin1,
bin2, CN), e.g.  (bin1=4960,bin2=5120,CN=3).

### Breakpoint ###

A *breakpoint* is the location of a transition between CNVs.  It is
represented by a single bin identifier and corresponds to the region
at the _beginning_ of the bin. For example, given two CNVs (binA,binB)
and (binC,binD), where binC=binB+1, then the breakpoint between these
two CNVs is binC.

Often it's useful to annotate CNVs as a collection of breakpoints
instead of segments. This is because breakpoints are shifted up or
downstream when applying different methods. By using the breakpoint
instead of the segment, the position is changed in one place instead
of changing the end and start position within two segments.  In such
cases, a left (upstream) and right (downstream) copy number is
associated with each breakpoint.

For example, given the CNVs (bin1=4960,bin2=5120,CN=2),
(bin1=5121,bin2=5234,CN=1), and (bin1=5235,bin2=5411,CN=2), then the
CNVs can also be described by two breakpoints at 5121 and
5235. Including the flanking CNs yields the tuples: (bin=5121, CNL=2,
CNR=1) and (bin=5235, CNL=1, CNR=2).

Instead of annotating specific CNs for the left and right of a
breakpoint, sometimes its useful to simply indicate whether the
breakpoint corresponds to an increase or decrease in copy number,
reading from left to right. In the previous example, the second
segment corresponds to a deletion, implying a loss ("L") of copy
number at the first breakpoint followed by a gain ("G") at the
second. Thus, these breakpoints would be annotated as (bin=5121,
change=L) and (bin=5235, change=G).

## Functions on Core Data ##

(Include this? I've had to repeatedly deal with these issues in
different parts of my code. If I were to rewrite then I'd have
dedicated routines for these action.  But this is probably
unnecessary.)

There are some common actions that are often performed on segments and
breakpoints.  These include *transform*, *collapse*, *neighborhood*,
*filter* and *merge*.

*Transform* converts between segments and breakpoints.

*Collapse* creates a maximal segment by overlaying multiple segments
and choosing the smallest start and largest end positions.

*Neighborhood* creates a segment from a set of breakpoints that are
within some distance and meets certain criteria. For example, a
neighborhood might correspond to a region containing "gain"
breakpoints that are less than 10 bins apart.

*Filter* refers to removing CNVs or breakpoints from a working set
because it fails to meet criteria. For example, a CNV may be
eliminated because it is too small, or the CNQ genotype quality values
are too low.  In some cases those segments should be replaced with a
NA value to indicate uncertainty or a region where no call should be
made. In other cases, we want to merge the resulting segments
together.

*Merge* combines two segments into a single segment if it is adjacent
(perhaps within some threshold of distance) and shares a common CN.


These functions may be chained. For example, given the segments
(100,200,CN=2),(201,209,CN=1),(210,300,CN=2), the middle segment may
be filtered out due to its small size and we may want to merge the
resulting segments together to achieve (100,300,CN=2).

## Data Interchange ##

The steps in the CNV pipeline will primarily use a standard input and
output data interchange. The standardized format described here is a
modified delimited format that will allow rapid loading into a
database table using the postgres COPY command, easy loading into R
using read.table, easy text processing, and use of *tabix* for random
access. Support for other databases such as SQLite may be provided in
the future.

A goal of the format is to allow the creation of a SQL table
definition from the metadata. This can also be used in R to coerce
data to the appropriate data types during reading. 

The proposed format is as follows:

* The first line contains a '#' in the first column
* Metadata in the first line is as follows:
  * A table name (without schema) followed by a colon
  * One or more column names, separated by commas and terminated by a semicolon
  * Each column name is optionally followed with an SQL
    data type in curly braces. If missing, then assumed to be TEXT.
  * One or more indices are specified as comma-separated column
    names. Indices within square braces are unique. Indices within
    parenthesis are not.
* Tab-delimited data. There is no quoting. Tabs are not allowed in
  data values.
* NA or NULL values are represnted as missing values, i.e. two
  adjacent tabs. 
   
For example, given a file cnv.txt with the following metadata on the
first line:

    # cnv: chrom{CHAR(2)}, start_bin{INT}, end_bin{INT}, sample{TEXT}, cn{INT}; [chrom,start_bin,sample] (sample)

results in an SQL table definition as:

    CREATE TABLE cnv AS 
	("chrom" CHAR(2), 
     "start_bin" INT, 
     "end_bin" INT,
     "sample" TEXT,
	 "cn" INT);
	
and a table load of:

    \COPY cnv FROM 'cnv.txt' FORMAT CSV DELIMITER E'\t' HEADER
	CREATE UNIQUE INDEX ON cnv(chrom,start_bin,sample);
	CREATE INDEX ON cnv(sample);
	VACUUM ANALYZE cnv;

In R, the metadata above would result in 

    cnv <- as.tibble(read.table('cnv.txt', skip=1, sep='\t', na.strings='',
	                            col.names=c('chrom','start_bin','end_bin','sample','cn'),
	                            colClasses=c('character','integer','integer','character','integer')))

I/O utility programs for SQL called pgRead and pgWrite, written in
Perl, will generate SQL that can be piped to psql. Example usage
is `pgRead cnv.txt | psql`, which will read the tab-delimited file
`cnv.txt`, pass SQL to `psql` that generates the table `cnv`.
`pgWrite cnv > cnv.txt` will read the table `cnv` from the database
and dump its contents with a metadata header to `cnv.txt`.

Two R function, tblRead, tblLoad and tblWrite, are provided in the
'cnv' library. Usage is simply `cnv <- tblRead('cnv.txt')` and
`tblWrite(cnv,'cnv.txt')`. tblLoad is similar to tblRead, but it does
not return a value. Instead it creates an object in the current
environment with the name of the table as specified in the metadata.

There are no helper functions for tabix or generic text processing
with perl, awk, java, etc., but usage is straightforward. For example,
a tabix index can be generated from cnv.txt as:

    bgzip cnv.txt > cnv.txt.gz
	tabix -b 1 -e 3 -s 4 -S 1 cnv.txt.gz

And the file can be queried with 

    tabix cnv.txt.gz 20:1000-2000

Both the Perl/SQL and R routines will automatically decompress an
input file or bgzip an output file that ends in '.gz'. Both tools also
include options to append data.

Currently supported data types map as follows, where types in
parentheses will translate to the default type on the other side of
the '<=>':

* `TEXT (VARCHAR, CHAR) <=> character (factor)`
* `INT <=> integer`
* `FLOAT (DOUBLE PRECISION) <=> numeric`
* `BOOLEAN <=> logical`

### Limitations ###

There is no support for dates or times, yet. 

Tabix requires sorted data, but there is currently no native
understanding of position fields and therefore no support for sorting
output files for tabix.

The read/write functions are lossy with respect to data types because
SQL has a richer set of data types than R. Therefore, when a table is
written to disk, read into R, written to disk and read into postgres,
then some data types are lost - namely fixed width characters.

The read/write functions do not support integrity constraints,
including foreign keys.

The database support is very specific to postgres. Support for other
database systems are possible, but there may be unforeseen
limitations.


