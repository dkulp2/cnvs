# CNV Pipeline Architecture #

## Core Data ##

### Block ###

The primary piece of data is a *block*. A block is conceptually a
subset of basic genomic data, which is a two-dimensional matrix of
chromosome *regions X samples*. Thus, a block corresponds to a
contiguous region on a chromosome and a set of samples. Each cell in
this matrix conceptually contain a sample's DNA.

### Region ###

A chromosome *region* is identified by a genome version, the
chromosome name, a start and end position along the chromosome, and a
bin identifier.

To avoid conflicts with reserved words, use "startp" and "endp" to
refer to positions in the chromosome. The interval is inclusive of the
start, but not the end position. Or you can think of the indexing as
"interbase".

Thus, the region represented by the tuple `(hg17,20,100,200)`
corresponds to the 100th through the 199th nucleotides on chromosome
20 of build hg17.  `(hg17,20,100,200)` and `(hg17,20,200,300)` are not
overlapping and their lengths are easily computed by subtracting.

### Bin ###

A *bin* is a region identifier. Chromosomes are partitioned into bins
of the same "effective" length although the genomic coordinate lengths
may vary.  Effective length refers to masking of low complexity
regions where alignments of sequenced reads tends to be
ambiguous. Each bin has approximately the same entropy.

For convenience finding adjacent features, bins are ordered
integers. (*Is this acceptable? Are bins unique across the genome?
This would make it easier to combine results. But what about gaps
(e.g. centromeres or across chromosomes) where bins may not be
adjacent? For now, I'm assuming that bins are unique within a block
and that merging data will include an additional block ID.*)

### Segment ###

A *segment* is a generic extent along the genome representing one or
more bins.  A segment is identified by a start and end bin,
i.e. `(bin1, bin2)`. Bins are *inclusive*.  Thus, `(bin1, bin1)`
refers to a segment that spans a single bin.

### CNV ###

A *CNV* is a copy number variant segment. The CNV has the same copy
number (*CN*) along the segment. For simplicity, sometimes the CN
value is the wildtype, usually 2, even though it is not a variant
value. A CNV is designated as a tuple of a segment and a CN, `(bin1,
bin2, CN)`, e.g.  `(bin1=4960,bin2=5120,CN=3)`.

### Breakpoint ###

A *breakpoint* is the location of a transition between CNVs.  It is
represented by a single bin identifier and corresponds to the region
at the _beginning_ of the bin. For example, given two CNVs
`(binA,binB)` and `(binC,binD)`, where `binC=binB+1`, then the
breakpoint between these two CNVs is binC.

Often it's useful to annotate CNVs as a collection of breakpoints
instead of segments. This is because breakpoints are shifted up or
downstream when applying different methods. By using the breakpoint
instead of the segment, the position is changed in one place instead
of changing the end and start position within two segments.  In such
cases, a left (upstream) and right (downstream) copy number is
associated with each breakpoint.

For example, given the CNVs `(bin1=4960,bin2=5120,CN=2)`,
`(bin1=5121,bin2=5234,CN=1)`, and `(bin1=5235,bin2=5411,CN=2)`, then the
CNVs can also be described by two breakpoints at 5121 and
5235. Including the flanking CNs yields the tuples: `(bin=5121, CNL=2,
CNR=1)` and `(bin=5235, CNL=1, CNR=2)`.

Instead of annotating specific CNs for the left and right of a
breakpoint, sometimes its useful to simply indicate whether the
breakpoint corresponds to an increase or decrease in copy number,
reading from left to right. In the previous example, the second
segment corresponds to a deletion, implying a loss ("L") of copy
number at the first breakpoint followed by a gain ("G") at the
second. Thus, these breakpoints would be annotated as `(bin=5121,
change=L)` and `(bin=5235, change=G)`.

### Profile ###

A *profile* corresponds to a genomic alignment of read fragments to a
reference assembly for a specific sample. A profile contains a count
of *aligned reads per bin*. 

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
`(100,200,CN=2),(201,209,CN=1),(210,300,CN=2)`, the middle segment may
be filtered out due to its small size and we may want to merge the
resulting segments together to achieve `(100,300,CN=2)`.

## Data Interchange ##

The steps in the CNV pipeline will primarily use a standard input and
output data interchange. The standardized format described here is
tentatively called *metaTbl*. It is a modified delimited format that
will allow rapid loading into a database table using the postgres COPY
command, easy loading into R using read.table, easy text processing,
and use of *tabix* for random access. Support for other databases such
as SQLite may be provided in the future.

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

Three R functions, tblRead, tblLoad and tblWrite, are provided in the
'metaTbl' library. Usage is simply `cnv <- tblRead('cnv.txt')` and
`tblWrite(cnv,'cnv.txt')`. tblLoad is similar to tblRead, but it does
not return a value. Instead it creates an object in the current
environment with the name of the table as specified in the metadata,
similar to load().

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
then some data types are lost -- namely fixed width characters.

The read/write functions do not support integrity constraints,
including foreign keys.

The database support is very specific to postgres. Support for other
database systems is possible, but there may be unforeseen
limitations.


## Pipeline Processes ##

Each *task* in the pipeline operates on a *block*. Most tasks can
operate independently of other blocks, but are typically dependent
on consistent blocks from task to task. For example, the
`ProfileGenotyper` can run on any region, independently, but the set
of samples used affects the output, so subsequent tasks that work on
the `ProfileGenotyper` output should work on the same samples on the
same region.

The expectation is that the input data is sliced into non-intersecting
blocks that are small enough to run a task in a reasonable time, but
are not too small such that the setup and breakdown of a job does not
dominate the running time. Blocks must include a reasonable number of
samples because some statistical operations require multiple samples
to generate distributions. But blocks cannot have too many samples or
a task will take too long to run.

Blocks are non-overlapping. Genome partitioning may be chosen smartly
based on large genomic features, such as centromeres, that are
unlikely to include CNVs of interest. *(Is this OK? Resolving
conflicts from overlapping blocks would be a hard problem to solve
generally and isn't worth the extra work for modest gain, in my
opinion.)*

A further expectation is that tasks are performed on separate
processors in distinct environments with fast local storage and a
slower network file system. It may optionally have access to a local
database and/or a networked database that is shared by multiple
processors.

When a task starts, it reads data from a metaTbl file on shared storage
or from a metaTbl table on a networked database. Similarly, when a job
completes, it writes its results to shared storage or a networked
database as a metaTbl. In this way, each task provides a permanent
checkpoint and merges its results into a larger storage system
containing results from multiple blocks.

To simplify analysis software, each block has a unique schema in the
shared database and a unique directory in the shared file
system. Since blocks are independent, then analysis software and
results can operate in separate database schemas. Visualization and
final collection processes can easily retrieve data from multiple
blocks by creating "meta" schema that contain views of unioned
tables across multiple schemas.

For computational speed, it's likely that multiple metaTbls will
sometimes be copied locally. In particular, it's likely that
multi-table database joins will run faster on local DBMS tables that
are smaller (because they contain only a single block) and are
operating without the contention of a shared DBMS.

Each task may be composed of multiple, sequential *steps*. Multiple
steps can be performed on a single machine using only local file
system and/or local DBMS storage.

*Tasks* and *steps* are both kinds of *jobs*, with the simple
distinction that *steps* are smaller jobs with only local
dependencies, but otherwise their *job descriptions* are the same.

Jobs have data and parameter dependencies. Data dependencies are
requirements that files or tables must be generated for job input.
Parameter dependencies refer to parameter settings in a job
that is used in a previous job.

*I don't know the details of the job description and its queue and
dependency system, yet!*


## CNV Pipeline ##

The CNV pipeline is a set of serial jobs for annotating the
breakpoints of CNVs. The initial input are *profiles* in a
*block*. The ProfileGenotyper estimates genotypes per sample per
bin. The genotypes are aggregated into rolling windows of
bins. Consecutive regions of the same CN are merged into CNVs. The
results are filtered and adjusted. Probability distributions of
breakpoint positions are derived. Breakpoints are further adjusted and
the final predictions reported.

The following is a breakdown of the jobs involved, including their
dependencies and parameters. No effort is made, yet, to group steps
into tasks.

### Job: Block Definition ###

* Parameters:
  * Sample Size: 100
  * Genomic Size: 6e7
  * Reference Genome: hg17
* Dependencies:
  * Reference Genome
* Action: Partition genome
* Output: 
  * `Genome_Regions` (GR_ID, Chrom, StartP, EndP) - physical spans along chromosomes
  * `Sample_Sets` (Set1, sampleA), (Set1, sampleB) - groups of samples processed independently
  * `Block_IDs` (ID1, ID2, etc.) - each ID corresponds to a pair <Sample_Sets, Genome_Regions>

### Job: Bin Definition ###

* Parameters:
  * Effective Bin Size (in nt): 100
* Dependencies:
  * `Genome_Regions`
* Action: Create non-overlapping bins with equal effective length
* Output:
  * `Bin_Map` (GR_ID, Bin, Chrom, StartP, EndP, Effective Length, G+C,
    G+C Total)
  
### Job: Profile Generation ###

* Parameters:
* Dependencies:
  * `Bin_Map`
  * `Sample_Sets`
  * External Profile Data
* Action: Aggregate counts using `Bin_Map`.
* Output: 
  * `Profile_Counts` (_Bin_, _Sample_, Observed, Expected) - summed
    read counts
	
### Job: Window Generation ###

* Parameters:
  * Window Size (bins): 12
* Dependencies:
  * `Bin_Map`
* Action: Create overlapping windows of *Window Size* bins
* Output:
  * `Window_Map` (StartB, EndB) - segments
	
### Job: ProfileGenotyper ###

* Parameters:
* Dependencies:
  * `Profile_Counts`
  * `Sample_Sets`
  * `Window_Map`
* Action: Run org.broadinstitute.sv.apps.ProfileGenotyper on
  `Window_Map` segments using `Profile_Counts` in each segment. Convert
  VCF output into simpler table, one row per sample per window.  Each
  window is identified by the middle bin.
* Outputs:
  * VCF
  * `Geno` (_StartB_, _Sample_, CN, CNQ) - genotype calls and quality per
    sample per segment.

### Job: Merge_CNV ###

* Parameters:
  * Call Percentage Threshold: 0.80 
  * CNQ Threshold: 13 
  * Variance Filter: (W=10, V=0.3, X=3)
  * Span Threshold (bins): 5 
  * NA Span Threshold (bins): 20
* Dependencies:
  * `Geno`
* Actions: 
  * Filter `Geno` segments where fraction of samples calls is below
    *Call Percentage Threshold*.
  * Mask `Geno` segments by setting calls to NA where
    * CNQ is below *CNQ Threshold* or
	* Variance of median CN in W segments is greater than V or
	* Median is greater than X
  * Generate contiguous CNVs as runs of segments with the same CN
  * Remove small CNVs and join new flanking CNVs when (length
    (in bins) of small CNV is less than *Span Threshold* or length of
    small CN=NA CNV is less than *NA Span Threshold*) and flanking
    CNVs have the same CN.
* Output:
  * `CNV_mrg` (_StartB_, _EndB_, _Sample_, CN)
  
### Job: Staircase ###

* Parameters:
  * Mid Segment Size (bins): 12
* Dependencies:
  * `CNV_mrg`
* Action: Remove CNVs that are less than *Mid Segment Size* when
  flanked by CNVs such that CNs are relatively {-1,0,+1} or {+1,0,-1}.
* Output:
  * `CNV_str` (_StartB_, _EndB_, _Sample_, CN)

### Job: Scanning Window Generation ###

* Parameters:
  * Scanning Window Size (bins): 5
* Dependencies:
  * `Bin_Map`
* Action: Create overlapping windows of *Scanning Window Size*
  bins. This is the same as "Job: Window Generation", but a different
  sized window is used for computing joint poissons.
* Output:
  * `Scanning_Window_Map` (StartB, EndB) - segments

### Job: Likelihoods ###

* Parameters:
  * CN Range: (0.1, 1, 2, 3)
* Dependencies:
  * `Profile_Counts`
  * `Sample_Sets`
  * `Scanning_Window_Map`
* Actions: 
  * For each segment per sample, compute the poisson probability for
    the observed counts in `Profile_Counts` given each CN in *CN
    Range*.
  * Combine adjacent pairs of segments (e.g. window of 5+5) to
    generate joint probabilities for each possible combination of CN
    transition. Sum over probabilities to derive gain and loss
    probabilities.
* Output:
  * `Pois` (_StartB_, _Sample_, cnL0.1, cnL1, cnL2, cnL3, cnR0.1, cnR1,
    cnR2, cnR3, Gain, Loss) - pois probabilities

### Job: Maximum Likelihood Estimate of Breakpoint ###

* Parameters:
* Dependencies:
  * `CNV_str`
  * `Pois`
* Actions:
  * For each pair of CNVs in `CNV_str`, setting cnL and cnR to the CN
    of the left and right CNVs, at each possible window of size
    2x*Scanning Window Size* that overlaps the initially annotated
    breakpoints, compute the joint poisson from `Pois` of a transition
    from cnL to cnR. 
  * Normalize across all considered positions. 
  * Compute the maximum likelihood position, the 90% confidence
    interval
  * Store also the number of considered positions and the tail
    probabilities on left and right.
  * Adjust breakpoints to use the MLE.
* Outputs:
  * `CNV_mle` (_StartB_, _EndB_, _Sample_, CN, [CNVID])
  * `Bkpt_mle` (_Bin_, _Sample_, BinRangeL, BinRangeR, BinCIL, BinCIR,
    CNVID\_L, CNVID\_R, LossGain)

### Job: Collapse CNVs ###

* Parameters:
* Dependencies:
  * `CNV_mle`
  * `Bkpt_mle`
* Action:
  * For all gain (loss) breakpoints in `Bkpt_mle`, generate bounds of
    possible gain (loss), by collapsing all overlapping non-wildtype
    CNVs.
* Output:
  * `Bkpt_Range` (_StartB_, _EndB_, [RangeID]) - a range of possible
    breakpoint positions
  * `Bkpt_Range_CNV` (RangeID, CNVID) - the CNVs that contributed to the
    computed range.

### Job: Prior Regions ###

* Parameters:
  * Extend Window Size: 10
* Dependencies:
  * `Bkpt_Range`
  * `Bkpt_Range_CNV`
  * `Bkpt_mle`
  * `Pois`
* Actions:
  * For each candidate breakpoint per sample in `Bkpt_mle`, if there
    is an overlapping `Bkpt_Range`, then for each sample in the range
    (from `Bkpt_Range_CNV`), retrieve the gain (loss) probabilities
    from `Pois` for every position in that range (extending +/-
    *Extend Window Size*).
  * Compute the mean normalized probability of gain (loss) at each
	position in the range.
* Outputs:
  * `Prior` (RangeID, Bin, Loss.u, Gain.u, NC.u)

### Job: Posterior ###

* Parameters:
  * Breakpoint Search Window Size: 10
  * Prior Blend: N/A
* Dependencies:
  * `Bkpt_mle`
  * `Pois`
  * `Prior`
* Actions:
  * For each candidate breakpoint sample in `Bkpt_mle`, retrieve the
    gain (loss) probabilities for that sample in a window of +/-
    *Breakpoint Search Window Size* around position
  * Retrieve the corresponding prior probability for gain (loss) and
    pairwise multiply. (Blend "external" and "internal" priors.)
  * Choose the maximum a posterior position and the 90% confidence
    interval
  * Adjust CNV intervals according to MAP
* Outputs:
  * `Posterior_Dist` (_Bin_, _Sample_, Loss.p, Gain.p, NC.p)
  * `CNV_map` (_StartB_, _EndB_, _Sample_, CN)
