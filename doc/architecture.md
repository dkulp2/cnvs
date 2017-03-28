# CNV Pipeline Architecture #

## Core Data ##

### Block ###

The primary piece of data is a *block*. A block is conceptually a
subset of basic genomic data, which is a two-dimensional matrix of
chromosome *regions X samples*. Thus, a block corresponds to a
contiguous region on a chromosome and a set of samples of the
region. Each cell in this matrix conceptually contain a sample's DNA.

### Region ###

A chromosome *region* is identified by a genome version, the chromosome
name, a start and end position along the chromosome (without "chr" or
similar prefix), and a bin identifier.  To avoid conflicts with
reserved words, use "startp" and "endp" to refer to positions in the
chromosome.

The interval is inclusive of the start, but not the end position,
i.e. in math notation it's often written
[startp,endp). Or you can think of the indexing as "interbase".

Thus, the region represented by the tuple (hg17,20,100,200) corresponds 
to the 100th through the 199th nucleotides on chromosome 20 of build hg17.
(hg17,20,100,200) and (hg17,20,200,300) are not overlapping and their
lengths are easily computed by subtracting. 

### Bin ###

A *bin* is a region identifier. Chromosomes are partitioned into regions 
of the same "effective" length although the genomic coordinate lengths may vary. 
Effective length refers to masking of low complexity regions where alignments
of sequenced reads tends to be ambiguous. Each bin has approximately the
same entropy.

### Foo ###

  * list item
  * another
  
