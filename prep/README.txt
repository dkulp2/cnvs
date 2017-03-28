Currently we're working on chromosome 1, GPC wave2.
We're working on just 100 samples.
Bob provided a deletion pipeline output gpc_wave2_batch1/gs_dels.genotypes.vcf.gz
and a CNV pipeline output gpc_wave2/gs_cnv.genotypes.vcf.gz

Bob used a sites.txt file to indicate the windows to use. I used this
to generate maximal ranges per chromosome for selection within the
vcf.

Convert to tab-delimited with

zcat ...vcf.gz | gs_del2tab -s .../sites.txt > ...txt

NOTE: gs_del2tab removes any calls that are FT=LQ. Regarding this, Bob says,
   
"This is fine for now, but for a more thorough analysis, it might also
be worth doing the following:

Keep the LQ sample-segs, but mark them as "copy number uncertain" (and
maybe record the most likely CN estimate?).  When you check overlaps
and take the union, also mark any confidently-discordant intervals as
"copy number uncertain".  That way, you can compute specificity
with/without these "uncertain" intervals.  The uncertain intervals
still have some evidence of copy number variation, they just aren't as
solid or there is conflicting evidence.  In general, I'm probably
happy to treat them leniently (i.e. not to penalize for calling them
and not to penalize for not calling them).

This will only matter, of course, when the calling quality is getting
to be very good."

====================

The known sets are flattened for evaluation purposes by manually running overlap.R.
See doc in file for more info. It creates _flt and _xlft files for evaluation on a
per sample and a per site basis.
Discordant calls are marked with a CN=99.

====================

Load database with profile data by first creating two text files from the profile:

 run profile2db.pl -i /home/dkulp/data/gpc_wave2_batch1/profile_seq_20_100.dat.gz -s /home/dkulp/data/gpc_wave2_batch1/profile_seq_20_100.seg.txt.gz -p /home/dkulp/data/gpc_wave2_batch1/profile_seq_20_100.pro.txt.gz

And then editing profile2db.sql to set the correct filenames and then executing, e.g.

 psql < profile2db.sql

====================
