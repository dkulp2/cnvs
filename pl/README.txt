Data prep:
Retrieve profile*.gz(.tbi) files
Retrieve irs_matrix_OMNI25.dat
Retrieve gs_dels.genotypes.vcf.gz(.tbl)

Create a chrom 20 version of gs_dels with
 tabix2 -o ./gs_dels.chr20.vcf.gz gs_dels.genotypes.vcf.gz 20:1-200000000

???Do I need to create a .bed version of the gs_dels.genotypes.vcf???

Prediction:

Set params in cnv_seg.sh and run


Evaluation:

See header in delSn.sh to prep data
Generate depth data - or remove from code
Set params within and run delSn.sh to generate lots of tallies

Then run DelSiteSnPlot.R to load the selected coverage file (.../tmp/1_coverage.txt) 

-------------------------------------------------------------------------------------------------------------------------------------------------------

Replacement SnSp:

Run SnSp.sh to prepare known data sets
SnSp.sql to load knowns and selected prediction into db and compute
Edges.sql for edged
SnSpPlot.R to read db tables and plot results for both SnSp and Edges

