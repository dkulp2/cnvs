After running a prediction set (param settings in cnv.conf and site.conf):

- Execute SnSpSQL.sh, which reloads ALL data - known and predicted.
- Run SnSpPlot.R to generate Sn, Sp and edge plots

One time analysis:

breakPoints.pl/breakPointsPlot.R - generated a profile boxplot of CNQ or O/E around predicted DEL edges
E.g. breakPoints.pl ~/mccarroll/cnv_seg.20.500.90/gs_dels.txt -g ~/mccarroll/cnv_seg.20.500.90/windows.vcf.gz.txt > ~/mccarroll/cnv_seg.20.500.90/breakpoints.txt
     breakPoints.pl ~/mccarroll/cnv_seg.20.500.90/gs_dels.txt -p ~/mccarroll/profiles/profile_seq_20_100.dat.gz > ~/mccarroll/cnv_seg.20.500.90/breakpoints_oe.txt
Then edit breakPointsPlot.R to specify input file

Orphaned:
Edges.sql?
