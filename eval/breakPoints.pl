#!/usr/bin/env perl
#
# breakPoints - given a set of prediction (gs_dels.txt), generate
#   metrics for K windows upstream and downstream of CN and CNQ from profileGenotyper or
#   Observed/Expected ratio from raw profile data
#
# breakPoints.pl gs_dels.txt [-g windows.vcf.gz.txt | -p profile_seq_20_100.dat.gz ]
#
# where gs_dels.txt is the output from gs_del2tab
# and windows.vcf.gz.txt is the output of vcf2tab from profileGenotyper
#
# If there's no windows.vcf.gz.txt.tbi, then create a tabix index.
#
# If -g: for each deletion, report the CNQ value for K windows upstream and
# downstream of the left and right edges of the deletion.
#
# Output:
#
# DELETION SAMPLE L/R [-k,...,-1,1,...k] PRED_CN GENO_CN GENO_CNQ DEL_FREQ
#
# If -p: for each deletion, report the expected, observed 
#
# Output:
#
# DELETION SAMPLE L/R [-k,...,-1,1,...k] PRED_CN EXPECTED OBSERVED DEL_FREQ

use strict;
use FileHandle;
use Getopt::Long;

my $genoFn;
my $proFn;

my $minDel = 1200;
my $quiet = 0;

GetOptions("geno|g=s" => \$genoFn,
	   "profile|p=s" => \$proFn,
	   "mindel|m=i" => \$minDel,
	   "quiet|q" => \$quiet)
);

if (defined $genoFn && defined $proFn) {
    die "Can only run with genotype or profile input\n";
}

# create tabix index for genotype table, if it doesn't already exist
if ($genoFn && $genoFn !~ /.gz$/ && !-e "$genoFn.gz.tbi") {
    system("bgzip -c $genoFn > $genoFn.gz");
    system("tabix -s 2 -b 3 -e 4 -S 1 $genoFn.gz");
}
$genoFn = "$genoFn.gz" if defined $genoFn && $genoFn !~ /.gz$/;
# tabix index should already exist for profile


my $inputFn = $genoFn || $proFn;

my $winK = 20; # number of windows upstream and downstream of pos

# determine sample column positions and window increments from vcf2tab genotype file by reading lines 3 and 4
my $inp = new FileHandle "zcat $inputFn|" or die $!;
my $header = <$inp>;
my @cols = split(/\t/,$header);
my $col_i = 0;
my %colMap;
foreach my $col (@cols) {
    $colMap{$col} = $col_i++;
}
my $line = <$inp>;
$line = <$inp>; # line 3
@cols = split(/\t/,$line);
my $start3 = $cols[2];
$line = <$inp>; #line 4
@cols = split(/\t/,$line);
my $start4 = $cols[2];
my $winInc = $start4 - $start3;
$inp->close();

my $dels = new FileHandle $ARGV[0] or die $!;
my @dels; # [ [ sample, id, chrom, start, end, gt, gq, readcount ], [ ... ], .... ]
# suck into memory, ordered by position
my $header = <$dels>;
my %cnvFreq;
while (<$dels>) {
    my @cols = split(/\t/);
    next if $cols[4] - $cols[3] < $minDel;
    push(@dels, [ @cols ]);
    $cnvFreq{$cols[1]}++;
}

my @delsSrt = sort { $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4] } @dels;

# Run tabix for each deletion and collection data.

sub tabix {
    my ($id, $sample, $chr, $pos, $cn, $dir) = @_;
    my $start = $pos - $winK*$winInc;
    my $end = $pos + $winK*$winInc;
    my @rows = `tabix $inputFn $chr:$start-$end`;
    foreach my $row (@rows) {
	my @cols = split(/\t/,$row);
	my $pos2 = $cols[2] + ($cols[3] - $cols[2]) / 2;
	my $offset = int(($pos2 - $pos) / $winInc);

	if (defined $genoFn) {
	    my $geno = $cols[$colMap{$sample}];
	    my $cnq = $cols[$colMap{"${sample}_CNQ"}];
	    print "$id\t$sample\t$dir\t$offset\t$cn\t$geno\t$cnq\t$cnvFreq{$id}\n";
	} else {
	    my ($exp,$obs) = split(/,/,$cols[$colMap{$sample}]);
	    print "$id\t$sample\t$dir\t$offset\t$cn\t$exp\t$obs\t$cnvFreq{$id}\n";
	}
    }
}

my $delIdx = 1;
foreach my $del (@delsSrt) {
    my @fields = @$del;
    my ($sample, $id, $chr, $start, $end, $gt) = @fields[0..5];

    warn sprintf("$delIdx/%d\t%d%%\t$sample\t$id\n",scalar @delsSrt,$delIdx/@delsSrt*100) unless $quiet; $delIdx++;
    tabix($id, $sample, $chr, $start, $gt, "L");
    tabix($id, $sample, $chr, $end, $gt, "R");
}
    
