#!/usr/bin/env perl
#
# For any bin, set expected and observed depth to zero for all samples when a bin
# has a observed/expected ratio that is more than K units greater than the flanking bins.
#
# Usage: outier_bin_filter.pl K < profileFile > profileFile

use strict;
use FileHandle;

my $K = shift(@ARGV);

if (@ARGV == 0) { unshift(@ARGV,"-") }
my $profileFh = new FileHandle "zcat $ARGV[0] |" or die $!;

sub colratio {
    my ($obs_tot, $exp_tot) = (0,0);
    map { my ($obs, $exp) = split(/,/, $_); $obs_tot += $obs; $exp_tot += $exp } @_;
    if ($exp_tot == 0) {
	return "inf";
    } else {
	return $obs_tot/$exp_tot;
    }
}

#BIN     CHR     START   END     ELENGTH GCRATIO SSC00003        SSC02852        SSC00827        SSC01958        SSC00542        SSC00081        SSC00904        SSC02643        SSC00727        SSC0219
#B1      20      60001   60145   100     79,145  16,4.8890       22,6.1481       20,8.5634       19,6.0003       25,6.4144       27,5.2435       11,5.3249       20,6.4552       22,5.5929       18,5.52
#B2      20      60146   60245   100     54,100  11,4.8648       12,6.1334       26,8.5392       13,5.9841       20,6.3812       16,5.2370       23,5.2914       19,6.4104       8,5.5636        8,5.508
#B3      20      60246   60345   100     54,100  11,4.8560       16,6.1013       16,8.5026       8,5.9604        10,6.3501       13,5.2131       17,5.2884       6,6.3793        17,5.5458       16,5.47
print scalar <$profileFh>;
my ($line, $line_m1, $cr_m1, $cr_m2, $nr);

while ($line = <$profileFh>) {

    my @cols = split(/\t/,$line);

    my $cr = colratio(@cols[6..$#cols]);

    if ($nr++ > 0) {
	if ($nr > 2 and $cr*$K < $cr_m1 and $cr_m2*$K < $cr_m1) {
	    my @pre_cols = split(/\t/,$line_m1);
	    print join("\t",@pre_cols[0..5]).("\t0,0" x ($#cols-5))."\n";
	} else {
	    print $line_m1;
	}
	$cr_m2 = $cr_m1;
    }

    $cr_m1 = $cr;
    $line_m1 = $line;

}
print $line;
