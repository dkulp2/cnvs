#!/usr/bin/env perl
#
# profle2db.pl -i profile.gz -s outseg.txt.gz -p outpro.txt.gz
#
# read compressed profile, write the data to two tables:
# outseg.txt.gz - tab delimited contains BIN, CHR, START, END, ELENGTH, GCRATIO
# outpro.txt.gz - tab delimited contains BIN, SAMPLE, OBSERVED, EXPECTED

use strict;
use FileHandle;
use Getopt::Long;

my $profileFn;
my $outsegFn;
my $outproFn;
my $quiet;

GetOptions("profile|p=s" => \$outproFn,
	   "s|seg=s" => \$outsegFn,
	   "i|input=s" => \$profileFn,
	   "quiet|q" => \$quiet);


die unless defined $outproFn && defined $outsegFn && defined $profileFn;

my $profileFh = new FileHandle "zcat $profileFn |" or die $!;
my $outsegFh;
if ($outsegFn =~ /.gz$/) {
    $outsegFh = new FileHandle "| gzip -c >$outsegFn" or die $!;
} else {
    $outsegFh = new FileHandle ">$outsegFn" or die $!;
}
my $outproFh;
if ($outproFn =~ /.gz$/) {
    $outproFh = new FileHandle "| gzip -c >$outproFn" or die $!;
} else {
    $outproFh = new FileHandle ">$outproFn" or die $!;
}


my $header = <$profileFh>;
chomp $header;
my @hCols = split(/\t/,$header);
shift @hCols foreach (1..6);  # remove first 6 columns

while (<$profileFh>) {
    chomp;
    my ($bin, $chr, $start, $end, $elength, $gcratio, @cols) = split;
    my ($gc, $gctotal) = split(/,/,$gcratio);
    
    $bin =~ s/^B//;
    
    $outsegFh->print("$bin\t$chr\t$start\t$end\t$elength\t$gc\t$gctotal\n");
    foreach my $i (0..$#cols) {
	my ($obs,$exp) = split(/,/,$cols[$i]);
	$outproFh->print("$bin\t$hCols[$i]\t$obs\t$exp\n");
    }
}

