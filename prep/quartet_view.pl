#!/usr/bin/env perl
#
# Combines tables in *a, *b, *c and *d schemas into a single table to be used for web app

use strict;

my $base="data_sfari_batch1";
my @tables=qw(bkpt cnv_mle cnv_mrg cnv_post pois posterior posterior_dist prior prior_region profile_counts profile_segment);
my @families=qw(a b c d);
foreach my $T (@tables) {
    print "CREATE VIEW $base.$T AS\n";
    print(join("UNION ALL\n", map { "SELECT * FROM $base$_.$T\n" } @families),";\n");
}

