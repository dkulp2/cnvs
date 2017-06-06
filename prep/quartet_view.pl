#!/usr/bin/env perl
#
# Combines tables in *a, *b, *c and *d schemas into a single table to be used for web app

use strict;

my $base="data_sfari_batch1";
my $schema="sfari_$ARGV[0]";
my $ext="_$ARGV[0]";

print "CREATE SCHEMA ${schema};\n";

my @tables=qw(bkpt cnv_mle cnv_mrg cnv_post pois posterior posterior_dist prior prior_region profile_counts geno);
my @families=qw(a b c d);

foreach my $T (@tables) {
    print "DROP VIEW IF EXISTS $schema.$T;\n";
    print "CREATE VIEW $schema.$T AS\n";
    print(join("UNION ALL\n", map { "SELECT * FROM $base$_$ext.$T\n" } @families),";\n");
}

print "DROP VIEW IF EXISTS $schema.profile_segment;\n";
print "CREATE VIEW $schema.profile_segment AS SELECT * FROM ${base}a$ext.profile_segment;\n";


