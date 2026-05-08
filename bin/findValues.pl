#!/usr/bin/perl

use strict;
use Getopt::Long;

my ($vcfFile, $outputFile, $sample);

GetOptions(
    "vcfFile|i=s"  => \$vcfFile,
    "outputFile|o=s" => \$outputFile,
    "sample|s=s"   => \$sample,
);

open(O, ">$outputFile") or die "Cannot open $outputFile: $!";
open(I, "zcat $vcfFile |") or die "Unable to open $vcfFile";

while (<I>) {
    next if /^#/;

    # Require simple (non-symbolic) REF and ALT alleles
    next unless /^(.+)\t(\d+)\t\.\t(\w+)\t(\w+)\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t(\S+)/;

    my ($chrom, $pos, $ref, $alt, $format, $sample_data) = ($1, $2, $3, $4, $5, $6);

    # Parse GT; skip heterozygous calls — they are masked as N×len(REF) in the
    # consensus FASTA and therefore introduce no coordinate shift.
    my @fmt_keys = split(/:/, $format);
    my @fmt_vals = split(/:/, $sample_data);
    my %fmt;
    @fmt{@fmt_keys} = @fmt_vals;
    my $gt = $fmt{GT} // next;

    my @alleles = split(/[\/|]/, $gt);
    my %uniq = map { $_ => 1 } @alleles;
    next if keys(%uniq) > 1;   # heterozygous
    next if $alleles[0] eq '0'; # homozygous ref

    my $shift_count = length($alt) - length($ref);
    print O "$sample\t$chrom\t$pos\t$shift_count\n";
}

close I;
close O;
