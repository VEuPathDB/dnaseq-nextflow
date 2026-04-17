#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ($vcfFile, $outputFile, $sample);

&GetOptions(
    "vcfFile|i=s"    => \$vcfFile,
    "outputFile|o=s" => \$outputFile,
    "sample|s=s"     => \$sample,
);

die "Usage: findValues.pl -i <vcf> -s <sample> -o <output>\n"
    unless defined $vcfFile && defined $outputFile && defined $sample;

my $cmd = ($vcfFile =~ /\.gz$/) ? "zcat \Q$vcfFile\E" : "cat \Q$vcfFile\E";
open(my $OUT, '>', $outputFile) or die "Cannot open $outputFile: $!";
open(my $IN, "$cmd |") or die "Cannot open $vcfFile: $!";

while (<$IN>) {
    next if /^#/;
    chomp;
    my @f = split(/\t/, $_);
    next unless @f >= 10;

    my $seqId     = $f[0];
    my $pos       = $f[1];
    my $refAllele = $f[3];
    my $altAllele = $f[4];
    my $refLen    = length($refAllele);
    my $altLen    = length($altAllele);
    next if $refLen == $altLen;

    my $shift = $altLen - $refLen;

    my @fmt_keys = split(/:/, $f[8]);
    my @fmt_vals = split(/:/, $f[9]);
    my %fmt;
    @fmt{@fmt_keys} = @fmt_vals;

    my $gt = $fmt{GT} // '.';
    my $zygosity;
    if ($gt =~ m|^(\d+)[/\|](\d+)|) {
        $zygosity = ($1 == $2) ? 'hom' : 'het';
    } else {
        $zygosity = 'hom';
    }

    print $OUT join("\t", $sample, $seqId, $pos, $shift, $zygosity, $refAllele, $altAllele), "\n";
}

close $IN;
close $OUT;
