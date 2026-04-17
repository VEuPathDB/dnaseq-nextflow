#!/usr/bin/perl
use strict;
use warnings;
use Test2::V0;
use File::Temp qw(tempfile);

my $vcf_content = <<'VCF';
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	100	.	ATGC	A	50	PASS	DP=30	GT:DP	1/1:30
chr1	200	.	A	ATGC	50	PASS	DP=25	GT:DP	0/1:25
chr1	300	.	A	G	50	PASS	DP=20	GT:DP	1/1:20
VCF

my ($in_fh, $in_file)   = tempfile(SUFFIX => '.vcf');
my ($out_fh, $out_file) = tempfile(SUFFIX => '.tsv');
print $in_fh $vcf_content;
close $in_fh;
close $out_fh;

system("perl bin/findValues.pl -i $in_file -s strainA -o $out_file") == 0
    or die "findValues.pl failed";

open(my $fh, '<', $out_file) or die "Cannot open $out_file";
my @lines = grep { chomp; length($_) } <$fh>;
close $fh;

is(scalar @lines, 2, 'two indels emitted (SNP skipped)');

my @f1 = split(/\t/, $lines[0]);
is($f1[0], 'strainA', 'sample name');
is($f1[1], 'chr1',    'seq_id');
is($f1[2], '100',     'position');
is($f1[3], '-3',      'shift (ATGC->A = -3)');
is($f1[4], 'hom',     'zygosity hom (GT=1/1)');
is($f1[5], 'ATGC',    'ref_allele');
is($f1[6], 'A',       'alt_allele');

my @f2 = split(/\t/, $lines[1]);
is($f2[3], '3',    'shift (A->ATGC = +3)');
is($f2[4], 'het',  'zygosity het (GT=0/1)');
is($f2[5], 'A',    'ref_allele');
is($f2[6], 'ATGC', 'alt_allele');

done_testing();
