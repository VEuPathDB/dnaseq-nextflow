#!/usr/bin/perl

use strict;
use warnings;
use VCF;
use Getopt::Long;
use Data::Dumper;

my ($vcfFile, $snpFile);
&GetOptions("vcf=s"=> \$vcfFile,
            "output=s"=> \$snpFile,
           );


my $vcf = VCF->new(file=>$vcfFile);

$vcf->parse_header();

my $firstLoc=$vcf->next_data_hash();
my @strains = keys %{$$firstLoc{gtypes}};

open(SNP, ">$snpFile");

my ($snps, $base, $strain, $percent, $reference, $sequence_source_id, $location, $matches_ref, $quality, $coverage, $pval, $snpId);
while (my $loc=$vcf->next_data_hash()) {
    for my $gt (keys %{$$loc{gtypes}}) {
        my ($alt,$sep,$alt2) = $vcf->parse_alleles($loc,$gt);
	next if ($alt eq '.' && $alt2 eq '.');
	my $snp;
	my $genotype = $$loc{gtypes}->{$gt}->{'GT'};
        if ($genotype eq '1/1' || $genotype eq '1/0') {
            $base = $alt;
        }
        else {
            $base = $alt2;
	}
	$strain = $gt;
	$percent = $$loc{gtypes}->{$gt}->{'FREQ'};
	$percent = $percent =~ s/%//gr;
	$reference = $$loc{REF};
	$sequence_source_id = $$loc{CHROM};
	$location = $$loc{POS};
	$matches_ref = $reference eq $base ? 1 : 0;
	$quality = $$loc{gtypes}->{$gt}->{'GQ'};
	$coverage = $$loc{gtypes}->{$gt}->{'SDP'};
	$pval = $$loc{gtypes}->{$gt}->{'PVAL'};
	$snpId = "NGS_SNP.$sequence_source_id.$location";
	print SNP join("\t", ($sequence_source_id, $location, $strain, $reference, $base, $coverage, $percent, $quality, $pval, $snpId)) . "\n";

    }
}

close SNP;

$vcf->close(); 




