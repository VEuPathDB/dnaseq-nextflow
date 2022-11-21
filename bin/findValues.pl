#!/usr/bin/perl

use strict;
use Getopt::Long;

# Creating Variable
my ($vcfFile,$outputFile,$sample,$refPos,$shift_count,$refAllele,$altAllele, $refLen, $altLen);

# Creating Arguments
&GetOptions("vcfFile|i=s" => \$vcfFile,
            "outputFile|o=s" => \$outputFile,
	    "sample|s=s" => \$sample,
           );

# Opening outputFile
open(O,">$outputFile");

# Starting output fasta file with defline. Retrieving total length of sequence.
open(I,"$vcfFile") || die "Unable to open $vcfFile";

while(<I>){
    if (/^(.+)\t(\d+)\t\.\t(\w+)\t(\w+)/) {
        $refPos = $2;
	$refAllele = $3;
	$altAllele = $4;
	$refLen = length($refAllele);
	$altLen = length($altAllele);
	$shift_count = $altLen - $refLen;
	print O "$sample\t$1\t$2\t$shift_count\n";
    }
}

close I;
close O;
