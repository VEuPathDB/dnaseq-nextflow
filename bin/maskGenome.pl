#!/usr/bin/perl

use strict;
use Getopt::Long;

# Creating Variable
my ($pileupFile,$faiFile,$depthCutoff,$outputFile,$seqlen,$currentPos,$filePos,$nuc,$cov);

# Creating Arguments
&GetOptions("pileupFile|p=s" => \$pileupFile,
            "faiFile|f=s" => \$faiFile, 
            "depthCutoff|dc=i"=> \$depthCutoff,
	    "outputFile|o=s" => \$outputFile,
           );

# Opening outputFile
open(O,">$outputFile");

# Starting output fasta file with defline. Retrieving total length of sequence.
open(F,"$faiFile") || die "Unable to open $faiFile";

while(<F>){
    if (/^(\S+)\t(\d+)/) {
        print O ">$1\n";
        $seqlen = $2;
    }
    else {
        die "faiFile not in correct format.";
    }
}

close F;

# Finding first position in pileupFile
open (P,"$pileupFile") || die "Unable to open $pileupFile";

while(<P>) {
    if(/^(\S+)\t(\d+)/){
        $currentPos=$2;
    }
    last;
}

close P;

# Printing N's for uncovered bases before first covered position
# 2 because this indexes at 1 not zero, and I do not want one for the first covered position
for my $i (2..$currentPos){
  print O "N";
}

# Using pileupFile, currentPos, and seqlen to create masked genome
open (P,"$pileupFile") || die "Unable to open $pileupFile";

while(<P>) {
    if(/^\S+\t(\d+)\t(\w+)\t(\d+)/){
        $filePos = $1;
        $nuc = $2;
        $cov = $3;
        if ($currentPos eq $filePos) {
            if ($cov >= $depthCutoff) {
	        print O $nuc;
	        $currentPos++;  
            }
            else {
	        print O "N";
	        $currentPos++;
            }
	}
        else {
            while($currentPos < $filePos){
	        print O "N";
	        $currentPos++;
            } 
            if ($cov >= $depthCutoff){
	        print O $nuc;
	        $currentPos++;  
            }
            else {
	        print O "N";
	        $currentPos++;  
            }
        }
  }
  else {
      die "$pileupFile in wrong format";
  }
}

# Print N for every uncovered base remaining in genome
while($currentPos ne $seqlen+1){
    print O "N";
    $currentPos++;
}

close P;
close O;
