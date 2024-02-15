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

# Starting output fasta file with defline. Retrieving total length of sequence.
open(F,"$faiFile") || die "Unable to open $faiFile";

while(<F>){
    if (/^(\S+)\t(\d+)\t/) {
	&maskSeq($1,$2,$depthCutoff,$pileupFile,$outputFile);
    }
    else {
        die "faiFile not in correct format.";
    }
}

close F;

sub maskSeq {
    my ($seq,$length,$depthCutoff,$pileup,$outputFile) = @_;
    open (P,"$pileup") || die "Unable to open $pileup";
    open(O,">>$outputFile") || die "Unable to open $outputFile";
    my $currentPos = 1;
    my $processedSeq = 0;
    my ($coveredPos,$nuc,$cov);
    while(<P>) {
	# Get Pile up data
        if(/^(\S+)\t(\d+)\t(\w)\t(\d)/){
            if ($1 eq $seq) {
	        print O ">$1\n" if ($processedSeq == 0);
		$processedSeq = 1;
		$coveredPos=$2;
		$nuc = $3;
		$cov = $4;
		until ($currentPos == $coveredPos) {
		    print O "N";
		    $currentPos++;
	        }
	        if ($cov >= $depthCutoff) {
		    print O "$nuc";
		    $currentPos++;
		}
		else {
		    print O "N";
		    $currentPos++;
		}
	    }
	    else {
	        if ($processedSeq == 1) {
		    until ($currentPos > $length) {
		        print O "N";
		        $currentPos++;
		    }
		    last;
	    	}
            }
	}
        else {
	    die "Improper pileup format";
        }
    }
    close P;
    close O;
}
