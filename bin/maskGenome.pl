#!/usr/bin/perl

use strict;
use Getopt::Long;

# Creating Variables
my ($pileupFile,$faiFile,$depthCutoff,$outputFile,$seqlen,$currentPos,$filePos,$nuc,$cov);

# Creating Arguments
&GetOptions("pileupFile|p=s" => \$pileupFile,
            "faiFile|f=s" => \$faiFile, 
            "depthCutoff|dc=i"=> \$depthCutoff,
	    "outputFile|o=s" => \$outputFile,
           );

# Starting output fasta file with defline. Retrieving total length of sequence.
open(F,"$faiFile") || die "Unable to open $faiFile";

while(<F>){
    # SequenceId and total length
    if (/^(\S+)\t(\d+)\t/) {
	&maskSeq($1,$2,$depthCutoff,$pileupFile,$outputFile);
    }
    else {
        die "faiFile not in correct format.";
    }
}

close F;

#========================================================================================================

sub maskSeq {
    my ($seq,$length,$depthCutoff,$pileup,$outputFile) = @_;
    
    open (P,"$pileup") || die "Unable to open $pileup";
    open(O,">>$outputFile") || die "Unable to open $outputFile";
    
    my $currentPos = 1;
    my $processedSeq = 0;
    my $maskedUntilEnd = 0;   
    my ($coveredPos,$nuc,$cov);
    
    while(<P>) {
	# Get Pile up data. SequenceId, Position, Nucleotide, Coverage
        if(/^(\S+)\t(\d+)\t(\w)\t(\d)/){

	    if ($1 eq $seq) {
		# Print Defline
	        print O ">$1\n" if ($processedSeq == 0);

		$processedSeq = 1;
		$coveredPos=$2;
		$nuc = $3;
		$cov = $4;

		# We know we have 0 coverage if the position isn't in the pileup file.
		until ($currentPos == $coveredPos) {
		    print O "N";
		    $currentPos++;
	        }

		# If coverage meets our requirement, print the nucleotide
	        if ($cov >= $depthCutoff) {
		    print O "$nuc";
		    $currentPos++;
		}

		# Not enough coverage, mask
		else {
		    print O "N";
		    $currentPos++;
		}
	    }

	    # If this line is referencing a different sequence
	    else {

		# We have processed a sequence and moved onto a different one in the pileup
	        if ($processedSeq == 1) {
		    
		    # Mask from last position in pileup to total length of sequence. Not in pileup means no coverage
		    until ($currentPos > $length) {
		        print O "N";
		        $currentPos++;
		    }

		    # Used to indicate masking has been complete for this sequence. Needed so we know if we have to still mask to the end for the last sequence in pileup.
                    $maskedUntilEnd = 1;

		    # We have done all processing for the required sequence. Masked sequence has been generated.
		    last;
	    	}
            }

	}

	else {
	    die "Improper pileup format";
        }

    }

    # Masking until end of last sequence. No more lines the in pileup, so we have no coverage, mask.
    if ($processedSeq == 1 && $maskedUntilEnd == 0) {
        until ($currentPos > $length) {
            print O "N";
            $currentPos++;
        }

	# Done with this sequence. Print a new line
        print O "\n";
    }	

    close P;
    close O;
}
