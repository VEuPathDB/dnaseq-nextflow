#!/usr/bin/perl

use strict;
use Getopt::Long;

# Creating Variable
my ($inputFile,$outputFile,$sample);

# Creating Arguments
&GetOptions("inputfile|i=s" => \$inputFile,
            "outputFile|o=s" => \$outputFile,
	    "sample|s=s" => \$sample,
           );

# Opening outputFile
open(O,">$outputFile");

# Starting output fasta file with defline. Retrieving total length of sequence.
open(I,"$inputFile") || die "Unable to open $inputFile";

while(my $line = <I>){
    if ($line =~ /^(>)(.*)/) {
        print O "$1$sample\n";
    }
    else {
        print O $line;
    }
}

close I;
close O;
