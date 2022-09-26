#!/usr/bin/perl

use strict;
use Getopt::Long;

# Creating Variable
my ($vcfFile,$outputFile,$dataFile,$seqId,$number,$info);

# Creating Arguments
&GetOptions("vcfFile|i=s" => \$vcfFile,
            "outputFile|o=s" => \$outputFile,
	    "dataFile|d=s" => \$dataFile
           );

# Opening outputFile
open(O,">$outputFile")  || die "Unable to open $outputFile";

open(D,"<$dataFile")  || die "Unable to open $dataFile";

# Starting output fasta file with defline. Retrieving total length of sequence.
open(I, "$vcfFile") || die "Unable to open $vcfFile";

while( my $line = <D>) {
    ($seqId, $number) = split(/\t/, $line);
}

close D;

while( my $line = <I>) {
    if ($line =~ /^$number(\t.+)/) {
        $info = $1;
	print O "$seqId$info\n";
    } else {
        print O "$line";
    }
}

close I;
close O;


