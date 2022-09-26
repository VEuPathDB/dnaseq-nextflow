#!/usr/bin/perl

use strict;
use Getopt::Long;

# Creating Variable
my ($vcfFile,$outputFile,$dataFile,$seqId,$sourceId,$info);

# Creating Arguments
&GetOptions("vcfFile|i=s" => \$vcfFile,
            "outputFile|o=s" => \$outputFile,
	    "dataFile|d=s" => \$dataFile
           );

# Opening outputFile
open(O,">$outputFile")  || die "Unable to open $outputFile";

open(D,">$dataFile")  || die "Unable to open $dataFile";

# Starting output fasta file with defline. Retrieving total length of sequence.
open(I, "$vcfFile") || die "Unable to open $vcfFile";

while( my $line = <I>) {
    if ($line =~ /^(\w+\.\d+)(\t.+)/) {
        $seqId = $1;
        $sourceId = $seqId;
	$info = $2;
	$seqId =~ s/[a-zA-Z\.]//g;
	if ($seqId =~ /^0/) {
	    $seqId = substr($seqId, -1);
	    print O "$seqId$info\n";
	    print D "$sourceId\t$seqId\t";
	}
	else {
            print O "$seqId$info\n";
	    print D "$sourceId\t$seqId\t";
	}
    }
    else {
        print O "$line";
    }
}

close I;
close O;
close D;

