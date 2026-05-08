#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Data::Dumper;
use CalculationsForTPM qw(calcTPMForCNVs);

my ($verbose,$geneFootprintFile,$countFile,$tpmFile);
&GetOptions("verbose!"=>\$verbose,
            "geneFootprintFile=s"=> \$geneFootprintFile,
            "countFile=s" => \$countFile,
            "tpmFile=s" => \$tpmFile,
    );

if(!$geneFootprintFile || !$countFile || !$tpmFile){
	die "usage: makeTpmFromHtseqCountsCNV.pl --geneFootprintFile <geneFootprintFile> --countFile <count file from htseq-count> --tpmFile <output file for TPM values>\n";
}

my $geneLengths;
open(IN, "<$geneFootprintFile");
my $line = <IN>;
while ($line=<IN>) {
    chomp($line);
    my ($project, $gene, $length, @rest) = split(/\t/, $line);
    $geneLengths->{$gene} = $length;
}
close(IN);

&calcTPMForCNVs ($geneLengths, $countFile, $tpmFile);

