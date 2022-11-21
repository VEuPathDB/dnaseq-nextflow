#!/usr/bin/perl

use strict;
use Getopt::Long;

# Creating Variable
my ($inputDir);

# Creating Arguments
&GetOptions("inputDir|i=s" => \$inputDir);

opendir(DIR, $inputDir) || die "can't opendir $inputDir: $!";
my @files = grep { /.fa/ } readdir(DIR);
closedir DIR;

my %uniqueids = {};

foreach my $f (@files) {
    my $filePath = "$inputDir/$f";
    open(IN,"$filePath") || die;
    while(<IN>) {
        if (/^>(.+)/) {
	    if (exists $uniqueids{$1}) {
                die "Duplicate defline $1 in $f and $uniqueids{$1}";
	    }
	    else {
                $uniqueids{$1} = $f;
	    }
        }
    }
    close IN;
}

print "No Repeats";
