#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my ($varscanConFile);

&GetOptions("varscanConFile=s"=> \$varscanConFile);

# my $minPval = $minPval ? $minPval : 1e-5;
# my $minCov = $minCov ? $minCov : 5;
# my $minPercent = $minPercent ? $minPercent : 20;

my $strain = $varscanConFile;
$strain =~ s/\.coverage\.txt//g;

my $minPval = 1e-5;
my $minCov = 5;
my $minPercent = 20;

open(my $data, '<', $varscanConFile) || die "Could not open file $varscanConFile: $!";
open(OUT,">output.txt");

my $previousChromosome = '';
my $coverages = '';
my $percents = '';
my @positionArray = ();

while (my $line = <$data>) {
    chomp $line;
    if ($line =~ /^Chrom/) {
	next;
    }
    my ($chrom,$position,$ref,$var,$values,$otherValues,$samplesRef,$samplesHet, $samplesHom, $samplesNC, $lastValues) = split(/\t/, $line);
    my @valueArray = split(":", $values);
    my $cons = $valueArray[0];
    my $cov = $valueArray[1];
    my $percent = $valueArray[4];
    my $pval = $valueArray[5];
    $percent =~ s/%//g;
    if ($percent =~ /-/) {
	if (scalar @positionArray == 0) {
	    $previousChromosome = $chrom;
	    $coverages = '';
	    $percents = '';
	    @positionArray = ();
	    next;
        }
	my $minPos = $positionArray[0];
	my $maxPos = $positionArray[-1];
	print OUT "$previousChromosome\t$minPos\t$maxPos\t$coverages\t$percents\n";
	$previousChromosome = $chrom;
	$coverages = '';
	$percents = '';
	@positionArray = ();
    }
    else {
        my $frequency = '100' - "$percent";
	if($cons eq 'N' || $cov < $minCov || $frequency < $minPercent || $pval < $minPval) {
	    if (scalar @positionArray == 0) {
	        $previousChromosome = $chrom;
	        $coverages = '';
	        $percents = '';
	        @positionArray = ();
	    next;
            }
	    my $minPos = $positionArray[0];
            my $maxPos = $positionArray[-1];
	    print OUT "$previousChromosome\t$minPos\t$maxPos\t$coverages\t$percents\n";
	    $previousChromosome = $chrom;
	    $coverages = '';
	    $percents = '';
	    @positionArray = ();
	}
	else {
	    if (length($coverages) == 0) {
                $coverages = $coverages . "$cov";
	        $percents = $percents . "$frequency";
	    }
	    else {
		$coverages = $coverages . ",$cov";
	        $percents = $percents . ",$frequency";
	    }
	    push(@positionArray, $position);
	    $previousChromosome = $chrom;
	}
    }
}
