#/usr/bin/perl

use strict;
use Getopt::Long;

my $file; 
my $out;
my $percentCutoff = 60;
my $coverageCutoff = 5;

&GetOptions("file|f=s" => \$file, 
            "percentCutoff|pc=i"=> \$percentCutoff,
            "coverageCutoff|cc=i"=> \$coverageCutoff,
            "outputFile|o=s"=> \$out,
            );

if (! -e $file || !$out){
die &getUsage();
}

if($file =~ /\.gz$/) {
  open(F, "zcat $file|") || die "unable to open file $file for reading: $!";
}
else {
  open(F, "$file") || die "unable to open file $file for reading: $!";
}


open(O,"|sort -k 1,1 -k 2,2n > $out") or die "Cannot open file for writing: $!";


my ($spanStart, $prevSeq, $prevLoc);

my @coverages;
my @percents;
while(<F>){
  next if /^Chrom\s+Position/;
  chomp;

  my @a = split(/\t/, $_);

  my @valueArray = split(":", $a[4]);

  chop $valueArray[4]; # chop off the % sign
  my $varPercent = 100 - $valueArray[4];
  my $coverage = $valueArray[1];

  my $hasCoverage = $a[3] eq '.' && $varPercent >= $percentCutoff && $coverage >= $coverageCutoff;

  my $isSameSequence = $prevSeq ? $prevSeq eq $a[0]  : 1;

  # start span (doesn't matter what the sequence is 
  if(!$spanStart && $hasCoverage) {
    $spanStart = $a[1];
    push @coverages, $coverage;
    push @percents, $varPercent;
  }
  # end span (sequence transition)
  elsif($spanStart && $hasCoverage && (!$isSameSequence || $prevLoc + 1 != $a[1])) {
    print O "$prevSeq\t$spanStart\t$prevLoc\t" . join(',', @coverages) . "\t" . join(',', @percents) . "\n";
    $spanStart = $a[1];;
    @coverages = ();
    @percents = ();
    push @coverages, $coverage;
    push @percents, $varPercent;
  }
  # end span (no coverge )
  elsif($spanStart && !$hasCoverage) {
    print O "$prevSeq\t$spanStart\t$prevLoc\t" . join(',', @coverages) . "\t" . join(',', @percents) . "\n";
    $spanStart = undef;
    @coverages = ();
    @percents = ();
  }
  elsif($hasCoverage) { 
    push @coverages, $coverage;
    push @percents, $varPercent;
  }
  else {}
  $prevSeq = $a[0];
  $prevLoc = $a[1];
}
if($spanStart) {
  print O "$prevSeq\t$spanStart\t$prevLoc\t" . join(',', @coverages) . "\t" . join(',', @percents) . "\n";
}

close F;
close O;



sub getUsage {
return <<endOfUsage;
parseVarscanToCoverage.pl usage:

  parseVarscanToCoverage.pl --file|f <varscan file> --percentCutoff|pc <frequency percent cutoff [60]> --outputFile|o <output File for coverage> 
endOfUsage
}
