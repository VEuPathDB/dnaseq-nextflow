#!/usr/bin/perl

use File::Basename;
use Data::Dumper;
use Getopt::Long;
use CBIL::Bio::SequenceUtils;
use Bio::Seq;
use Bio::Tools::GFF;
use Bio::Coordinate::GeneMapper;
use Bio::Coordinate::Pair;
use Bio::Location::Simple;
use Bio::Tools::CodonTable;
use VEuPath::SnpUtils  qw(sampleCacheFileColumnNames snpFileColumnNames alleleFileColumnNames productFileColumnNames);
use VEuPath::MergeSortedSeqVariations;
use ApiCommonData::Load::FileReader;
use locale;
use Sort::Naturally;
use Set::CrossProduct;
use experimental 'smartmatch';

my ($newSampleFile, $cacheFile, $cleanCache, $organismAbbrev, $undoneStrainsFile, $varscanDirectory, $referenceStrain, $help, $debug, $isLegacyVariations, $forcePositionCompute, $consensusFasta, $genomeFasta, $indelFile, $gtfFile);
&GetOptions("new_sample_file=s"=> \$newSampleFile,
            "cache_file=s"=> \$cacheFile,
            "clean_cache"=> \$cleanCache,
            "undone_strains_file=s" => \$undoneStrainsFile,
            "varscan_directory=s" => \$varscanDirectory,
            "is_legacy_variations" => \$isLegacyVariations,
            "organism_abbrev=s" =>\$organismAbbrev,
            "reference_strain=s" => \$referenceStrain,
            "force_position_compute" => \$forcePositionCompute,
            "debug" => \$debug,
            "help|h" => \$help,
	    "consensus=s"=> \$consensusFasta,
	    "genome=s"=> \$genomeFasta,
	    "indelFile=s"=> \$indelFile,
	    "gtfFile=s"=> \$gtfFile
            );

if($help) {
  &usage();
}

my $cacheFileExists = -e $cacheFile;

if(!$cacheFileExists || $cleanCache) {

    if (!$cacheFileExists){
        $|=1; #autoflush 
        print "\033[JCreating empty cache because none exists."."\033[G";        
    }

    if ($cleanCache){
	$|=1; #autoflush 
        print "\033[JCreating empty cache because --clean_cache option was specified"."\033[G";
        open(CACHE, ">$cacheFile") or die "Cannot create a cache file: $!";
        close CACHE;
    }

}

unless(-d $varscanDirectory) {
  &usage("Required Directory Missing") unless($isLegacyVariations);
}

open(UNDONE, $undoneStrainsFile) or die "Cannot open file $undoneStrainsFile for reading: $!";
my @undoneStrains =  map { chomp; $_ } <UNDONE>;
close UNDONE;

my $CODON_TABLE = Bio::Tools::CodonTable->new( -id => 1); #standard codon table

my $dirname = dirname($cacheFile);

my $tempCacheFile = $dirname . "/cache.tmp";
my $snpOutputFile = $dirname . "/snpFeature.dat";
my $alleleOutputFile = $dirname . "/allele.dat";
my $productOutputFile = $dirname . "/product.dat";

my ($snpFh, $cacheFh, $alleleFh, $productFh);
open($cacheFh, "> $tempCacheFile") or die "Cannot open file $tempCacheFile for writing: $!";
open($snpFh, "> $snpOutputFile") or die "Cannot open file $snpOutputFile for writing: $!";
open($alleleFh, "> $alleleOutputFile") or die "Cannot open file $alleleOutputFile for writing: $!";
open($productFh, "> $productOutputFile") or die "Cannot open file $productOutputFile for writing: $!";

my $strainVarscanFileHandles = &openVarscanFiles($varscanDirectory, $isLegacyVariations);

my @allStrains = keys %{$strainVarscanFileHandles};

my $transcriptSummary = &makeTranscriptSummary($gtfFile, $indelFile);

my $currentShifts = &createCurrentShifts($indelFile);

if($forcePositionCompute) {
  push @undoneStrains, $referenceStrain;
}

my $snpFileLine = 1;
my $cacheFileLine = 1;
my ($variations, $variationCount, $greatestLocation, $lastCacheLocation, $lastSnpLocation);

my $cacheFileCount = `wc -l < $cacheFile`;
my $snpFileCount = `wc -l < $newSampleFile`;

chomp $cacheFileCount;
chomp $snpFileCount;

if ($cacheFileCount == 0) {
    $lastCacheLocation = 0;
}

else {
    my $lastCacheLine = `sed "${cacheFileCount}q;d" $cacheFile`;
    my @cacheLineArray = split("\t", $lastCacheLine);
    $lastCacheLocation = $cacheLineArray[1];
}

my $lastSnpLine = `sed "${snpFileCount}q;d" $newSampleFile`;
my @snpLineArray = split("\t", $lastSnpLine);

$lastSnpLocation = $snpLineArray[1];

$greatestLocation = $lastCacheLocation >= $lastSnpLocation ? $lastCacheLocation : $lastSnpLocation;

my $strainFrame;
my $count = 0;
my $lastSnp = 0;

my ($prevSequenceId, $prevTranscriptMaxEnd, $prevTranscripts, $counter, $prevTranscript);

my $i = 1;

while ($i <= $greatestLocation) {
   
  $variations = ();
 
  if ($i == $greatestLocation) {
      $lastSnp = 1;
  }
      
  ($variations, $snpFileLine, $cacheFileLine, $i) = &getVariations($i, $snpFileLine, $cacheFileLine, $newSampleFile, $snpFileCount, $cacheFile, $cacheFileCount, \@undoneStrains);
  
  $i = $lastSnp == 1 ? $greatestLocation+1 : $i;

  next unless ($variations->[0]->{sequence_source_id});

  my ($sequenceId, $location) = &snpLocationFromVariations($variations);

  my ($referenceAllele, $positionsInCds, $positionsInProtein, $referenceVariation, $isCoding);
  
  ($variations, $strainFrame) = &addShiftedLocation($variations, $strainFrame, $currentShifts);

  $variations = &addTranscript($variations, $location, $transcriptSummary);
  
  # clear the transcripts cache once we move past the current transcript
  if($prevTranscript && $variations->[0]->{transcript} && $prevTranscript ne $variations->[0]->{transcript}) {
      # Uncomment to see how this works
      #print "$prevTranscript now onto $variations->[0]->{transcript}\n";
      #print "$transcriptSummary->{$prevTranscript}->{cache}->{ref_cds}\n";
      &cleanCdsCache($transcriptSummary, $prevTranscript, \@allStrains);
      #print "Now is $transcriptSummary->{$previousTranscript}->{cache}->{ref_cds}\n";
  }

  my $transcript = $variations->[0]->{transcript};
  
  if ($transcript) {
      $prevTranscript = $transcript;
  }
  
  my $referenceAllele = $variations->[0]->{reference};
  my $isCoding = $variations->[0]->{is_coding};
    
  my ($refProduct, $refCodon, $refPositionInCds, $refPositionInProtein, $adjacentSnpCausesProductDifference);
  
  if ($isCoding == 1) {
      my $cds_number = $variations->[0]->{cds_number};
      ($refProduct, $refCodon, $refPositionInCds, $refPositionInProtein, $transcriptSummary) = &refProduct($transcriptSummary, $location, $transcript, $genomeFasta, $organismAbbrev, $cds_number);
      ($transcriptSummary, $variations) = &variationProduct($transcriptSummary, $variations, $consensusFasta, $refProduct, $refCodon, $refPositionInCds, $refPositionInProtein);
  }

  foreach my $variation (@$variations) {
    if($variation->{product} ne $refProduct) {
        $adjacentSnpCausesProductDifference = 1;
	$variation->{has_nonsynonomous} = 1;
    }
  }

  $referenceVariation = {'base' => $referenceAllele,
                         'reference' => $referenceAllele,    
                         'location' => $location,
                         'sequence_source_id' => $sequenceId,
                         'matches_reference' => 1,
                         'position_in_cds' => $refPositionInCds,
                         'strain' => $referenceStrain,
                         'product' => $refProduct,
                         'position_in_codon' => $refPositionInProtein,
                         'is_coding' => $isCoding,
                         'has_nonsynonomous' => $adjacentSnpCausesProductDifference,
      			 'ref_codon' => $refCodon
                         };

  push @$variations, $referenceVariation;

  # No need to continue if there is no variation at this point:  Important for when we undo!!
  if(!&hasVariation($variations) && !$isLegacyVariations) {
      next;
  }

  # loop over all strains add coverage vars
  my @variationStrains = map { $_->{strain} } @$variations;
    
  unless($isLegacyVariations) {
      my $coverageVariations = &makeCoverageVariations(\@allStrains, \@variationStrains, $strainVarscanFileHandles, $referenceVariation);
      my @coverageVariationStrains = map { $_->{strain} } @$coverageVariations;
      print STDERR "HAS COVERAGE VARIATIONS FOR THE FOLLOWING:  " . join(",", @coverageVariationStrains) . "\n" if($debug);
      push @$variations, @$coverageVariations;
  }
    
  # loop through variations and print
  foreach my $variation (@$variations) {

      # Adding a field to show if snp is downstream of a frameshift
      if ($variation->{transcript}) {
	  my $downstreamOfFrameshift = &checkForDownstreamOfFrameshift($transcriptSummary, $variation);
	  $variation->{downstream_of_frameshift} = $downstreamOfFrameshift;
      }
      else {
          $variation->{downstream_of_frameshift} = 0;
      }
      
      if($variation->{strain} eq $referenceStrain) {
          next;
      }

      $variation->{matches_reference} = ($variation->{base} eq $referenceAllele) ? 1 : 0; 
      
      if ($variation->{strain} ~~ @variationStrains) {
	  &printVariation($variation, $cacheFh);
      }   

  }

  my $snpFeature = &makeSNPFeatureFromVariations($variations, $referenceVariation);

  #print Dumper $snpFeature;
  
  &printSNPFeature($snpFeature, $snpFh);

  if ($referenceVariation->{is_coding} == 1) {
      my $productFeature = &makeProductFeatureFromVariations($variations, $referenceVariation);
      my $alleleFeature = &makeAlleleFeatureFromVariations($variations);
      &printProductFeature($productFeature, $productFh);
      &printAlleleFeature($alleleFeature, $alleleFh);
  }

  $count += 1;
  if ($count % 1000 == 0) {
      print "Processed $count Snps\n";
  }

}

close $cacheFh;
close $snpFh;
close $alleleFh;
close $productFh;
close OUT;
&closeVarscanFiles($strainVarscanFileHandles);

# compare file sizes of old and new cache file                                                                                                                                                             
my $newCacheCount = `cat $tempCacheFile | wc -l`;
chomp($newCacheCount);

# Rename the output file to full cache file                                                                                                                                                                
unlink $cacheFile or warn "Could not unlink $cacheFile: $!";
rename $tempCacheFile, $cacheFile;

# overwrite existing sample file w/ empty file unless we are skipping coverage                                                                                                                             
open(TRUNCATE, ">$newSampleFile") or die "Cannot open file $newSampleFile for writing: $!";
close(TRUNCATE);

# overwrite existing UndoneStrains file w/ empty file                                                                                                                                                      
open(TRUNCATE, ">$undoneStrainsFile") or die "Cannot open file $undoneStrainsFile for writing: $!";
close(TRUNCATE);

#--------------------------------------------------------------------------------
# BEGIN SUBROUTINES
#--------------------------------------------------------------------------------

sub usage {
  my ($m) = @_;

  if($m) {
    print STDERR $m . "\n";
    die "Error running program";
  }

  print STDERR "usage:  processSequenceVariations.pl --new_sample_file=<FILE> --cache_file=<FILE> [--gusConfigFile=<GUS_CONFIG>] --undone_strains_file=<FILE> --varscan_directory=<DIR> --transcript_extdb_spec=s --organism_abbrev=s --reference_strain=s\n";
  exit(0);
}


sub printVariation {
  my ($variation, $fh) = @_;
  my $keys = &sampleCacheFileColumnNames();
  my $products = $variation->{product};
  my $productString = "";
  foreach $product (@$products) {
      $productString = $productString . "$product:";
  }
  chop($productString);
  $variation->{product} = $productString;
  print $fh join("\t", map {$variation->{$_}} @$keys) . "\n";
}


sub makeCoverageVariations {
  my ($allStrains, $variationStrains, $strainVarscanFileHandles, $referenceVariation) = @_;

  my @rv;

  foreach my $strain (@$allStrains) {
    my $hasVariation;

    foreach my $varStrain (@$variationStrains) {
      if($varStrain eq $strain) {
        $hasVariation = 1;
        last;
      }
    }

    unless($hasVariation) {
      my $fileReader = $strainVarscanFileHandles->{$strain} ;

      my $variation = &makeCoverageVariation($fileReader, $referenceVariation, $strain);

      if($variation) {
        push @rv, $variation;
      }

    }

  }
  return \@rv;
}


sub makeCoverageVariation {
  my ($fileReader, $referenceVariation, $strain) = @_;

  my $rv;

  my $location = $referenceVariation->{location};
  my $sequenceId = $referenceVariation->{sequence_source_id};
  my $referenceAllele = $referenceVariation->{base};

  while($fileReader->hasNext()) {
    # look at the line in memory to see if my sequence and location are inside;  if so, then last
    my @p = $fileReader->getPeek();
    my $pSequenceId = $p[0];
    my $pStart = $p[1];
    my $pEnd = $p[2];

    if($pSequenceId eq $sequenceId && $location >= $pStart && $location <= $pEnd) {

      unless(defined($fileReader->{_coverage_array})) {
        my @coverageArray = split(",", $p[3]);
        my @percentsArray = split(",", $p[4]);
        
        $fileReader->{_coverage_array} = \@coverageArray;
        $fileReader->{_percents_array} = \@percentsArray;
      }

      my $index = $location - $pStart;

      $rv = {'base' => $referenceAllele,
             'location' => $location,
             'sequence_source_id' => $sequenceId,
             'matches_reference' => 1,
             'strain' => $strain,
             'coverage' => $fileReader->{_coverage_array}->[$index],
             'percent' => $fileReader->{_percents_array}->[$index],

      };
      last;
    }

    # stop when the location from the line in memory is > the refLoc
    if($pSequenceId gt $sequenceId || ($pSequenceId eq $sequenceId && $pStart > $location)) {
      last;
    }

    # read the next line into memory
    $fileReader->nextLine();
    $fileReader->{_coverage_array} = undef;
    $fileReader->{_percents_array} = undef;
  }

  return $rv;
}

sub hasVariation {
  my ($variations) = @_;

  my %alleles;
  foreach(@$variations) {
    my $base = $_->{base};
    $alleles{$base}++;
  }

  return scalar(keys(%alleles)) > 1;
}

sub snpLocationFromVariations {
  my ($variations) = @_;

  my $sequenceIdRv;
  my $locationRv;

  foreach(@$variations) {
    my $sequenceSourceId = $_->{sequence_source_id};
    my $location = $_->{location};

    die "sequenceSourceId and location required for every variation" unless($sequenceSourceId && $location);

    if(($sequenceIdRv && $sequenceIdRv ne $sequenceSourceId) || ($locationRv && $locationRv != $location)) {
      print STDERR Dumper $variations;
      die "Multiple variation locations found for a snp";
    }

    $sequenceIdRv = $sequenceSourceId;
    $locationRv = $location;
  }
  return($sequenceIdRv, $locationRv);
}


sub closeVarscanFiles {
  my ($fhHash) = @_;

  foreach(keys %$fhHash) {
    $fhHash->{$_}->closeFileHandle();
  }
}

sub openVarscanFiles {
  my ($varscanDirectory, $isLegacyVariations) = @_;

  my %rv;

  return \%rv if($isLegacyVariations);

  opendir(DIR, $varscanDirectory) or die "Cannot open directory $varscanDirectory for reading: $!";

  while(my $file = readdir(DIR)) {
    my $reader;
    my $fullPath = $varscanDirectory . "/$file";

    if($file =~ /(.+)\.coverage\.txt$/) {
      my $strain = $1;

      if($file =~ /\.gz$/) {
        print STDERR "OPEN GZ FILE: $file for Strain $strain\n" if($debug);

        $reader = ApiCommonData::Load::FileReader->new("zcat $fullPath |", [], qr/\t/);
    } 
      else {
        $reader = ApiCommonData::Load::FileReader->new($fullPath, [], qr/\t/);
      }

      $rv{$strain} = $reader;
    }
  }

  return \%rv;
}


sub cleanCdsCache {
  my ($transcriptSummary, $prevTranscript, $allStrains) = @_;  
  $transcriptSummary->{$prevTranscript}->{cache}->{ref_cds} = undef;
  foreach my $strain (@$allStrains) {
      #print "$transcriptSummary->{$prevTranscript}->{cache}->{$strain}->{consensus_cds}\n";
      $transcriptSummary->{$prevTranscript}->{cache}->{$strain}->{consensus_cds} = undef;
      #print "Now is $transcriptSummary->{$prevTranscript}->{cache}->{$strain}->{consensus_cds}\n";
  }
}

sub getValues {
    my ($file) = @_;
    my @values = ();
    open(my $data, '<', $file) || die "Could not open file $file: $!";
    my $counter = 0;
    while (my $line = <$data>) {
        chomp $line;
        my ($seqName, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, $line);
	my ($transcriptId, $geneId, $geneName) = split(/;/, $attribute);
        $transcriptId = $transcriptId =~ s/transcript_id //gr;
	$transcriptId = $transcriptId =~ s/"//gr;
        $geneId = $geneId =~ s/gene_id //gr;
        $geneId = $geneId =~ s/"//gr;		
        $geneName = $geneName =~ s/gene_name //gr;
        $geneName = $geneName =~ s/"//gr;		
	push ( @{$values[$counter]}, ($seqName, $feature, $start, $end, $strand, $transcriptId, $geneId, $geneName));
	$counter++;
    }
    return @values;
}

sub makeTranscriptSummary {
    my ($gtfFile,$indelFile) = @_;
    my $currentShifts = &createCurrentShifts($indelFile);
    my @values = &getValues($gtfFile);
    my %transcriptSummary;
    my $valueLen = @values;
    my $valueLenIndex = $valueLen - 1;
    my $currentFrame = 0;
    my ($seqName, $strand, $transcriptId, $geneId, $geneName, $cdsFrame, $cdsStart, $cdsEnd, $cdsStrand);
    until ($currentFrame > $valueLenIndex) {
        my $currentTranscriptId = $values[$currentFrame][5];
        my $cdsCounter = 1;
	while ($values[$currentFrame][5] eq $currentTranscriptId) {
            $seqName = $values[$currentFrame][0];
            $strand = $values[$currentFrame][4];
            $transcriptId = $values[$currentFrame][5];
            $geneId = $values[$currentFrame][6];
            $geneName = $values[$currentFrame][7];
            if ($strand eq "-") {
                $strand = -1;
            }
            else {
                $strand = 1;
            }
	    
            if ($values[$currentFrame][1] eq "CDS") {
                my $cdsStart = $values[$currentFrame][2];
                my $cdsEnd = $values[$currentFrame][3];
		my $cdsStartField = "cds_start_$cdsCounter";
		my $cdsEndField = "cds_end_$cdsCounter";
	        $transcriptSummary{$transcriptId}->{$cdsStartField} = $cdsStart;
                $transcriptSummary{$transcriptId}->{$cdsEndField} = $cdsEnd;
                $cdsCounter++;
            }
	    $transcriptSummary{$transcriptId}->{sequence_source_id} = $seqName;
            $transcriptSummary{$transcriptId}->{transcript_id} = $transcriptId;
            $transcriptSummary{$transcriptId}->{cds_strand} = $strand;
            $currentFrame++;
        }
    }
    my $transcriptSummaryShifted = &addStrainCDSShiftsToTranscriptSummary($currentShifts, \%transcriptSummary);
    my $transcriptSummaryWithFrameShifts = &addFrameShifts($currentShifts, $transcriptSummaryShifted);
    #my $transcriptSummaryWithTranscriptomicLocations = &addTranscriptomicLocation($transcriptSummaryShifted, $currentShifts);
    return $transcriptSummaryWithFrameShifts;
}

sub calculateAminoAcidPosition {
  my ($codingPosition) = @_;

  my $aaPos = ($codingPosition % 3 == 0) ? int($codingPosition % 3) + 1 : int($codingPosition % 3);

  return($aaPos);
}


sub getAminoAcidSequenceOfSnp {
  my ($cdsSequence, $positionInCds) = @_;
  my $codon;

  my $positionInProtein = &calculateAminoAcidPosition($positionInCds);
  
  if ($positionInProtein == 1) {
      $codon = substr($cdsSequence, $positionInCds-1, 3);
  }
  elsif ($positionInProtein == 2) {
      $codon = substr($cdsSequence, $positionInCds-2, 3);
  }
  elsif ($positionInProtein == 3) {
      $codon = substr($cdsSequence, $positionInCds-3, 3);
  }
  else {
      print "$positionInProtein not in range\n";
      die;
  }
  
  my $codons = &calculatePossibleCodons($codon);

  my $products;
  my $productsLen = scalar @$codons;
  $productsLen=$productsLen-1;

  foreach my $i (0..$productsLen) {
      my $current_codon = $codons->[$i];
      my $product = $CODON_TABLE->translate($current_codon);
      push @{ $products }, $product; 
  }

  return $codon, $products;
}

sub createCurrentShifts {
    my ($file) = @_;
    my @strains = &getDistinctStrains($file);
    my $currentShifts;
    foreach my $strain (@strains) {
	my @locationshifts = ();
        my $counter = 0;
        my $currentShift = 0;
	my $sourceId;
        open(my $data, '<', $file) || die "Could not open file $file: $!";
	while (my $line = <$data>) {
	    chomp $line;
	    my ($name, $seqId, $refpos, $shift) = split(/\t/, $line);
	    if ( $name eq $strain ) {
		push ( @{$locationshifts[$counter]}, ($refpos, $shift + $currentShift));
		$counter++;
		$sourceId = $seqId;
                $currentShift = $shift + $currentShift;
	    }	
	}
	push @{ $currentShifts->{$strain}->{$sourceId}}, \@locationshifts;
    }
    return $currentShifts;
}

sub getDistinctStrains {
    my ($file) = @_;
    my @strains = ();
    open(my $data, '<', $file) || die "Could not open file $file: $!";
    while (my $line = <$data>) {
        chomp $line;
        my ($name, $seqId, $refpos, $shift) = split(/\t/, $line);
	if (!grep(/^$name$/,@strains)) {
            push (@strains, $name);
        }
    }
    return @strains;
}


sub addStrainCDSShiftsToTranscriptSummary {
    my ($currentShifts, $transcriptSummary) = @_;
    my ($oldShift, $shiftFrame);

    foreach my $strain (keys %{ $currentShifts }) {

	my $shiftArray = $currentShifts->{$strain};

	foreach my $chromosome (keys %{ $shiftArray }) {	   

	    my @chromosomeShiftArray = $shiftArray->{$chromosome};
	    my $indexedArray = $chromosomeShiftArray[0][0];
	    my $shiftArrayLen = scalar @{ $indexedArray };
	    my $shiftFrameLimit = $shiftArrayLen - 1;
	    my $oldCdsShift = 0;
	    my $cdsShiftFrame = 0;
	    my ($cds_start, $cds_end);
            my $startIndicator = "start";
            my $endIndicator = "end";
	    my @sorted_keys = nsort keys %{ $transcriptSummary };

	    foreach my $transcript (@sorted_keys) {

		my ($shifted_cds_start, $shifted_cds_end);
		my $cds_count = 0;

		foreach my $key (keys %{ $$transcriptSummary{$transcript} }) {
		    if ($key =~ /cds/) {
			$cds_count++;
		    }
		}

		$cds_count = $cds_count/2;
		
		foreach my $number (1 .. $cds_count) {

		    my $cdsStartField = "cds_start_$number";
		    my $cdsEndField = "cds_end_$number";
		    my $cdsShiftedStartField = "cds_shifted_start_$number";
		    my $cdsShiftedEndField = "cds_shifted_end_$number";
                    		    
		    if ($$transcriptSummary{$transcript}->{sequence_source_id} !~ /$chromosome/) {
		        if ($$transcriptSummary{$transcript}->{$strain}->{$cdsShiftedStartField}) {
                            next;
                        }
                        else {
                            $$transcriptSummary{$transcript}->{$strain}->{$cdsShiftedStartField} = $$transcriptSummary{$transcript}->{$cdsStartField};
                            $$transcriptSummary{$transcript}->{$strain}->{$cdsShiftedEndField} = $$transcriptSummary{$transcript}->{$cdsEndField};
                            next;
                        }
		    }

	            $cds_start = $$transcriptSummary{$transcript}->{$cdsStartField};
                    $cds_end = $$transcriptSummary{$transcript}->{$cdsEndField};	    

		    ($shifted_cds_start, $cdsShiftFrame, $oldCdsShift) = &calcCoordinates($cdsShiftFrame, $shiftFrameLimit, $oldCdsShift, $cds_start, $startIndicator, \@chromosomeShiftArray);
		    ($shifted_cds_end, $cdsShiftFrame, $oldCdsShift) = &calcCoordinates($cdsShiftFrame, $shiftFrameLimit, $oldCdsShift, $cds_end, $endIndicator, \@chromosomeShiftArray);

                    $$transcriptSummary{$transcript}->{$strain}->{$cdsShiftedStartField} = $shifted_cds_start;
		    $$transcriptSummary{$transcript}->{$strain}->{$cdsShiftedEndField} = $shifted_cds_end;

               }		
	    }
        }
    }
    return $transcriptSummary;
}


sub calcCoordinates {
    my ($shiftFrame, $shiftFrameLimit, $oldShift, $coordinate, $indicator, $shiftArray) = @_;
    my $shiftedLocation;
    my $oldFrame;

    if ($coordinate < $shiftArray->[0][0][$shiftFrame][0]) {
	$shiftedLocation = $oldShift + $coordinate;
    }

    elsif ($shiftArray->[0][0][$shiftFrame][0] == $coordinate) {
        my $currentShift = $shiftArray->[0][$shiftFrame][1];

	if ($currentShift == 0) {
	    $shiftedLocation = $coordinate;
        }

	elsif ($indicator eq 'start') {
	    $shiftedLocation = $oldShift + $coordinate;
        }

	elsif ($indicator eq 'end' && $currentShift > 0) {
	    $shiftedLocation = $currentShift + $coordinate;   
        }

	elsif ($indicator eq 'end' && $currentShift < 0) {
	    $shiftedLocation = $oldShift + $coordinate;     
        }
    }

    elsif ($coordinate > $shiftArray->[0][0][$shiftFrame][0] || $shiftFrame == $shiftFrameLimit) {

	until ($shiftArray->[0][0][$shiftFrame][0] >= $coordinate || $shiftFrame == $shiftFrameLimit) {
	    $oldShift = $shiftArray->[0][0][$shiftFrame][1];
	    $shiftFrame++;
	}

	if ($shiftFrame == $shiftFrameLimit && $coordinate < $shiftArray->[0][0][$shiftFrame][0]) {
	    $shiftedLocation = $coordinate + $shiftArray->[0][0][$shiftFrame-1][1];
	}

	elsif ($shiftFrame == $shiftFrameLimit && $coordinate > $shiftArray->[0][0][$shiftFrame][0]) {
	    $shiftedLocation = $coordinate + $shiftArray->[0][0][$shiftFrame][1];
	}

	elsif ($shiftArray->[0][0][$shiftFrame][0] == $coordinate) {

	    if ($shiftArray->[0][0][$shiftFrame][1] == 0) {
                $shiftedLocation = $coordinate;
            }

	    elsif ($indicator eq 'start') {
		$oldFrame = $shiftFrame - 1;
                $shiftedLocation = $shiftArray->[0][0][$oldFrame][1] + $coordinate;
            }

	    elsif ($indicator eq 'end' && $shiftArray->[0][0][$shiftFrame][1] > 0) {
                $shiftedLocation = $shiftArray->[0][0][$shiftFrame][1] + $coordinate;     
            }

	    elsif ($indicator eq 'end' && $shiftArray->[0][0][$shiftFrame][1] < 0) {
                $oldFrame = $shiftFrame - 1;
                $shiftedLocation = $shiftArray->[0][0][$oldFrame][1] + $coordinate;     
            }
	}

	else {
	    $shiftedLocation = $oldShift + $coordinate;
	}
    }
    return ($shiftedLocation, $shiftFrame, $oldShift);   
}


sub addShiftedLocation {
    my ($variations, $strainFrame, $currentShifts) = @_;
    my ($location, $strain, $indexedArray, $shiftArrayLen, $shiftFrameLimit, $oldShift, $shiftFrame, $shiftedLocation, $chromosome);

    foreach my $variation (@$variations) {

	$location = $variation->{location};
        $strain = $variation->{strain};
	$chromosome = $variation->{sequence_source_id};
        my @shiftArray = $currentShifts->{$strain}->{$chromosome};
        $indexedArray = $shiftArray[0][0];
        $shiftArrayLen = scalar @{ $indexedArray };
        $shiftFrameLimit = $shiftArrayLen - 1;

	if ($currentShifts->{$strain}->{$chromosome}) {

	    if ($strainFrame->{$strain}->{$chromosome}->{shiftFrame}) {
                $shiftFrame = $strainFrame->{$strain}->{$chromosome}->{shiftFrame};
                $oldShift = $strainFrame->{$strain}->{$chromosome}->{oldShift};
            }

	    else {
                $shiftFrame = 0;
                $oldShift = 0;
            }

	    until ($shiftArray[0][0][$shiftFrame][0] >= $location || $shiftFrame == $shiftFrameLimit) {
                $oldShift = $shiftArray[0][0][$shiftFrame][1];
                $shiftFrame++;
            }

	    if ($shiftFrame == $shiftFrameLimit && $location <= $shiftArray[0][0][$shiftFrame][0]) {
                $shiftedLocation = $location + $shiftArray[0][0][$shiftFrame-1][1];
            }

	    elsif ($shiftFrame == $shiftFrameLimit && $location > $shiftArray[0][0][$shiftFrame][0]) {
                $shiftedLocation = $location + $oldShift;
            }

	    elsif ($location <= $shiftArray[0][0][$shiftFrame][0]) {
                $shiftedLocation = $location + $oldShift;
	    }

	    else {
                $shiftedLocation = $location + $oldShift;
            }
	}

	else {
            $shiftedLocation = $location;
	    $oldShift = 0;
	}

	$variation->{shifted_location} = $shiftedLocation;
        $variation->{current_shift} = $oldShift;
        $strainFrame->{$strain}->{$chromosome}->{oldShift} = $oldShift;
        $strainFrame->{$strain}->{$chromosome}->{shiftFrame} = $shiftFrame;
    }
    return ($variations, $strainFrame);
}


sub getCodingSequence {
    my ($defline, $start, $end, $strand, $fasta) = @_;
    my $seq;
    $seq = `samtools faidx $fasta $defline:$start-$end`;
    $seq = $seq =~ s/>.+\n//gr;
    $seq = $seq =~ s/\n//gr;
    if ($strand == -1) {
	$seq = &reverse_complement_IUPAC($seq)
    }
    return $seq;
}


sub calculatePossibleCodons {
    my ($codon) = @_;
    my $codonList;

    if (length($codon) != 3) {
	$codon = "NNN";
    }
    $codon = uc($codon);
    my @codonArray=split(//, $codon);

    my %translate = (A => ['A'],
    		     G => ['G'],
		     C => ['C'],
		     T => ['T'],
                     R => ['A','G'],
		     Y => ['C','T'],
		     K => ['G','T'],
		     M => ['A','C'],
		     S => ['G','C'],
		     W => ['A','T'],
		     B => ['G','T','C'],
		     D => ['G','A','T'],
		     H => ['A','C','T'],
		     V => ['G','C','A'],
		     N => ['A','G','C','T']
                     );
    my @expanded = map { $translate{$_} } @codonArray;

    my $iterator = Set::CrossProduct->new(\@expanded);

    foreach my $codon ($iterator->combinations) {
        my $string = join(",", @$codon);
        $string = $string =~ s/,//gr;
        push @{ $codonList }, $string;
    }

    return $codonList;
}


sub printSNPFeature {
    my ($snp, $snpFh) = @_;
    my $keys = VEuPath::SnpUtils::snpFileColumnNames();
    print $snpFh join("\t", map {$snp->{$_}} @$keys) . "\n";
}

sub printProductFeature {
    my ($products, $productFh)= @_;
    my $keys = VEuPath::SnpUtils::productFileColumnNames();
    foreach my $product ($products) {
	print $productFh join("\t", map {$product->[0]->{$_}} @$keys) . "\n";
    }
}

sub printAlleleFeature {
    my ($alleles, $alleleFh)= @_;
    my $keys = VEuPath::SnpUtils::alleleFileColumnNames();
    foreach my $allele ($alleles) {	
	print $alleleFh join("\t", map {$allele->[0]->{$_}} @$keys) . "\n";
    }
}


sub makeAlleleFeatureFromVariations {
  my ($variations) = @_;
  my %alleleCounts;
  my %strains;
  my $alleles;
  foreach my $variation (@$variations) {
    my $allele = $variation->{base};
    $alleleCounts{$allele} ++;
  }
  foreach my $allele (keys %alleleCounts) {
      my $distinctStrainCount;
      my %strains;
      my $count;
      my $avg_read_percent=0;
      my $avg_coverage=0;
      foreach my $variation (@$variations) {
	  next unless($variation->{base} eq $allele);
	      $count++;
	      my $strain = $variation->{strain};
	      $strains{$strain}++;
	      my $percent = $variation->{percent};
	      my $coverage = $variation->{coverage};
	      $avg_read_percent += $percent;
	      $avg_coverage += $coverage;
      }
      $distinctStrainCount = scalar keys %strains;
      $avg_read_percent = $avg_read_percent / $count;
      my $rounded_avg_read_percent = sprintf("%.2f", $avg_read_percent);
      $avg_coverage = $avg_coverage / $count;
      my $rounded_avg_coverage = sprintf("%.2f", $avg_coverage);
      my $all = { "allele" => $allele,
                  "distinct_strain_count" => $distinctStrainCount,
	          "allele_count" => $count,
	          "average_coverage" => $rounded_avg_coverage,
	          "average_read_percent" => $rounded_avg_read_percent
                };
      push @$alleles, $all;
  }
  return $alleles;
}


sub makeProductFeatureFromVariations {
  my ($variations, $referenceVariation) = @_;
  my %productCounts;
  my $products;
  my $refLocationCds = $referenceVariation->{position_in_cds};
  my $refLocationProtein = $referenceVariation->{position_in_codon};
  foreach my $variation (@$variations) {
      next unless($variation->{product});
      my $productsArray = $variation->{product};
      my $productsLen = scalar @$productsArray;
      $productsLen = $productsLen-1;
      foreach my $i (0..$productsLen) {
          my $product = $productsArray->[$i];
	  $productCounts{$product}++;
      }
  }
  foreach my $variation (@$variations) {
      next unless($variation->{codon});
      my $transcript = $variation->{transcript};
      my $downstreamOfFrameshift = $variation->{downstream_of_frameshift};
      my $position_in_codon = $variation->{position_in_codon};
      my $codon = $variation->{codon};
      my $codons =  &calculatePossibleCodons($codon);
      my $codonsLen = scalar @$codons;
      $codonsLen=$codonsLen-1;
      foreach my $i (0..$codonsLen) {
          $codon = $codons->[$i];
          my $product = $CODON_TABLE->translate($codon);
          my $pro = { "product" => $product,
                      "transcript" => $transcript,
		      "count" => $productCounts{$product},
                      "codon" => $codon,
	              "position_in_codon" => $position_in_codon,
	              "ref_location_cds" => $refLocationCds,
		      "ref_location_protein" => $refLocationProtein,
		      "downstream_of_frameshift" => $downstreamOfFrameshift	  
	             };
	  push @$products, $pro;      
      }
  }
  return $products;
}


sub makeSNPFeatureFromVariations {
  my ($variations, $referenceVariation) = @_;
  my $location = $referenceVariation->{location};
  my $sourceId = $referenceVariation->{sequence_source_id};
  my $referenceStrain = $referenceVariation->{strain};
  my %alleleCounts;
  my %productCounts;
  my %strains;
  my $totalAlleleCount = scalar @$variations;
  my $hasStopCodon = 0;
  my $transcript = $variations->[0]->{transcript};
  foreach my $variation (@$variations) {
    my $allele = $variation->{base};
    my $strain = $variation->{strain};
    $alleleCounts{$allele} ++;
    $strains{$strain}++; 
    my $products = $variation->{product};
    my @productsArray = split(":", $products);
    my $productsLen = scalar @productsArray;
    $productsLen = $productsLen-1;
    foreach my $i (0..$productsLen) {
	my $product = $productsArray[$i];
        $productCounts{$product}++;
	$hasStopCodon = 1 if($product eq '*');
    }
  }
  my $distinctStrainCount = scalar keys %strains;
  my $distinctAlleleCount =  scalar keys %alleleCounts;
  my $hasNonSynonymousAllele = scalar keys %productCounts > 1 ? 1 : 0;
  my @sortedAlleles = sort { ($alleleCounts{$b} <=> $alleleCounts{$a}) || ($a cmp $b) } keys %alleleCounts;
  my @sortedProducts = sort { ($productCounts{$b} <=> $productCounts{$a}) || ($a cmp $b) } keys %productCounts;
  my @sortedAlleleCounts = map {$alleleCounts{$_}} @sortedAlleles;
  my $majorAllele = $sortedAlleles[0];
  my $minorAllele = $sortedAlleles[1];
  my $majorProduct = $sortedProducts[0];
  my $minorProduct = $sortedProducts[1];
  my $majorAlleleCount = $sortedAlleleCounts[0];
  my $minorAlleleCount = $sortedAlleleCounts[1];
  my $snp = {     "source_id" => $sourceId,
	          "location" => $location,
	          "reference_strain" => $referenceStrain,
	          "reference_na" => $referenceVariation->{base},
	          "has_nonsynonymous_allele" => $hasNonSynonymousAllele,
		  "major_allele" => $majorAllele,
		  "minor_allele" => $minorAllele,
		  "major_allele_count" => $majorAlleleCount,
		  "minor_allele_count" => $minorAlleleCount,
		  "major_product" => $majorProduct,
		  "minor_product" => $minorProduct,
		  "distinct_strain_count" => $distinctStrainCount,
		  "distinct_allele_count" => $distinctAlleleCount,
		  "has_coding_mutation" => $referenceVariation->{is_coding},
		  "total_allele_count" => $totalAlleleCount,
		  "has_stop_codon" => $hasStopCodon,
		  "ref_codon" => $referenceVariation->{ref_codon},
		  "transcript_id" => $transcript    
  };
  return $snp;
}

sub addTranscript {
    my ($variations, $location, $transcriptSummary) = @_;
    my @sorted_keys = nsort keys %{ $transcriptSummary };
    foreach my $variation ( @$variations ) {
	my $isCoding = 0;
	my $variationTranscript;
        foreach my $transcript (@sorted_keys) {
	    my ($shifted_cds_start, $shifted_cds_end);
	    my $cds_count = 0;
	    foreach my $key (keys %{ $$transcriptSummary{$transcript} }) {
	        if ($key =~ /cds/) {
		    $cds_count++;
	        }
	    }
            $cds_count = $cds_count/2;
	    foreach my $number (1 .. $cds_count) {
		my $cdsStartField = "cds_start_$number";
		my $cdsEndField = "cds_end_$number";
		my $cds_start = $$transcriptSummary{$transcript}->{$cdsStartField};
	        my $cds_end = $$transcriptSummary{$transcript}->{$cdsEndField};
		if ($location >= $cds_start && $location <= $cds_end) {
		    $isCoding = 1;
                    $variationTranscript = $transcript;
                }
		if ($isCoding == 1) {
		    $variation->{cds_number} = $number;
		    last;
		}
            }
	    if ($variationTranscript) {
		last;
	    }   
        }
	$variation->{is_coding} = $isCoding;
	$variation->{transcript} = $variationTranscript;
    }
    return $variations;
}

sub addTranscriptomicLocation {
    my ($transcriptSummary, $currentShifts) = @_;
    my @sorted_keys = nsort keys %{ $transcriptSummary };
    my $refSeqSourceId;
    foreach my $strain (keys %{ $currentShifts }) {
        my $cdsStartLocation = 1;
        foreach my $transcript (@sorted_keys) {
	    my $cds_count = 0;
            $refSeqSourceId = $$transcriptSummary{$transcript}->{sequence_source_id};
            foreach my $key (keys %{ $$transcriptSummary{$transcript} }) {
                if ($key eq 'cds_strand') {
                    next;
                }
                elsif ($key =~ /cds/) {
                    $cds_count++;
                }
                else {
                    next;
                }
            }
            $cds_count = $cds_count/2;
            foreach my $number (1 .. $cds_count) {
                my $cdsShiftedStartField = "cds_shifted_start_$number";
                my $cdsShiftedEndField = "cds_shifted_end_$number";
		my $transcriptomicShiftedStartField = "transcript_shifted_start_$number";
                my $transcriptomicShiftedEndField = "transcript_shifted_end_$number";
                my $cds_shifted_start = $$transcriptSummary{$transcript}->{$strain}->{$cdsShiftedStartField};
                my $cds_shifted_end = $$transcriptSummary{$transcript}->{$strain}->{$cdsShiftedEndField};
                my $cdsLen = $cds_shifted_end - $cds_shifted_start;
		my $transcriptomicShiftedStart = $cdsStartLocation;
                if ($cdsStartLocation == 1) {
		    my $transcriptomicShiftedEnd = $cdsLen + 1;
		    $$transcriptSummary{$transcript}->{$strain}->{$transcriptomicShiftedStartField} = $transcriptomicShiftedStart;
		    $$transcriptSummary{$transcript}->{$strain}->{$transcriptomicShiftedEndField} = $transcriptomicShiftedEnd;
		    $cdsStartLocation = $transcriptomicShiftedEnd + 1;
		}
                else {
                    my $transcriptomicShiftedEnd = $cdsLen + 1 + $transcriptomicShiftedStart;
		    $$transcriptSummary{$transcript}->{$strain}->{$transcriptomicShiftedStartField} = $transcriptomicShiftedStart;
		    $$transcriptSummary{$transcript}->{$strain}->{$transcriptomicShiftedEndField} = $transcriptomicShiftedEnd;
		    $cdsStartLocation = $transcriptomicShiftedEnd + 1;
		}
            }
        }
    }
    my $cdsStartLocation = 1;
    foreach my $transcript (@sorted_keys) {
        my $cds_count = 0;
        foreach my $key (keys %{ $$transcriptSummary{$transcript} }) {
            if ($key eq 'cds_strand') {
                next;
            }
            elsif ($key =~ /cds/) {
                $cds_count++;
            }
            else {
                next;
            }
        }
        $cds_count = $cds_count/2;
        foreach my $number (1 .. $cds_count) {
            my $cdsStartField = "cds_start_$number";
            my $cdsEndField = "cds_end_$number";
	    my $transcriptomicStartField = "transcript_start_$number";
            my $transcriptomicEndField = "transcript_end_$number";
            my $cds_start = $$transcriptSummary{$transcript}->{$cdsStartField};
            my $cds_end = $$transcriptSummary{$transcript}->{$cdsEndField};
            my $cdsLen = $cds_end - $cds_start;
            my $transcriptomicStart = $cdsStartLocation;
	    if ($cdsStartLocation == 1) {
		my $transcriptomicEnd = $cdsLen + 1;
	        $$transcriptSummary{$transcript}->{$transcriptomicStartField} = $transcriptomicStart;
	        $$transcriptSummary{$transcript}->{$transcriptomicEndField} = $transcriptomicEnd;
	        $cdsStartLocation = $transcriptomicEnd + 1;
	    }
            else {
                my $transcriptomicEnd = $cdsLen + 1 + $transcriptomicStart;
	        $$transcriptSummary{$transcript}->{$transcriptomicStartField} = $transcriptomicStart;
	        $$transcriptSummary{$transcript}->{$transcriptomicEndField} = $transcriptomicEnd;
	        $cdsStartLocation = $transcriptomicEnd + 1;
	    }
        }
    }
    return $transcriptSummary;
}

sub reverse_complement_IUPAC {
    my ($dna) = @_;
        
    # reverse the DNA sequence
    my $revcomp = reverse($dna);

    # complement the reversed DNA sequence
    $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return $revcomp;
}


sub addFrameShifts {
    my ($currentShifts,$transcriptSummary) = @_;
    my @sorted_keys = nsort keys %{ $transcriptSummary };
    my $refSeqSourceId;
    foreach my $strain (keys %{ $currentShifts }) {
        foreach my $transcript (@sorted_keys) {
	    my $hasFrameShift = 0;
	    my $frameShiftLocation;
	    my $cds_count = 0;
            foreach my $key (keys %{ $$transcriptSummary{$transcript} }) {
                if ($key eq 'cds_strand') {
                    next;
                }
                elsif ($key =~ /cds/) {
                    $cds_count++;
                }
                else {
                    next;
                }
            }
            $cds_count = $cds_count/2;
            foreach my $number (1 .. $cds_count) {
                my $cdsStartField = "cds_start_$number";
                my $cdsEndField = "cds_end_$number";
		my $cds_start = $$transcriptSummary{$transcript}->{$cdsStartField};
                my $cds_end = $$transcriptSummary{$transcript}->{$cdsEndField};
		my $shiftArray = $currentShifts->{$strain};
	        foreach my $chromosome (keys %{ $shiftArray }) {	   
	            my @chromosomeShiftArray = $shiftArray->{$chromosome};
	            my $indexedArray = $chromosomeShiftArray[0][0];
		    my $shiftArrayLen = scalar @{ $indexedArray };
                    my $shiftFrameLimit = $shiftArrayLen - 1;
		    my $shiftFrame = 0;
		    until ($shiftFrame == $shiftFrameLimit) {
			my $location = $indexedArray->[$shiftFrame][0];
			if ($location >= $cds_start && $location <= $cds_end) {
			    #print "Location $location start $cds_start end $cds_end\n";
			    $hasFrameShift = 1;
			    $frameShiftLocation = $location;
			    last;
			}
			else {
			    $shiftFrame+=1;
			}
		    }
		}
	    }
	    $$transcriptSummary{$transcript}->{$strain}->{has_frameshift} = $hasFrameShift;
	    $$transcriptSummary{$transcript}->{$strain}->{frameshift_location} = $frameShiftLocation;
	}
    }
    return $transcriptSummary;
}


sub getVariations {
  my ($location, $snpFileLine, $cacheFileLine, $snpFile, $snpFileCount, $cacheFile, $cacheFileCount, $undoneStrains) = @_;
  my (@snpsArray, @cacheStrains, @variations, $nextCacheLocation, $nextSnpLocation, $nextLocation);
  my $currentLocation = $location;
  if ($cacheFileCount > 0) {
    while ($currentLocation == $location) {
      my $cacheLine = `sed "${cacheFileLine}q;d" $cacheFile`;
      chomp $cacheLine;
      my @cacheLineArray = split("\t", $cacheLine);
      $currentLocation = $cacheLineArray[1];
      if ($currentLocation == $location) {
        if ($cacheLineArray[2] ~~ $undoneStrains) {
          $cacheFileLine += 1;
        }
        else {
          push @snpsArray, \@cacheLineArray;
          $cacheFileLine += 1;
        }
      }
      else {
	  if ($location < $currentLocation) {
	      $nextCacheLocation = $currentLocation;
	  }
	  else {
	      $cacheFileLine += 1;
	  }
      }
    }
  }
  $currentLocation = $location;
  @cacheStrains = &getCacheStrains(\@snpsArray);
  if ($snpFileCount > 0) {
    while ($currentLocation == $location) {
      my $snpLine = `sed "${snpFileLine}q;d" $snpFile`;
      chomp $snpLine;
      my @snpLineArray = split("\t", $snpLine);
      $currentLocation = $snpLineArray[1];
      if ($currentLocation == $location) {
        if ($snpLineArray[2] ~~ @cacheStrains) {
          $snpFileLine += 1;
        }
        else {
          push @snpsArray, \@snpLineArray;
          $snpFileLine += 1;
        }    
      }
      else {
        if ($location < $currentLocation) {
          $nextSnpLocation = $currentLocation;
        }
        else {
          $snpFileLine += 1;
        }
      }
    }
  }
  if ($nextCacheLocation) {
    if ($nextSnpLocation) {  
      $nextSnpLocation = $nextCacheLocation <= $nextSnpLocation ? $nextCacheLocation : $nextSnpLocation;
    }
    else {
      $nextSnpLocation = $nextCacheLocation;
    }
  }
  else {
    $nextSnpLocation = $nextSnpLocation;
  }
  @variations = &makeVariationsFromArray(@snpsArray);
  return(\@variations, $snpFileLine, $cacheFileLine, $nextSnpLocation);
}

sub getCacheStrains {
  my ($snpsArray) = @_;
  my @cacheStrains;
  foreach my $snp (@$snpsArray) {
      push @cacheStrains, $snp->[2];
  }
  return \@cacheStrains;
}

sub makeVariationsFromArray {
  my (@snps) = @_;
  my @variations;
  foreach my $snp (@snps) {
    my $variation;
    $variation->{sequence_source_id} = $snp->[0];
    $variation->{location} = $snp->[1];
    $variation->{strain} = $snp->[2];
    $variation->{reference} = $snp->[3];
    $variation->{base} = $snp->[4];
    $variation->{coverage} = $snp->[5];
    $variation->{percent} = $snp->[6];
    $variation->{quality} = $snp->[7];
    $variation->{pvalue} = $snp->[8];
    $variation->{snp_source_id} = $snp->[9];
    $variation->{is_coding} = $snp->[10];
    $variation->{position_in_cds} = $snp->[11];
    $variation->{position_in_codon} = $snp->[12];
    $variation->{downstream_of_frameshift} = $snp->[13];
    $variation->{transcript} = $snp->[14];
    my $productString = $snp->[15];
    my @productsArray = split(":", $productString);
    $variation->{product} = ();
    $variation->{product} = \@productsArray;
    $variation->{reference_codon} = $snp->[16];
    $variation->{codon} = $snp->[17];
    $variation->{has_nonsynonomous} = $snp->[18];
    $variation->{shifted_location} = $snp->[19];
    $variation->{cds_number} = $snp->[20];
    $variation->{current_shift} = $snp->[21];
    push @variations, $variation;
  }
  return @variations;
}


sub refProduct {
    my ($transcriptSummary, $location, $transcript, $genomeFasta, $organismAbbrev, $cds_number) = @_;
    my ($refProduct, $refCodon, $ref_pos_in_cds, $ref_pos_in_pro, $referenceCodingSequence);
    my $cds_count = 0;
    my $prior_cds_len = 0;
    my $prior_ref_cds_len = 0;
    my ($pos_in_cds, $pos_in_pro, $product, $codon);
    foreach my $key (keys %{ $$transcriptSummary{$transcript} }) {
        if ($key eq 'cds_strand') {
	    next;
	}
        elsif ($key =~ /cds/) {
           $cds_count++;
	}
	else {
	    next;
	}
    }
    $cds_count = $cds_count/2;
    if($transcriptSummary->{$transcript}->{cache}->{ref_cds}) { # We already have retrieved the reference coding sequence for this transcript
        $referenceCodingSequence = $transcriptSummary->{$transcript}->{cache}->{ref_cds};
	foreach my $number (1 .. $cds_count) {
	    my $cdsStartField = "cds_start_$number";
	    my $cdsEndField = "cds_end_$number";
	    my $cds_start = $$transcriptSummary{$transcript}->{$cdsStartField};
	    my $cds_end = $$transcriptSummary{$transcript}->{$cdsEndField};
	    if ($number != $cds_number) {
	        $prior_ref_cds_len = $prior_ref_cds_len + $cds_end - $cds_start;
	    }
	    elsif ($number == $cds_number) {
	        $ref_pos_in_cds = $location - $cds_start + $prior_ref_cds_len + 1;
	    }
	}
    }
    else { # first time through for this transcript. Use same functionality for retrieving the reference coding sequence
        foreach my $number (1 .. $cds_count) {
            my $cdsStartField = "cds_start_$number";
	    my $cdsEndField = "cds_end_$number";
	    my $cds_start = $$transcriptSummary{$transcript}->{$cdsStartField};
	    my $cds_end = $$transcriptSummary{$transcript}->{$cdsEndField};
	    my $cds_ref_sequence_chunk = &getCodingSequence($organismAbbrev, $cds_start, $cds_end, $strand, $genomeFasta);
	    $referenceCodingSequence = $referenceCodingSequence . $cds_ref_sequence_chunk;
	    if ($number != $cds_number) {
	        $prior_ref_cds_len = $prior_ref_cds_len + $cds_end - $cds_start;
	    }
	    elsif ($number == $cds_number) {
	        $ref_pos_in_cds = $location - $cds_start + $prior_ref_cds_len + 1;
	    }
        }
	$transcriptSummary->{$transcript}->{cache}->{ref_cds} = $referenceCodingSequence;
    }

    ($refCodon, $refProduct) = &getAminoAcidSequenceOfSnp($referenceCodingSequence, $ref_pos_in_cds);
    $ref_pos_in_pro = &calculateAminoAcidPosition($ref_pos_in_cds);
    $refProduct = $refProduct->[0];
    return($refProduct, $refCodon, $ref_pos_in_cds, $ref_pos_in_pro, $transcriptSummary);
}


sub variationProduct {
    my ($transcriptSummary, $variations, $consensusFasta, $refProduct, $refCodon, $ref_pos_in_cds, $ref_pos_in_pro) = @_;
    foreach my $variation (@$variations) {
	if ($variation->{position_in_codon}) {
	    next;
	}
	my $consensusCodingSequence;
	my $strain = $variation->{strain};
	my $strand = $variation->{strand};
	my $transcript = $variation->{transcript};
	my $cds_number = $variation->{cds_number};
	my $shifted_location = $variation->{shifted_location};
	my $location = $variation->{location};
	my $cds_count = 0;
        my $prior_cds_len = 0;
	my $prior_ref_cds_len = 0;
	my ($pos_in_cds, $pos_in_pro, $product, $codon);
        foreach my $key (keys %{ $$transcriptSummary{$transcript} }) {
	    if ($key eq 'cds_strand') {
		next;
	    }
            elsif ($key =~ /cds/) {
                $cds_count++;
	    }
	    else {
	        next;
	    }
	}
	$cds_count = $cds_count/2;
        if($transcriptSummary->{$transcript}->{cache}->{$strain}->{consensus_cds}) {
            $consensusCodingSequence = $transcriptSummary->{$transcript}->{cache}->{$strain}->{consensus_cds};
	    foreach my $number (1 .. $cds_count) {
                my $cdsShiftedStartField = "cds_shifted_start_$number";
	        my $cdsShiftedEndField = "cds_shifted_end_$number";
                my $cds_shifted_start = $$transcriptSummary{$transcript}->{$strain}->{$cdsShiftedStartField};
                my $cds_shifted_end = $$transcriptSummary{$transcript}->{$strain}->{$cdsShiftedEndField};
	        if ($number != $cds_number) {
	            $prior_cds_len = $prior_cds_len + $cds_shifted_end - $cds_shifted_start;
	        }
	        elsif ($number == $cds_number) {
	            $pos_in_cds = $shifted_location - $cds_shifted_start + $prior_cds_len + 1;
	        }
            }
        }
        else { # first time through for this transcript. Get coding sequence using samtools faidx and consensus sequence using shifted cds
            foreach my $number (1 .. $cds_count) {
                my $cdsShiftedStartField = "cds_shifted_start_$number";
	        my $cdsShiftedEndField = "cds_shifted_end_$number";
	        my $cds_shifted_start = $$transcriptSummary{$transcript}->{$strain}->{$cdsShiftedStartField};
	        my $cds_shifted_end = $$transcriptSummary{$transcript}->{$strain}->{$cdsShiftedEndField};
	        my $cds_sequence_chunk = &getCodingSequence($strain, $cds_shifted_start, $cds_shifted_end, $strand, $consensusFasta);
	        $consensusCodingSequence = $consensusCodingSequence . $cds_sequence_chunk;
	        if ($number != $cds_number) {
	            $prior_cds_len = $prior_cds_len + $cds_shifted_end - $cds_shifted_start;
	        }
	        elsif ($number == $cds_number) {
	            $pos_in_cds = $shifted_location - $cds_shifted_start + $prior_cds_len + 1;
	        }
            }
	    $transcriptSummary->{$transcript}->{cache}->{$strain}->{consensus_cds} = $consensusCodingSequence;
	}

	($codon, $product) = &getAminoAcidSequenceOfSnp($consensusCodingSequence, $pos_in_cds);
        $pos_in_pro = &calculateAminoAcidPosition($pos_in_cds);
	
	$refProduct = $refProduct->[0];
        $variation->{product} = $product;
	$variation->{position_in_codon} = $pos_in_pro;
	$variation->{position_in_cds} = $pos_in_cds;
        $variation->{codon} = $codon;
	$variation->{reference_codon} = $refCodon;
    }
    return($transcriptSummary, $variations);
}

sub checkForDownstreamOfFrameshift {
    my ($transcriptSummary, $variation) = @_;
    my $downstreamOfFrameshift;
    my $transcript = $variation->{transcript};
    my $strain = $variation->{strain};
    if ($transcriptSummary->{$transcript}->{$strain}->{has_frameshift} == 1) {
        my $frameshiftLocation = $transcriptSummary->{$transcript}->{$strain}->{frameshift_location};
        if ($transcriptSummary->{$transcript}->{cds_strain} == 1 && $variation->{shifted_location} > $frameshiftLocation) {
            $downstreamOfFrameshift = 1;
        }
        elsif ($transcriptSummary->{$transcript}->{cds_strain} == -1 && $variation->{shifted_location} < $frameshiftLocation) {
            $downstreamOfFrameshift = 1;
        }
        else {
            $downstreamOfFrameshift = 0;
        }
    }
    else {
        $downstreamOfFrameshift = 0;
    }
    return $downstreamOfFrameshift;
}
