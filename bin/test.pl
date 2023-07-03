#!/usr/bin/perl

use File::Basename;
use Data::Dumper;
use Getopt::Long;
use GUS::ObjRelP::DbiDatabase;
use GUS::Supported::GusConfig;
use CBIL::Bio::SequenceUtils;
use Bio::Seq;
use Bio::Tools::GFF;
use Bio::Coordinate::GeneMapper;
use Bio::Coordinate::Pair;
use Bio::Location::Simple;
use Bio::Tools::CodonTable;
use GUS::Community::GeneModelLocations;
use VEuPath::SnpUtils  qw(sampleCacheFileColumnNames snpFileColumnNames alleleFileColumnNames productFileColumnNames);
use DBI;
use DBD::Oracle;
use VEuPath::MergeSortedSeqVariations;
use ApiCommonData::Load::FileReader;
use locale;
use Sort::Naturally;
use Set::CrossProduct;

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
	    "gtfFile=s"=> \$gtfFile,
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

my $initialCacheCount = `cat $cacheFile | wc -l`;
chomp($initialCacheCount);

unless(-d $varscanDirectory) {
  &usage("Required Directory Missing") unless($isLegacyVariations);
}

unless(-e $undoneStrainsFile) {
  open(FILE, "> $undoneStrainsFile") or die "Could not open file $undoneStrainsFile for writing: $!";
  close FILE;
}

my $CODON_TABLE = Bio::Tools::CodonTable->new( -id => 1); #standard codon table

my $totalTime;
my $totalTimeStart = time();

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

open(UNDONE, $undoneStrainsFile) or die "Cannot open file $undoneStrainsFile for reading: $!";
my @undoneStrains =  map { chomp; $_ } <UNDONE>;
close UNDONE;

if($forcePositionCompute) {
  push @undoneStrains, $referenceStrain;
}

my $merger = VEuPath::MergeSortedSeqVariations->new($newSampleFile, $cacheFile, \@undoneStrains, qr/\t/);

my ($prevSequenceId, $prevTranscriptMaxEnd, $prevTranscript);

my $strainFrame;
my $count = 0;

my ($prevSequenceId, $prevTranscriptMaxEnd, $prevTranscripts, $counter);

while($merger->hasNext()) {

  my $variations = $merger->nextSNP();
  my ($sequenceId, $location) = &snpLocationFromVariations($variations);

  print STDERR "SEQUENCEID=$sequenceId\tLOCATION=$location\n" if($debug);

  my ($referenceAllele, $positionsInCds, $positionsInProtein, $referenceVariation, $isCoding);
  
  ($variations, $strainFrame) = &addShiftedLocation($variations, $strainFrame, $currentShifts);

  $variations = &addTranscript($variations, $location, $transcriptSummary);

  my $cachedReferenceVariation = &cachedReferenceVariation($variations, $referenceStrain);

  if($cachedReferenceVariation && !$isLegacyVariations) {
      print STDERR "HAS_CACHED REFERENCE VARIATION\n" if($debug);
      $referenceVariation = $cachedReferenceVariation;
  }
  else {
      $referenceAllele = $variations->[0]->{reference};
      $isCoding = $variations->[0]->{is_coding};
      my ($refProduct, $refCodon, $refPositionInCds, $refPositionInProtein, $adjacentSnpCausesProductDifference, $reference_aa_full);
      if ($isCoding == 1) {
	  ($refProduct, $refCodon, $refPositionInCds, $refPositionInProtein, $adjacentSnpCausesProductDifference, $reference_aa_full) = &variationAndRefProduct($transcriptSummary, $variations, $consensusFasta, $genomeFasta);
      }
      $referenceVariation = {'base' => $referenceAllele,
      			     'reference' => $referenceAllele,    
                             'location' => $location,
                             'sequence_source_id' => $sequenceId,
                             'matches_reference' => 1,
                             'position_in_cds' => $refPositionInCds,
                             'strain' => $referenceStrain,
                             'product' => $refProduct,
                             'position_in_protein' => $refPositionInProtein,
                             'is_coding' => $isCoding,
                             'has_nonsynonomous' => $adjacentSnpCausesProductDifference,
      			     'ref_codon' => $refCodon,
      			     'reference_aa_full' => $reference_aa_full	 
      };
      push @$variations, $referenceVariation;
      if ($isCoding == 1) {
	  #print Dumper $variations;
      }
  }
  #print Dumper $variations;
}
#die;

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

sub cachedReferenceVariation {
    my ($variations, $referenceStrain) = @_;

    foreach(@$variations) {
      return $_ if($_->{strain} eq $referenceStrain);
    }
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
  my ($transcriptSummary, $transcripts) = @_;

  foreach my $transcript (@$transcripts) {
    $transcriptSummary->{$transcript}->{cache} = undef;
  }
}

sub lookupByLocation {
  my ($sequenceId, $l, $geneLocs) = @_;

  return(undef) unless(ref ($geneLocs->{$sequenceId}) eq 'ARRAY');

  my @locations = @{$geneLocs->{$sequenceId}};

  my $startCursor = 0;
  my $endCursor = scalar(@locations) - 1;
  my $midpoint;

  return(undef) if($l < $locations[$startCursor]->{start} || $l > $locations[$endCursor]->{end});

  while ($startCursor <= $endCursor) {
    $midpoint = int(($endCursor + $startCursor) / 2);

    my $location = $locations[$midpoint];

    if ($l > $location->{start}) {
      $startCursor = $midpoint + 1;
    } 
    elsif ($l < $location->{start}) {
      $endCursor = $midpoint - 1;
    }
    else {  }

    if($l >= $location->{start} && $l <= $location->{end}) {
      return($location);
    }
  }
  return(undef);
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
    return \%transcriptSummary;
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

  print "$cdsSequence\t$positionInCds\t$positionInProtein\t$codon\n";
  
  my $codons = &calculatePossibleCodons($codon);

  my $products;
  my $productsLen = scalar @$codons;
  $productsLen=$productsLen-1;

  foreach my $i (0..$productsLen) {
      $current_codon = $codons->[$i];
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
	    $oldCdsShift = 0;
	    $cdsShiftFrame = 0;
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
    my ($location, $strain, $indexedArray, $shiftArrayLen, $shiftFrameLimit, $oldShift, $shiftFrame, $shiftedLocation);

    foreach my $variation (@$variations) {

	$location = $variation->{location};
        $strain = $variation->{strain};
	$chromosome = $variation->{sequence_source_id};
        my @shiftArray = $currentShifts->{$strain}->{$chromosome};
        $indexedArray = $shiftArray[0][0];
        $shiftArrayLen = scalar @{ $indexedArray };
        $shiftFrameLimit = $shiftArrayLen - 1;
        $oldShift;
        $shiftFrame;

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
	    $currentShift = 0;
	    $oldShift = 0;
	}

	$variation->{shifted_location} = $shiftedLocation;
        $variation->{current_shift} = $oldShift;
        $strainFrame->{$strain}->{$chromosome}->{oldShift} = $oldShift;
        $strainFrame->{$strain}->{$chromosome}->{shiftFrame} = $shiftFrame;
    }
    return ($variations, $strainFrame);
}




sub variationAndRefProduct {
    my ($transcriptSummary, $variations, $consensusFasta, $genomeFasta) = @_;
    my ($product, $refProduct, $codon, $refCodon);
    my $adjacentSnpCausesProductDifference = 0;
    my $refConsensusCodingSequence;
    foreach my $variation (@$variations) {
	my $refSeqSourceId = $variation->{sequence_source_id};
	my $consensusCodingSequence;
	my $strain = $variation->{strain};
        my $strand = $transcriptSummary->{$transcript}->{cds_strand};
	my $strand = $variation->{strand};
	my $transcript = $variation->{transcript};
	my $cds_number = $variation->{cds_number};
	my $shifted_location = $variation->{shifted_location};
	my $location = $variation->{location};
	my $cds_count = 0;
        my $prior_cds_len = 0;
	my $prior_ref_cds_len = 0;
	my ($pos_in_cds, $pos_in_protein, $ref_pos_in_cds, $ref_pos_in_protein);
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
	            $prior_cds_len = $prior_cds_len + $cds_end - $cds_start;
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
	            $prior_cds_len = $prior_cds_len + $cds_end - $cds_start;
	        }
	        elsif ($number == $cds_number) {
	            $pos_in_cds = $shifted_location - $cds_shifted_start + $prior_cds_len + 1;
	        }
            }
	    $transcriptSummary->{$transcript}->{cache}->{$strain}->{consensus_cds} = $consensusCodingSequence;
	}
	if($transcriptSummary->{$transcript}->{cache}->{ref_cds}) { # We already have retrieved the reference coding sequence for this transcript
	    $refConsensusCodingSequence = $transcriptSummary->{$transcript}->{cache}->{ref_cds};
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
		my $cds_ref_sequence_chunk = &getCodingSequence($refSeqSourceId, $cds_start, $cds_end, $strand, $genomeFasta);
		$refConsensusCodingSequence = $refConsensusCodingSequence . $cds_ref_sequence_chunk;
		if ($number != $cds_number) {
		    $prior_ref_cds_len = $prior_ref_cds_len + $cds_end - $cds_start;
		}
		elsif ($number == $cds_number) {
		    $ref_pos_in_cds = $location - $cds_start + $prior_ref_cds_len + 1;
		}
            }
	    $transcriptSummary->{$transcript}->{cache}->{ref_cds} = $refConsensusCodingSequence;
        }

	print "Reference\n";
	($refCodon, $refProduct) = &getAminoAcidSequenceOfSnp($refConsensusCodingSequence, $ref_pos_in_cds);
	$ref_pos_in_pro = &calculateAminoAcidPosition($ref_pos_in_cds);
	print "Variation\n";
        ($codon, $product) = &getAminoAcidSequenceOfSnp($consensusCodingSequence, $pos_in_cds);
        $pos_in_pro = &calculateAminoAcidPosition($pos_in_cds);
        if($product ne $refProduct) {
            $adjacentSnpCausesProductDifference = 1;
        }
	$refProduct = $refProduct->[0];
        $variation->{product} = $product;
        $variation->{position_in_protein} = $pos_in_pro;
	$variation->{position_in_cds} = $pos_in_cds;
        $variation->{codon} = $codon;
	$variation->{reference_codon} = $refCodon;
	$variation->{has_nonsynonomous} = $adjacentSnpCausesProductDifference;
	$variation->{reference_aa} = $refProduct;
    }
    return($refProduct, $refCodon, $ref_pos_in_cds, $ref_pos_in_pro, $adjacentSnpCausesProductDifference, $refConsensusCodingSequence);
}


sub getCodingSequence {
    my ($defline, $start, $end, $strand, $fasta) = @_;
    my $seq;

    if ($strand == 1) {
	$seq = `samtools faidx $fasta $defline:$start-$end`;
    }

    else {
        $seq = `samtools faidx -i $fasta $defline:$start-$end`;
    }

    my $seq = $seq =~ s/>.+\n//gr;
    my $seq = $seq =~ s/\n//gr;
    return $seq;
}


sub calculatePossibleCodons {
    my ($codon) = @_;
    my $codonList;

    if (!$codon) {
        return $codonList;
    }
    else {
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
}


sub printSNPFeature {
    my ($snp, $snpFh) = @_;
    my $keys = VEuPath::SnpUtils::snpFileColumnNames();
    print $snpFh join("\t", map {$snp->{$_}} @$keys) . "\n";
}

sub printProductFeature {
    my ($products, $productFh)= @_;
    $keys = VEuPath::SnpUtils::productFileColumnNames();
    foreach my $product ($products) {
	print $productFh join("\t", map {$product->[0]->{$_}} @$keys) . "\n";
    }
}

sub printAlleleFeature {
    my ($alleles, $alleleFh)= @_;
    $keys = VEuPath::SnpUtils::alleleFileColumnNames();
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
      my $distinct_strain_count;
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
  my $refLocationProtein = $referenceVariation->{position_in_protein};
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
		      "ref_location_protein" => $refLocationProtein
	             };
	  push @$products, $pro;      
      }
  }
  return $products;
}


sub makeSNPFeatureFromVariations {
  my ($variations, $referenceVariation) = @_;
  my $location = $referenceVariation->{location};
  my $snpSourceId = $referenceVariation->{snp_source_id} ? $referenceVariation->{snp_source_id} : "NGS_SNP.$referenceVariation->{sequence_source_id}.$location";
  my $referenceStrain = $referenceVariation->{strain};
  my %alleleCounts;
  my %productCounts;
  my %strains;
  my $totalAlleleCount = scalar @$variations;
  my $hasStopCodon = 0;
  foreach my $variation (@$variations) {
    my $geneNaSequenceId = $variation->{gene_na_sequence_id};
    my $sequenceSourceId = $variation->{sequence_source_id};  
    my $allele = $variation->{base};
    my $strain = $variation->{strain};
    $alleleCounts{$allele} ++;
    $strains{$strain}++; 
    my $products = $variation->{product};
    my $productsLen = scalar @$products;
    $productsLen = $productsLen-1;
    foreach my $i (0..$productsLen) {
	my $product = $products->[$i];
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
  my $snp = {     "source_id" => $snpSourceId,
	          "gene_na_sequence_id" => $geneNaSequenceId,
	          "location" => $location,
	          "reference_strain" => $referenceStrain,
	          "reference_na" => $referenceVariation->{base},
	          "reference_aa" => $referenceVariation->{product}->[0],
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
		  "reference_aa_full" => $referenceVariation->{reference_aa_full}    
            };
  return $snp;
}

sub addTranscript {
    my ($variations, $location, $transcriptSummary) = @_;
    my @sorted_keys = nsort keys %{ $transcriptSummary };
    foreach $variation ( @$variations ) {
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
