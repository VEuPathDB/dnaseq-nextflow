#!/usr/bin/perl

use File::Basename;
use Data::Dumper;
use Getopt::Long;
use VEuPath::DbiDatabase;
use VEuPath::GusConfig;
use VEuPath::SequenceUtils;
use Bio::Seq;
use Bio::Tools::GFF;
use Bio::Coordinate::GeneMapper;
use Bio::Coordinate::Pair;
use Bio::Location::Simple;
use Bio::Tools::CodonTable;
use VEuPath::GeneModelLocations;
use VEuPath::SnpUtils  qw(fileColumnNames);
use DBI;
use DBD::Oracle;
use VEuPath::MergeSortedSeqVariations;
use VEuPath::FileReader;
use locale;

my ($newSampleFile, $cacheFile, $cleanCache, $transcriptExtDbRlsSpec, $organismAbbrev, $undoneStrainsFile, $gusConfigFile, $varscanDirectory, $referenceStrain, $help, $debug, $extDbRlsSpec, $isLegacyVariations, $forcePositionCompute);
&GetOptions("new_sample_file=s"=> \$newSampleFile,
            "cache_file=s"=> \$cacheFile,
            "clean_cache"=> \$cleanCache,
            "gusConfigFile|gc=s"=> \$gusConfigFile,
            "undone_strains_file=s" => \$undoneStrainsFile,
            "varscan_directory=s" => \$varscanDirectory,
            "is_legacy_variations" => \$isLegacyVariations,
            "transcript_extdb_spec=s" => \$transcriptExtDbRlsSpec,
            "extdb_spec=s" => \$extDbRlsSpec,
            "organism_abbrev=s" =>\$organismAbbrev,
            "reference_strain=s" => \$referenceStrain,
            "force_position_compute" => \$forcePositionCompute,
            "debug" => \$debug,
            "help|h" => \$help,
    );
if($help) {
  &usage();
}
$gusConfigFile = $ENV{GUS_HOME} . "/config/gus.config" unless(-e $gusConfigFile);
my $cacheFileExists = -e $cacheFile;
if(!$cacheFileExists || $cleanCache) {
  print STDERR "Creating empty cache because none exists.\n"
    if !$cacheFileExists;
  print STDERR "Creating empty cache because \"--clean_cache\" option was specified.\n"
    if $cleanCache;
  open(CACHE, ">$cacheFile") or die "Cannot create a cache file: $!";
  close CACHE;
}
my $initialCacheCount = `cat $cacheFile | wc -l`;
chomp($initialCacheCount);
print STDERR "Starting with cache file of $initialCacheCount records\n";

unless(-e $newSampleFile && -e $gusConfigFile) {
  &usage("Required File Missing");
}
unless(-d $varscanDirectory) {
  &usage("Required Directory Missing") unless($isLegacyVariations);
}
unless($transcriptExtDbRlsSpec && $organismAbbrev && $referenceStrain && $extDbRlsSpec) {
  &usage("Missing Required param value");
}
unless(-e $undoneStrainsFile) {
  open(FILE, "> $undoneStrainsFile") or die "Could not open file $undoneStrainsFile for writing: $!";
  close FILE;
}
my $CODON_TABLE = Bio::Tools::CodonTable->new( -id => 1); #standard codon table
my $totalTime;
my $totalTimeStart = time();
my $gusConfig = VEuPath::GusConfig->new($gusConfigFile);

my ($dbiDsn, $login, $password, $core);
$dbiDsn = $gusConfig->getDbiDsn();
$login = $gusConfig->getDatabaseLogin();
$password = $gusConfig->getDatabasePassword();
$core = $gusConfig->getCoreSchemaName();

my $dbh = DBI->connect('dbi:Oracle:database=tryp-inc;SERVICE_NAME=trypbl8n.upenn.edu;host=localhost;port=1528', $login, $password);

my $SEQUENCE_QUERY = "select substr(s.sequence, ?, ?) as base
                      from dots.nasequence s
                     where s.source_id = ?";
my $SEQUENCE_QUERY_SH = $dbh->prepare($SEQUENCE_QUERY);
my $dirname = dirname($cacheFile);
my $tempCacheFile = $dirname . "/cache.tmp";
my $snpOutputFile = $dirname . "/snpFeature.dat";
my ($snpFh, $cacheFh);
open($cacheFh, "> $tempCacheFile") or die "Cannot open file $tempCacheFile for writing: $!";
open($snpFh, "> $snpOutputFile") or die "Cannot open file $snpOutputFile for writing: $!";

my $strainVarscanFileHandles = &openVarscanFiles($varscanDirectory, $isLegacyVariations);

my @allStrains = keys %{$strainVarscanFileHandles};

my $currentShifts = &createCurrentShifts(\@allStrains, $dbh);
#print keys %{ $currentShifts };
#die; 

my $strainExtDbRlsAndProtocolAppNodeIds = &queryExtDbRlsAndProtocolAppNodeIdsForStrains(\@allStrains, $dbh, $organismAbbrev);

my $transcriptExtDbRlsId = &queryExtDbRlsIdFromSpec($dbh, $transcriptExtDbRlsSpec);
my $thisExtDbRlsId = &queryExtDbRlsIdFromSpec($dbh, $extDbRlsSpec);


my $geneModelLocations = VEuPath::GeneModelLocations->new($dbh, $transcriptExtDbRlsId, 1);
my $agpMap = $geneModelLocations->getAgpMap();

# NOTE:  The key in this hash is actually the aa_feature_id in dots.translatedaafeature because we are interested in cds coords
my $transcriptSummary = &getTranscriptLocations($dbh, $transcriptExtDbRlsId, $agpMap);

my $geneLocations = &getGeneLocations($transcriptSummary);

my @shiftArray = $currentShifts->{'LV39cl5_chr1'};

my $indexedArray = $shiftArray[0][0];
my $shiftArrayLen = scalar @{ $indexedArray };
#print $shiftArrayLen;
#print "$shiftArrayLen\n";
#die;

$transcriptSummary = &addStrainExonShiftsToTranscriptSummary($currentShifts, $transcriptSummary);
#print Dumper $transcriptSummary;
#die;

open(UNDONE, $undoneStrainsFile) or die "Cannot open file $undoneStrainsFile for reading: $!";
my @undoneStrains =  map { chomp; $_ } <UNDONE>;
close UNDONE;

if($forcePositionCompute) {
  push @undoneStrains, $referenceStrain;
}

my $naSequenceIds = &queryNaSequenceIds($dbh);

my $merger = VEuPath::MergeSortedSeqVariations->new($newSampleFile, $cacheFile, \@undoneStrains, qr/\t/);

while($merger->hasNext()) {
  my $variations = $merger->nextSNP();
  print Dumper $variations;
}
die;

#--------------------------------------------------------------------------------
# BEGIN SUBROUTINES
#--------------------------------------------------------------------------------


sub uniqueValueFromHashRef {
    my ($hashRef) = @_;
    my %seen = ();
    my @unique = grep { ! $seen{$_} ++ } values %$hashRef;

    if(scalar @unique == 1) {
	return $unique[0];
    }
    return undef;
}

sub makeHashRefFromString {
    my ($string) = @_;

    return unless($string);

    my $rv = eval $string;

    if($@) {
	die "Error evaluating string: $@";
    }

    unless(ref($rv) eq 'HASH') {
	die "Error evaluating string [$string]to hash";
    }
    return $rv;
}


sub queryNaSequenceIds {
    my ($dbh) = @_;

  my $sql = "select s.na_sequence_id, s.source_id
from dots.nasequence s, sres.ontologyterm o
where s.sequence_ontology_id = o.ontology_term_id
and o.name in ('random_sequence', 'chromosome', 'contig', 'supercontig','mitochondrial_chromosome','plastid_sequence','cloned_genomic','apicoplast_chromosome', 'variant_genome','maxicircle')
";
    print "Preparing\n";
    my $sh = $dbh->prepare($sql);
    print "Executing\n";
    $sh->execute();

    my %naSequences;
    my $counter = 0;
    while(my ($naSequenceId, $sourceId) = $sh->fetchrow_array()) {
	$naSequences{$sourceId} = $naSequenceId;
	print "$counter\n";
	$counter++;
    }
    $sh->finish();

    return \%naSequences;
}


sub calculateCdsPosition {
    my ($transcripts, $transcriptSummary, $sequenceId, $location) = @_;

    my %cdsPositions;
    my %proteinPositions;

    my $isCoding;

    foreach my $transcript (@$transcripts) {
	my $cdsStart = $transcriptSummary->{$transcript}->{min_cds_start};
	my $cdsEnd = $transcriptSummary->{$transcript}->{max_cds_end};
	my $cdsStrand = $transcriptSummary->{$transcript}->{cds_strand};

	next unless($cdsStart && $cdsEnd);

	next if($location < $cdsStart || $location > $cdsEnd);

    my $gene = Bio::Coordinate::GeneMapper->new(
      -in  => "chr",
      -out => "cds",
      -cds => Bio::Location::Simple->new(
         -start  => $cdsStart,
         -end  => $cdsEnd,
         -strand => $cdsStrand,
         -seq_id => $sequenceId,
      ),
	-exons => $transcriptSummary->{$transcript}->{exons}
	);

    my $loc =   Bio::Location::Simple->new(
      -start => $location,
      -end   => $location,,
      -strand => +1,
      -seq_id => $sequenceId,
	);

	my $map = $gene->map($loc);

	my $cdsPos = $map->start;

	if($cdsPos && $cdsPos > 1) {
	    my $positionsInProtein = &calculateAminoAcidPosition($cdsPos);

	    $cdsPositions{$transcript} = $cdsPos;
	    $proteinPositions{$transcript} = $positionsInProtein;

	    $isCoding = 1;
	}
    }

    return($isCoding, \%cdsPositions, \%proteinPositions);
}



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

    my $keys = &variationFileColumnNames();

    print $fh join("\t", map {$variation->{$_}} @$keys) . "\n";

}

sub printSNP {
    my ($snp, $fh) = @_;

    my $keys = &snpFileColumnNames();

    print $fh join("\t", map {$snp->{$_}} @$keys) . "\n";
}


sub makeSNPFeatureFromVariations {
    my ($variations, $referenceVariation, $geneNaFeatureId, $extDbRlsId) = @_;

    my $sequenceSourceId = $referenceVariation->{sequence_source_id};
    my $location = $referenceVariation->{location};
    my $snpSourceId = $referenceVariation->{snp_source_id} ? $referenceVariation->{snp_source_id} : "NGS_SNP.$sequenceSourceId.$location";

    my $referenceStrain = $referenceVariation->{strain};

    my %alleleCounts;
    my %productCounts;
    my %strains;

    my $totalAlleleCount = scalar @$variations;
    my $hasStopCodon = 0;

    foreach my $variation (@$variations) {
	my $allele = $variation->{base};
	my $product = $variation->{product};
	my $strain = $variation->{strain};

	$hasStopCodon = 1 if($product eq '*');

	$alleleCounts{$allele} ++;
	$strains{$strain}++;

	next if(uc($product) eq 'X');
	$productCounts{$product}++;
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

    return { "gene_na_feature_id" => $geneNaFeatureId,
           "source_id" => $snpSourceId,
	     "na_sequence_id" => $referenceVariation->{na_sequence_id},
           "location" => $location,
           "reference_strain" => $referenceStrain,
	     "reference_na" => $referenceVariation->{base},
	     "reference_aa" => $referenceVariation->{product},
	     "position_in_cds" => $referenceVariation->{position_in_cds},
	     "position_in_protein" => $referenceVariation->{position_in_protein},
           "external_database_release_id" => $extDbRlsId,
           "has_nonsynonymous_allele" => $hasNonSynonymousAllele,
           "major_allele" => $majorAllele,
           "minor_allele" => $minorAllele,
           "major_allele_count" => $majorAlleleCount,
           "minor_allele_count" => $minorAlleleCount,
           "major_product" => $majorProduct,
           "minor_product" => $minorProduct,
           "distinct_strain_count" => $distinctStrainCount,
           "distinct_allele_count" => $distinctAlleleCount,
	     "is_coding" => $referenceVariation->{is_coding},
	     "positions_in_cds" => $referenceVariation->{positions_in_cds},
	     "positions_in_protein" => $referenceVariation->{positions_in_protein},
	     "reference_aa_full" => $referenceVariation->{products},
           "total_allele_count" => $totalAlleleCount,
           "has_stop_codon" => $hasStopCodon
    };

}

sub hashRefToString {
    my ($hashRef) = @_;

    return unless(keys %$hashRef);

    my @a;
    while(my ($key, $value) = each %$hashRef) {
	push @a, "'$key'=>'$value'";
    }

    return "{" . join(",", @a) . "}";
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

      #print STDERR Dumper $fileReader->{_coverage_array};
      #print STDERR Dumper $fileReader->{_percents_array};

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

sub querySequenceSubstring {
    my ($dbh, $sequenceId, $start, $end) = @_;

    my $length = $end - $start + 1;

    $SEQUENCE_QUERY_SH->execute($start, $length, $sequenceId);
    my ($base) = $SEQUENCE_QUERY_SH->fetchrow_array();

    $SEQUENCE_QUERY_SH->finish();
    return $base;
}

sub variationProduct {
    my ($extDbRlsId, $transcripts, $transcriptSummary, $sequenceId, $location, $positionsInCds, $allele, $refExtDbRlsId, $refSequenceId) = @_;

    my %products;
    my %refProducts;

    my $adjacentSnpCausesProductDifference = 0;

    foreach my $transcript (@$transcripts) {

	my $consensusCodingSequence;
	if($transcriptSummary->{$transcript}->{cache}->{$extDbRlsId}->{consensus_cds}) {
	    $consensusCodingSequence = $transcriptSummary->{$transcript}->{cache}->{$extDbRlsId}->{consensus_cds};
	}
	else { # first time through for this transcript
	    $consensusCodingSequence = &getCodingSequence($dbh, $sequenceId, $transcriptSummary, $transcript, $location, $location, $extDbRlsId);
	    $transcriptSummary->{$transcript}->{cache}->{$extDbRlsId}->{consensus_cds} = $consensusCodingSequence;
	}


	my $refConsensusCodingSequence;
	if($transcriptSummary->{$transcript}->{cache}->{$refExtDbRlsId}->{consensus_cds}) {
	    $refConsensusCodingSequence = $transcriptSummary->{$transcript}->{cache}->{$refExtDbRlsId}->{consensus_cds};
	}
	else { # first time through for this transcript
	    $refConsensusCodingSequence = &getCodingSequence($dbh, $refSequenceId, $transcriptSummary, $transcript, $location, $location, $refExtDbRlsId);
	    $transcriptSummary->{$transcript}->{cache}->{$refExtDbRlsId}->{consensus_cds} = $refConsensusCodingSequence;
	}

	my $strand = $transcriptSummary->{$transcript}->{cds_strand};

	my $positionInCds = $positionsInCds->{$transcript};

	next unless($positionInCds);

	next if($positionInCds > length $consensusCodingSequence);

	my $product = &getAminoAcidSequenceOfSnp($consensusCodingSequence, $positionInCds);
	my $refProduct = &getAminoAcidSequenceOfSnp($refConsensusCodingSequence, $positionInCds);

	if($product ne $refProduct) {
	    $adjacentSnpCausesProductDifference = 1;
	}

	$products{$transcript} = $product;
    }

    return(\%products, $adjacentSnpCausesProductDifference);
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


sub queryExtDbRlsAndProtocolAppNodeIdsForStrains {
    my ($strains, $dbh, $organismAbbrev) = @_;

  my $sql = "select d.name, d.external_database_name, d.version
from apidb.organism o, apidb.datasource d 
where o.abbrev = ?
and o.taxon_id = d.taxon_id
and d.name like ?
";

    my $sh = $dbh->prepare($sql);

    my %rv;

    foreach my $strain (@$strains) {

	my $match = "\%_${strain}_HTS_SNPSample_RSRC";
	$sh->execute($organismAbbrev, $match);

	my $ct;
	while(my ($name, $extDbName, $extDbVersion) = $sh->fetchrow_array()) {
	    my $spec = "$extDbName|$extDbVersion";
	    my $extDbRlsId = &queryExtDbRlsIdFromSpec($dbh, $spec);

	    $rv{$strain}->{ext_db_rls_id} = $extDbRlsId;
	    my $protocolAppNodeId = &queryProtocolAppNodeIdFromExtDbRlsId($dbh, $extDbRlsId);
	    $rv{$strain}->{protocol_app_node_id} = $protocolAppNodeId;
	    $ct++;
	}

	$sh->finish();
	die "Expected Exactly one hts snp sample row for organism=$organismAbbrev and strain=$strain" unless $ct = 1;
    }

    return \%rv;
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

		$reader = VEuPath::FileReader->new("zcat $fullPath |", [], qr/\t/);
	    } 
	    else {
		$reader = VEuPath::FileReader->new($fullPath, [], qr/\t/);
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


sub queryExtDbRlsIdFromSpec {
    my ($dbh, $extDbRlsSpec) = @_;

  my $sql = "select r.external_database_release_id
from sres.externaldatabaserelease r, sres.externaldatabase d
where d.external_database_id = r.external_database_id
and d.name || '|' || r.version = '$extDbRlsSpec'";


    my $sh = $dbh->prepare($sql);
    $sh->execute();

    my ($extDbRlsId) = $sh->fetchrow_array();

    $sh->finish();
    die "Could not find ext db rls id for spec: $extDbRlsSpec" unless($extDbRlsId);

    return $extDbRlsId;
}

sub queryProtocolAppNodeIdFromExtDbRlsId {
    my ($dbh, $extDbRlsId) = @_;

    my $sql = "select protocol_app_node_id from study.protocolappnode where external_database_release_id = $extDbRlsId and name like '%Sequence Variation%'";

    my $sh = $dbh->prepare($sql);
    $sh->execute();

    my ($panId) = $sh->fetchrow_array();

    $sh->finish();

    return $panId;
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


sub getGeneLocations {
    my ($transcriptSummary) = @_;

    my %geneSummary;
    my %geneLocations;
    my %sortedGeneLocs;

    foreach my $transcriptId (keys %$transcriptSummary) {
	my $geneId = $transcriptSummary->{$transcriptId}->{gene_na_feature_id};
	my $sequenceSourceId = $transcriptSummary->{$transcriptId}->{sequence_source_id};
	my $transcriptMinStart = $transcriptSummary->{$transcriptId}->{min_exon_start};
	my $transcriptMaxEnd = $transcriptSummary->{$transcriptId}->{max_exon_end};

	if(!{$geneSummary{$geneId}->{max_end}} || $transcriptMaxEnd > $geneSummary{$geneId}->{max_end}) {
	    $geneSummary{$geneId}->{max_end} = $transcriptMaxEnd;
	}

	if(!{$geneSummary{$geneId}->{min_start}} || $transcriptMinStart > $geneSummary{$geneId}->{min_start}) {
	    $geneSummary{$geneId}->{min_start} = $transcriptMinStart;
	}

	$geneSummary{$geneId}->{sequence_source_id} = $sequenceSourceId;
	push @{$geneSummary{$geneId}->{transcripts}}, $transcriptId;
    }

    foreach my $geneId (keys %geneSummary) {
	my $sequenceSourceId = $geneSummary{$geneId}->{sequence_source_id};

	my $location = { transcripts => $geneSummary{$geneId}->{transcripts},
			 start => $geneSummary{$geneId}->{min_start},
			 end => $geneSummary{$geneId}->{max_end},
                     na_feature_id => $geneId,
	};

	push(@{$geneLocations{$sequenceSourceId}}, $location);

    }

    foreach my $seqId (keys %geneLocations) {
	my @sortedLocations = sort { $a->{start} <=> $b->{start} } @{$geneLocations{$seqId}};
	push @{$sortedGeneLocs{$seqId}}, @sortedLocations;
    }

    return \%sortedGeneLocs;
}
  

sub getTranscriptLocations {
    my ($dbh, $transcriptExtDbRlsId, $agpMap) = @_;

    my %transcriptSummary;

my $sql = "select listagg(taf.source_id, ',') WITHIN GROUP (ORDER BY taf.aa_feature_id) as transcripts,
       s.source_id, 
       tf.parent_id as gene_na_feature_id, 
       el.start_min as exon_start, 
       el.end_max as exon_end,
       decode(el.is_reversed, 1, afe.coding_end, afe.coding_start) as cds_start,
       decode(el.is_reversed, 1, afe.coding_start, afe.coding_end) as cds_end,
       el.is_reversed
from dots.transcript tf
   , dots.translatedaafeature taf
   , dots.aafeatureexon afe
   , dots.exonfeature ef
   , dots.nalocation el
   , dots.nasequence s
where tf.na_feature_id = taf.na_feature_id
 and taf.aa_feature_id = afe.aa_feature_id
 and afe.exon_feature_id = ef.na_feature_id
 and ef.na_feature_id = el.na_feature_id
 and ef.na_sequence_id = s.na_sequence_id
 and tf.external_database_release_id = $transcriptExtDbRlsId
group by s.source_id, tf.parent_id, el.start_min, el.end_max, afe.coding_start, afe.coding_end, el.is_reversed
order by s.source_id, el.start_min
";

    my $sh = $dbh->prepare($sql);
    $sh->execute();

    while(my ($transcripts, $sequenceSourceId, $geneNaFeatureId, $exonStart, $exonEnd, $cdsStart, $cdsEnd, $isReversed) = $sh->fetchrow_array()) {
	my @transcripts = split(",", $transcripts);

	my $strand = $isReversed ? -1 : +1;

	my $agpArray = $agpMap->{$sequenceSourceId};

	if($agpArray) {
	    my $agp = VEuPath::GeneModelLocations::findAgpFromPieceLocation($agpArray, $exonStart, $exonEnd, $sequenceSourceId);

      # if this sequence is a PIECE in another sequence... lookup the higher level sequence
	    if($agp) {
        my $exonMatch = Bio::Location::Simple->
            new( -seq_id => 'exon', -start => $exonStart  , -end => $exonEnd , -strand => $strand );

        if($cdsStart && $cdsEnd) {
          my $cdsMatch = Bio::Location::Simple->
              new( -seq_id => 'cds', -start => $cdsStart  , -end => $cdsEnd , -strand => $strand );
          my $cdsMatchOnVirtual = $agp->map( $cdsMatch );
          $cdsStart = $cdsMatchOnVirtual->start();
          $cdsEnd = $cdsMatchOnVirtual->end();
        }

        my $matchOnVirtual = $agp->map( $exonMatch );
     
        $sequenceSourceId = $matchOnVirtual->seq_id();
        $exonStart = $matchOnVirtual->start();
        $exonEnd = $matchOnVirtual->end();
        $strand = $matchOnVirtual->strand();
	    }
	}


	my $loc = Bio::Location::Simple->new( -seq_id => $sequenceSourceId, -start => $exonStart  , -end => $exonEnd , -strand => $strand);


	foreach my $transcriptId (@transcripts) {
	    push @{$transcriptSummary{$transcriptId}->{exons}}, $loc;
	    $transcriptSummary{$transcriptId}->{gene_na_feature_id} = $geneNaFeatureId;
	    $transcriptSummary{$transcriptId}->{sequence_source_id} = $sequenceSourceId;

	    if($cdsStart && $cdsEnd) {
		$transcriptSummary{$transcriptId}->{cds_strand} = $strand;
		if(!{$transcriptSummary{$transcriptId}->{max_cds_end}} || $cdsEnd > $transcriptSummary{$transcriptId}->{max_cds_end}) {
		    $transcriptSummary{$transcriptId}->{max_cds_end} = $cdsEnd;
		}
		if(!$transcriptSummary{$transcriptId}->{min_cds_start} || $cdsStart < $transcriptSummary{$transcriptId}->{min_cds_start}) {
		    $transcriptSummary{$transcriptId}->{min_cds_start} = $cdsStart;
		}
	    }

	    if(!$transcriptSummary{$transcriptId}->{max_exon_end} || $exonEnd > $transcriptSummary{$transcriptId}->{max_exon_end}) {
		$transcriptSummary{$transcriptId}->{max_exon_end} = $exonEnd;
	    }

	    if(!$transcriptSummary{$transcriptId}->{min_exon_start} || $exonStart < $transcriptSummary{$transcriptId}->{min_exon_start}) {
		$transcriptSummary{$transcriptId}->{min_exon_start} = $exonStart;
	    }

	}
    }

    $sh->finish();

    return \%transcriptSummary;
}


sub getCodingSequence {
    my ($dbh, $sequenceId, $transcriptSummary, $transcriptId, $snpStart, $snpEnd, $seqExtDbRlsId) = @_;

    return unless($transcriptId);

    my @exons = @{$transcriptSummary->{$transcriptId}->{exons}};
    my $minCodingStart = $transcriptSummary->{$transcriptId}->{min_cds_start};
    my $maxCodingEnd = $transcriptSummary->{$transcriptId}->{max_cds_end};

    unless (@exons) {
	die ("Transcript with na_feature_id = $transcriptId had no exons\n");
    }

    my $codingSequence;

  # sort exons by start_min
    for my $exon (sort {$a->{_start} <=> $b->{_start} } @exons) {
	my $exonStart = $exon->start();
	my $exonEnd = $exon->end();
	my $exonIsReversed = $exon->strand() == -1 ? 1 : 0;

	my $codingMin = $minCodingStart >= $exonStart ? $minCodingStart : $exonStart;
	my $codingMax = $maxCodingEnd <= $exonEnd ? $maxCodingEnd : $exonEnd;

	my $chunk =   &querySequenceSubstring($dbh, $sequenceId, $codingMin, $codingMax);

	if($exonIsReversed) {
	    $chunk = CBIL::Bio::SequenceUtils::reverseComplementSequence($chunk);
	    $codingSequence = $chunk . $codingSequence;
	}
	else {
	    $codingSequence .= $chunk;
	}

    }

    return($codingSequence);
}


sub calculateAminoAcidPosition {
    my ($codingPosition) = @_;

    my $aaPos = ($codingPosition % 3 == 0) ? int($codingPosition / 3) : int($codingPosition / 3) + 1;

    return($aaPos);
}



sub getAminoAcidSequenceOfSnp {
    my ($cdsSequence, $positionInCds) = @_;


    my $codonLength = 3;
    my $modCds = ($positionInCds - 1)  % $codonLength;
    my $offset = $positionInCds - $modCds;

    my $codon = substr $cdsSequence, $offset - 1, $codonLength;

  # Assuming this is a simple hash lookup which is quick; 
    return $CODON_TABLE->translate($codon);
}

sub createCurrentShifts {
    my ($strains, $dbh) = @_;
    my $currentShifts;
    foreach my $strain (@$strains) {
        my $INDEL_QUERY = "SELECT i.location, i.shift FROM apidb.indel i WHERE sample_name = '$strain'";
        my $INDEL_QUERY_SH = $dbh->prepare($INDEL_QUERY);
        $INDEL_QUERY_SH->execute();
        my @locationshifts = ();
	my $counter = 0;
	my $currentShift = 0;
	while (my ($location, $shift) = $INDEL_QUERY_SH->fetchrow_array()) {
	    push ( @{$locationshifts[$counter]}, ($location, $shift + $currentShift));
	    $counter++;
	    $currentShift = $shift + $currentShift;
	}
        push @{ $currentShifts->{$strain} }, \@locationshifts;
    }
    return $currentShifts;
}


sub calcCoordinates {
    my ($shiftFrame, $shiftFrameLimit, $oldShift, $coordinate, $indicator, $shiftArray) = @_;
    my $shiftedLocation;
    my $oldFrame;
    #print "$coordinate\t$shiftArray[0][$shiftFrame][0]\t$shiftFrame\n";
    if ($coordinate < $shiftArray[0][0][$shiftFrame][0]) {
	$shiftedLocation = $oldShift + $coordinate;
        #print "Coordinate Less\n";
    }
    elsif ($shiftArray[0][0][$shiftFrame][0] == $coordinate) {
        my $currentShift = $shiftArray[$shiftFrame][1];
        #print "Coordinate Equal\n";
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
    elsif ($coordinate > $shiftArray[0][0][$shiftFrame][0] || $shiftFrame == $shiftFrameLimit) {
        #print "Coordinate Greater\t$coordinate\t$shiftArray[0][$shiftFrame][0]\n";
	until ($shiftArray[0][0][$shiftFrame][0] >= $coordinate || $shiftFrame == $shiftFrameLimit) {
	    $oldShift = $shiftArray[0][0][$shiftFrame][1];
	    $shiftFrame++;
            #print "DING\n";
	}
        #print "Processed\t$coordinate\t$shiftArray[0][$shiftFrame][0]\n";
	if ($shiftFrame == $shiftFrameLimit && $coordinate < $shiftArray[0][0][$shiftFrame][0]) {
	    $shiftedLocation = $coordinate + $shiftArray[0][0][$shiftFrame-1][1];
	}
	elsif ($shiftFrame == $shiftFrameLimit && $coordinate > $shiftArray[0][0][$shiftFrame][0]) {
	    $shiftedLocation = $coordinate + $shiftArray[0][0][$shiftFrame][1];
	}
	elsif ($shiftArray[0][0][$shiftFrame][0] == $coordinate) {
	    if ($shiftArray[0][0][$shiftFrame][1] == 0) {
                $shiftedLocation = $coordinate;
            }
            elsif ($indicator eq 'start') {
		$oldFrame = $shiftFrame - 1;
                $shiftedLocation = $shiftArray[0][0][$oldFrame][1] + $coordinate;
            }
            elsif ($indicator eq 'end' && $shiftArray[0][0][$shiftFrame][1] > 0) {
                $shiftedLocation = $shiftArray[0][0][$shiftFrame][1] + $coordinate;     
            }
	    elsif ($indicator eq 'end' && $shiftArray[0][0][$shiftFrame][1] < 0) {
                $oldFrame = $shiftFrame - 1;
                $shiftedLocation = $shiftArray[0][0][$oldFrame][1] + $coordinate;     
            }
        }
	else {
	    $shiftedLocation = $oldShift + $coordinate;
	}
    }
    return ($shiftedLocation, $oldShift, $shiftFrame);   
}

sub addStrainExonShiftsToTranscriptSummary {
    my ($currentShifts, $transcriptSummary) = @_;
    my ($oldShift, $shiftFrame);
    foreach my $strain (keys %{ $currentShifts }) {
        my @shiftArray = $currentShifts->{$strain};
	my $indexedArray = $shiftArray[0][0];
        my $shiftArrayLen = scalar @{ $indexedArray };
        my $shiftFrameLimit = $shiftArrayLen - 1;
        $oldShift = 0;
        $shiftFrame = 0;
        my ($exon_start, $exon_end);
        my $startIndicator = "start";
        my $endIndicator = "end";
	foreach my $transcript (keys %{ $transcriptSummary }) {    
            my ($shifted_start, $shifted_end);
	    $exon_start = $$transcriptSummary{$transcript}->{min_exon_start};
	    $exon_end = $$transcriptSummary{$transcript}->{max_exon_end};
	    ($shifted_start, $shiftFrame, $oldShift) = &calcCoordinates($shiftFrame, $shiftFrameLimit, $oldShift, $exon_start, $startIndicator, $shiftArray);
	    ($shifted_end, $shiftFrame, $oldShift) = &calcCoordinates($shiftFrame, $shiftFrameLimit, $oldShift, $exon_end, $endIndicator, $shiftArray);    
	    #print "Values $exon_start\t$exon_end\t$shifted_start\t$shifted_end\t$oldShift\t$shiftFrame\n";
	    $$transcriptSummary{$transcript}->{$strain}->{shifted_start} = $shifted_start;
	    $$transcriptSummary{$transcript}->{$strain}->{shifted_end} = $shifted_end;
	}
    }
    return $transcriptSummary;
}

