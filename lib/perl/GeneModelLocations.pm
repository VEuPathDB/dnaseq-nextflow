package VEuPath::GeneModelLocations;

use strict;
use Bio::Location::Simple;
use Bio::Coordinate::Pair;
use Bio::Coordinate::GeneMapper;
use Data::Dumper;

BEGIN {

  my $geneSh;
  my $seqSh;
  my $virtualSeqSh;

  sub setGeneSh {
    my ($sh) = @_;

    $geneSh = $sh;
  }

  sub getGeneSh {
    return $geneSh;
  }

  sub setSeqSh {
    my ($sh) = @_;

    $seqSh = $sh;
  }

  sub getSeqSh {
    return $seqSh;
  }

  sub setVirtualSeqSh {
    my ($sh) = @_;

    $virtualSeqSh = $sh;
  }

  sub getVirtualSeqSh {
    return $virtualSeqSh;
  }
}

sub new {
  my ($class, $dbh, $geneExtDbRlsId, $wantTopLevel, $optionalSoTerm, $soExclude) = @_;

  my $self = bless {'_database_handle' => $dbh, 
                    '_gene_external_database_release_id' => $geneExtDbRlsId, 
                    '_want_top_level' => $wantTopLevel,
                    '_sequence_ontology_term' => $optionalSoTerm,
                    '_sequence_ontology_exclude' => $soExclude,
  }, $class;

  my $agpMap = {};

  $agpMap = $self->queryForAgpMap($dbh);


  $self->setAgpMap($agpMap);

  $self->_initAllModelsHash();
  $self->_setTranscriptLocations();

  return $self;
}

#--------------------------------------------------------------------------------
# public methods start
#--------------------------------------------------------------------------------

sub getDatabaseHandle { $_[0]->{_database_handle} }
sub getGeneExternalDatabaseReleaseId { $_[0]->{_gene_external_database_release_id} }        
sub getWantTopLevel { $_[0]->{_want_top_level} }

sub getAgpMap { $_[0]->{_agp_map} }
sub setAgpMap { $_[0]->{_agp_map} = $_[1] }

sub getSequenceOntologyTerm { $_[0]->{_sequence_ontology_term} }
sub setSequenceOntologyTerm { $_[0]->{_sequence_ontology_term} = $_[1] }

sub getSequenceOntologyExclude { $_[0]->{_sequence_ontology_exclude} }
sub setSequenceOntologyExclude { $_[0]->{_sequence_ontology_exclude} = $_[1] }


sub getAllGeneIds {
  my ($self) = @_;

  my @rv = keys %{$self->_getAllModelsHash()};

  return \@rv;
}

sub bioperlFeaturesFromGeneSourceId {
  my ($self, $geneSourceId) = @_;

  my @rv;

  my $sourceTag = 'GUS';

  my %modToPhase = (0 => 0,
                    1 => 2,
                    2 => 1
      );


  my $geneModelHash = $self->getGeneModelHashFromGeneSourceId($geneSourceId);

  my $geneFeature = GUS::Community::GeneModelLocations::Gene->new( -start => $geneModelHash->{start}, 
                                                   -end => $geneModelHash->{end}, 
                                                   -seq_id => $geneModelHash->{sequence_source_id},
                                                   -strand => $geneModelHash->{strand}, 
                                                   -primary => $geneModelHash->{sequence_ontology_term},
                                                   -source_tag => $sourceTag, 
                                                   -tag    => { ID => $geneSourceId,
                                                                NA_FEATURE_ID => $geneModelHash->{na_feature_id},
                                                                NA_SEQUENCE_ID => $geneModelHash->{na_sequence_id},
                                                                SEQUENCE_IS_PIECE => $geneModelHash->{sequence_is_piece}, 
                                                                EBI_BIOTYPE => $geneModelHash->{ebi_biotype},
                                                   });

  push @rv, $geneFeature;

  my @cdsFeatures;
  my @utrFeatures;

  my %transcriptMap = ();

  foreach my $transcriptSourceId (keys %{$geneModelHash->{transcripts}}) {
    my $transcriptHash = $geneModelHash->{transcripts}->{$transcriptSourceId};

#    my $exonSourceIds = $transcriptHash->{exonSourceIds};

#    my $transcriptStart = $geneModelHash->{end};
 #   my $transcriptEnd = $geneModelHash->{start};

 #   foreach my $exonSourceId (@$exonSourceIds) {
 #     my $exonStart = $geneModelHash->{exons}->{$exonSourceId}->{start};
 #     my $exonEnd = $geneModelHash->{exons}->{$exonSourceId}->{end};
 #     $transcriptStart =  $exonStart if($exonStart < $transcriptStart);
 #     $transcriptEnd =  $exonEnd if($exonEnd > $transcriptEnd);
 #   }

    my $transcriptFeature = GUS::Community::GeneModelLocations::Transcript->new ( -start => $transcriptHash->{start},
                                                           -end => $transcriptHash->{end},
                                                   -seq_id => $geneModelHash->{sequence_source_id},
                                                   -strand => $transcriptHash->{gene_strand}, 
                                                   -primary => $transcriptHash->{sequence_ontology_term},
                                                   -source_tag => $sourceTag, 
                                                   -tag    => { ID => $transcriptSourceId,
                                                                NA_FEATURE_ID => $transcriptHash->{na_feature_id},
                                                                PARENT => $geneModelHash->{source_id},
                                                                PARENT_NA_FEATURE_ID => $geneModelHash->{na_feature_id},
                                                                NA_SEQUENCE_ID => $geneModelHash->{na_sequence_id},
                                                                GENE_EBI_BIOTYPE => $geneModelHash->{ebi_biotype},
                                                   });

    $transcriptMap{$transcriptSourceId} = $transcriptFeature;

    push @rv, $transcriptFeature;

    my $cdsCount = 1;

    foreach my $proteinSourceId (keys %{$transcriptHash->{proteins}}) {
      my $proteinHash = $transcriptHash->{proteins}->{$proteinSourceId};

      my $prevCdsLength;

      my $minCdsLoc = $transcriptHash->{end};
      my $maxCdsLoc = $transcriptHash->{start};


      my @utrs = ();

      foreach my $pExon (sort { $a->{exon_start} * $a->{strand} <=> $b->{exon_start} * $b->{strand}} @{$proteinHash->{exons}}) {
        if($pExon->{is_all_utr}) {
          push @utrs, [$pExon->{exon_start}, $pExon->{exon_end}, $pExon->{strand}];
          next;
        }

        my $phase;
        if($prevCdsLength) {
          my $mod = $prevCdsLength % 3;
          $phase = $modToPhase{$mod};
        }
        else {
          $phase = 0;
        }

        $prevCdsLength += $pExon->{cds_end} - $pExon->{cds_start} + 1;

        $minCdsLoc = $pExon->{cds_start} if($pExon->{cds_start}) < $minCdsLoc;
        $maxCdsLoc = $pExon->{cds_end} if($pExon->{cds_end}) > $maxCdsLoc;


        my $cdsId = "${proteinSourceId}-CDS$cdsCount";

        my $cdsFeature = GUS::Community::GeneModelLocations::CDS->new ( -start => $pExon->{cds_start},
                                                        -end => $pExon->{cds_end}, 
                                                        -seq_id => $geneModelHash->{sequence_source_id},
                                                        -strand => $pExon->{strand},
                                                        -frame => $phase,
                                                        -primary => 'CDS',
                                                        -source_tag => $sourceTag, 
                                                        -tag    => { ID => $cdsId,
                                                                     AA_FEATURE_ID => $proteinHash->{aa_feature_id},
                                                                     AA_SEQUENCE_ID => $proteinHash->{aa_sequence_id},
                                                                     PARENT => $transcriptSourceId,
                                                                     NA_SEQUENCE_ID => $geneModelHash->{na_sequence_id},
                                                                     PROTEIN_SOURCE_ID => $proteinSourceId,
                                                                     PARENT_NA_FEATURE_ID => $transcriptHash->{na_feature_id}
                                                        });


        if($pExon->{cds_start} > $pExon->{exon_start}) {
          push @utrs, [$pExon->{exon_start}, $pExon->{cds_start} - 1, $pExon->{strand}];
        }
        if($pExon->{cds_end} < $pExon->{exon_end}) {
          push @utrs, [$pExon->{cds_end} + 1, $pExon->{exon_end}, $pExon->{strand}];
        }


        push @cdsFeatures, $cdsFeature;
        $cdsCount++;
      }

      my $utrCount;
      foreach my $utr (sort { $a->[0] <=> $b->[0] } @utrs) {
        my $utrStart = $utr->[0];
        my $utrEnd = $utr->[1];
        my $utrStrand = $utr->[2];


        my $utrDirection = &utrDirection($minCdsLoc, $maxCdsLoc, $utrStart, $utrEnd, $utrStrand, $transcriptSourceId);

        $utrCount++;
        my $utrId = "utr_${transcriptSourceId}_$utrCount";
        my $utrFeature = GUS::Community::GeneModelLocations::UTR->new ( -start => $utr->[0],
                                                          -end => $utr->[1],
                                                          -strand => $utr->[2],
                                                          -source_tag => $sourceTag,
                                                          -primary => $utrDirection,
                                                          -seq_id => $geneModelHash->{sequence_source_id},
                                                          -tag    => { ID => $utrId,
                                                                       PARENT => $transcriptSourceId,
                                                                       PARENT_NA_FEATURE_ID => $transcriptHash->{na_feature_id},
                                                                       NA_SEQUENCE_ID => $geneModelHash->{na_sequence_id},
                                                          });
        push @utrFeatures, $utrFeature;
      }
    }
  }


  my @exonFeatures;

  foreach my $exonHash (sort { $a->{start} * $a->{strand} <=> $b->{start} * $b->{strand} } values %{$geneModelHash->{exons}}) {
    die "exon strand cannot be set to 0" if($exonHash->{strand} == 0);

    my $parent = join(",", @{$exonHash->{transcripts}});


    my $exonFeature = GUS::Community::GeneModelLocations::Exon->new ( -start => $exonHash->{start}, 
                                                   -end => $exonHash->{end}, 
                                                   -seq_id => $exonHash->{sequence_source_id},
                                                   -strand => $exonHash->{strand}, 
                                                     -primary => 'exon',
                                                   -source_tag => $sourceTag, 
                                                   -tag    => { ID => $exonHash->{source_id},
                                                                NA_FEATURE_ID => $exonHash->{na_feature_id},
                                                                PARENT => $parent,
                                                                GENE_NA_FEATURE_ID => $geneModelHash->{na_feature_id},
                                                                NA_SEQUENCE_ID => $geneModelHash->{na_sequence_id},
                                                   });


    foreach my $transcriptId (@{$exonHash->{transcripts}}) {
      my $transcriptFeature = $transcriptMap{$transcriptId};
      $transcriptFeature->add_exon($exonFeature);
    }

    push @exonFeatures, $exonFeature
  }

  my @sortedExonFeatures = sort { $a->start <=> $b->start} @exonFeatures;
  push @rv, @sortedExonFeatures;

  my @sortedCdsFeatures = sort { $a->start <=> $b->start} @cdsFeatures;
  push @rv, @sortedCdsFeatures;

  my @sortedUtrFeatures = sort { $a->start <=> $b->start} @utrFeatures;
  push @rv, @sortedUtrFeatures;


  return \@rv;
}


sub utrDirection {
  my ($minCdsLoc, $maxCdsLoc, $utrStart, $utrEnd, $utrStrand, $transcriptSourceId) = @_;


  # utr5prime and utr3prime are what the bioperl object expects for primary_tag (although not what gff3 wants argh;  must deal with that downstream)
  if($utrStrand == -1) {
    if($utrStart <= $minCdsLoc && $utrEnd <= $minCdsLoc) {
      return 'utr3prime';
    }
    if($utrStart >= $maxCdsLoc && $utrEnd >= $maxCdsLoc) {
      return 'utr5prime';
    }
  }
  else {
    if($utrStart <= $minCdsLoc && $utrEnd <= $minCdsLoc) {
      return 'utr5prime';
    }
    if($utrStart >= $maxCdsLoc && $utrEnd >= $maxCdsLoc) {
      return 'utr3prime';
    }
  }
  die "UTR [$utrStart, $utrEnd, $utrStrand] for transcript $transcriptSourceId fell w/in cds coords [$minCdsLoc, $maxCdsLoc]";
}


sub getGeneModelHashFromGeneSourceId {
  my ($self, $geneSourceId) = @_;

  return $self->_getAllModelsHash()->{$geneSourceId};
}

sub getTranscriptIdsFromGeneSourceId {
  my ($self, $geneSourceId) = @_;

  my @transcriptIds  = keys %{$self->_getAllModelsHash()->{$geneSourceId}->{transcripts}};

  return \@transcriptIds;
}

sub getTranscriptHashFromTranscriptId {
  my ($self, $transcriptId) = @_;

  my $geneSourceId = $self->_getTranscriptToGeneMap()->{$transcriptId};

  my $geneHash = $self->getGeneModelHashFromGeneSourceId($geneSourceId);

  return $geneHash->{transcripts}->{$transcriptId};
}


sub getProteinIdsFromTranscriptSourceId {
  my ($self, $transcriptSourceId) = @_;

  my $transcriptHash = $self->getTranscriptHashFromTranscriptId($transcriptSourceId);

  my @proteinIds  = keys %{$transcriptHash->{proteins}};

  return \@proteinIds;
}


sub getProteinHashFromProteinId {
  my ($self, $proteinId) = @_;

  my $transcriptId = $self->_getProteinToTranscriptMap->{$proteinId};

  my $transcriptHash = $self->getTranscriptHashFromTranscriptId($transcriptId);

  return $transcriptHash->{proteins}->{$proteinId};
}



sub getProteinToGenomicCoordMapper {
  my ($self, $proteinId) = @_;

  my ($exonLocs, $cdsRange, $exonRange) = $self->getExonLocsAndCdsRangeFromProteinId($proteinId);

  my $mapper = Bio::Coordinate::GeneMapper->new(
    -in    => 'peptide',
    -out   => 'chr',
    -exons => $exonLocs,
    -cds => $cdsRange
      );

  return $mapper;
}

sub getExonLocsAndCdsRangeFromProteinId {
  my ($self, $proteinId) = @_;

  my $proteinHash = $self->getProteinHashFromProteinId($proteinId);

  my $naSequenceSourceId = $proteinHash->{na_sequence_source_id};

  my ($minCds, $maxCds, $minExon, $maxExon, $strand);

  my @exonLocs;

  foreach my $exon (@{$proteinHash->{exons}}) {
    # init min and max w/ first value
    $minCds = $exon->{cds_start} unless($minCds);
    $maxCds = $exon->{cds_start} unless($maxCds);

    $minExon = $exon->{exon_start} unless($minExon);
    $maxExon = $exon->{exon_start} unless($maxExon);


    my ($min, $max) = sort {$a <=> $b} ($exon->{cds_start}, $exon->{cds_end});

    my ($minE, $maxE) = sort {$a <=> $b} ($exon->{exon_start}, $exon->{exon_end});

    $minCds = $min if($min < $minCds);
    $maxCds = $max if($max > $maxCds);

    $minExon = $minE if($minE < $minExon);
    $maxExon = $maxE if($maxE > $maxExon);


    #used for exon and cds
    $strand = $exon->{strand};

    my $exonLoc = Bio::Location::Simple->new( -seq_id => $naSequenceSourceId, 
                                              -start => $exon->{exon_start}, 
                                              -end => $exon->{exon_end}, 
                                              -strand => $strand);  

    push @exonLocs, $exonLoc;
  }

  my $cdsRange = Bio::Location::Simple->new( -seq_id => $naSequenceSourceId, 
                                             -start => $minCds, 
                                             -end => $maxCds, 
                                             -strand => $strand);  


  my $exonRange = Bio::Location::Simple->new( -seq_id => $naSequenceSourceId, 
                                              -start => $minExon, 
                                              -end => $maxExon, 
                                              -strand => $strand);  



  return(\@exonLocs, $cdsRange, $exonRange);
}

sub getCdsToGenomicCoordMapper {
  my ($self, $proteinId) = @_;

  my ($exonLocs, $cdsRange, $exonRange) = $self->getExonLocsAndCdsRangeFromProteinId($proteinId);

  my $mapper = Bio::Coordinate::GeneMapper->new(
    -in    => 'cds',
    -out   => 'chr',
    -exons => $exonLocs,
    -cds => $cdsRange
      );

  return $mapper;
}



sub getTranscriptToGenomicCoordMapper {
  my ($self, $proteinId) = @_;

  my ($exonLocs, $cdsRange, $exonRange) = $self->getExonLocsAndCdsRangeFromProteinId($proteinId);

  my $mapper = Bio::Coordinate::GeneMapper->new(
    -in    => 'cds',
    -out   => 'chr',
    -exons => $exonLocs,
    -cds => $exonRange
      );

  return $mapper;
}

sub getGenomicToTranscriptCoordMapper {
  my ($self, $proteinId) = @_;

  my ($exonLocs, $cdsRange, $exonRange) = $self->getExonLocsAndCdsRangeFromProteinId($proteinId);

  my $mapper = Bio::Coordinate::GeneMapper->new(
    -in    => 'chr',
    -out   => 'cds',
    -exons => $exonLocs,
    -cds => $exonRange
      );

  return $mapper;
}







#--------------------------------------------------------------------------------
# private methods
#--------------------------------------------------------------------------------

sub _setTranscriptLocations {
  my ($self) = @_;

  my $geneIds = $self->getAllGeneIds();

  foreach my $geneSourceId (@$geneIds) {
    my $geneModelHash = $self->getGeneModelHashFromGeneSourceId($geneSourceId);

    foreach my $transcriptSourceId (keys %{$geneModelHash->{transcripts}}) {
      my $transcriptHash = $geneModelHash->{transcripts}->{$transcriptSourceId};

      my $exonSourceIds = $transcriptHash->{exonSourceIds};

      my $transcriptStart = $geneModelHash->{end};
      my $transcriptEnd = $geneModelHash->{start};

      foreach my $exonSourceId (@$exonSourceIds) {
        my $exonStart = $geneModelHash->{exons}->{$exonSourceId}->{start};
        my $exonEnd = $geneModelHash->{exons}->{$exonSourceId}->{end};
        $transcriptStart =  $exonStart if($exonStart < $transcriptStart);
        $transcriptEnd =  $exonEnd if($exonEnd > $transcriptEnd);
      }

      $geneModelHash->{transcripts}->{$transcriptSourceId}->{start} = $transcriptStart;
      $geneModelHash->{transcripts}->{$transcriptSourceId}->{end} = $transcriptEnd;
      $geneModelHash->{transcripts}->{$transcriptSourceId}->{gene_strand} = $geneModelHash->{strand};
    }
  }
}


sub _getTranscriptToGeneMap { $_[0]->{_transcript_to_gene_map} }
sub _getProteinToTranscriptMap { $_[0]->{_protein_to_transcript_map} }
sub _getAllModelsHash { $_[0]->{_all_models_hash}}

sub _initAllModelsHash { 
  my ($self) = @_;

  my $wantTopLevel = $self->getWantTopLevel();

  my $dbh = $self->getDatabaseHandle();
  my $extDbRlsId = $self->getGeneExternalDatabaseReleaseId();

  my $sql = "select gf.source_id as gene_source_id
     , gf.na_feature_id as gene_na_feature_id
     , t.source_id as transcript_source_id
     , t.na_feature_id as transcript_na_feature_id
     , ef.source_id as exon_source_id
     , ef.na_feature_id as exon_na_feature_id
     , l.start_min as exon_start_min
     , l.end_max as exon_end_max
     , l.is_reversed as exon_is_reversed
     , p.translation_start
     , p.translation_stop
     , gfl.start_min as gene_start_min
     , gfl.end_max as gene_end_max
     , gfl.is_reversed as gene_is_reversed
     , s.source_id as sequence_source_id
     , p.protein_source_id
     , p.protein_aa_feature_id
     , p.protein_aa_sequence_id
     , p.coding_start
     , p.coding_end
     , gf.na_sequence_id 
     , so.name as gene_so_term
     , tso.name as transcript_so_term
     , s.taxon_id
     , gf.name
from dots.genefeature gf
   , dots.nalocation gfl
   , dots.nasequence s
   , dots.exonfeature ef
   , dots.nalocation l
   , dots.transcript t
   , dots.rnafeatureexon rfe
   , sres.ontologyterm so
   , sres.ontologyterm tso
   , (select taf.na_feature_id 
           , taf.translation_start
           , taf.translation_stop
           , aas.source_id as protein_source_id
           , taf.aa_feature_id as protein_aa_feature_id
           , aas.aa_sequence_id as protein_aa_sequence_id
           , afe.coding_start
           , afe.coding_end
           , afe.exon_feature_id
       from dots.translatedaafeature taf
          , dots.aafeatureexon afe
          , dots.translatedaasequence aas
       where taf.aa_feature_id = afe.aa_feature_id
       and taf.aa_sequence_id = aas.aa_sequence_id
       ) p
where gf.na_feature_id = ef.parent_id
and gf.na_feature_id = gfl.na_feature_id
and gf.na_sequence_id = s.na_sequence_id
and ef.na_feature_id = l.na_feature_id
and t.na_feature_id = rfe.rna_feature_id
and rfe.exon_feature_id = ef.na_feature_id
and t.parent_id = gf.na_feature_id
and rfe.exon_feature_id = p.exon_feature_id (+)
and rfe.rna_feature_id = p.na_feature_id (+)
and gf.external_database_release_id = ?
and gf.sequence_ontology_id = so.ontology_term_id
and t.sequence_ontology_id = tso.ontology_term_id
order by gf.na_feature_id, t.na_feature_id, l.start_min
";

  my $sh = getGeneSh();
  unless ($sh) {
    $sh = $dbh->prepare($sql);
    setGeneSh($sh);
  }
  $sh->execute($extDbRlsId);

  my $geneModels = {};

  my $transcriptToGeneMap = {};
  my $proteinToTranscriptMap = {};

  my %seenGenes;
  my %seenTranscripts;
  my %seenProteins;
  my %seenExons;
  my %seenTranscriptExons;

  my %seenGenomicSequences;

  my $agpMap = $self->getAgpMap();

  my $virtualSequenceIds = $self->{_virtual_sequence_ids};

  my $soTerm = $self->getSequenceOntologyTerm();
  my $soExclude = $self->getSequenceOntologyExclude();


  while(my $arr = $sh->fetchrow_arrayref()) {
    my $geneSourceId = $arr->[0];
    my $geneNaFeatureId = $arr->[1];
    my $transcriptSourceId = $arr->[2];
    my $transcriptNaFeatureId = $arr->[3];
    my $exonSourceId = "exon_" . $arr->[4];
    my $exonNaFeatureId = $arr->[5];
    my $exonStart = $arr->[6];
    my $exonEnd = $arr->[7];
    my $exonIsReversed = $arr->[8];
    my $translationStart = $arr->[9];
    my $translationStop = $arr->[10];
    my $geneStart = $arr->[11];
    my $geneEnd = $arr->[12];
    my $geneIsReversed = $arr->[13];
    my $sequenceSourceId = $arr->[14];
    my $proteinSourceId = $arr->[15];
    my $proteinAaFeatureId = $arr->[16];
    my $proteinAaSequenceId = $arr->[17];
    my $codingStart = $arr->[18];
    my $codingEnd = $arr->[19];
    my $naSequenceId = $arr->[20];
    my $gfSoTerm = $arr->[21];
    my $transcriptSoTerm = $arr->[22];
    my $taxonId = $arr->[23];
    my $ebiBiotype = $arr->[24];

    my $strand = $geneIsReversed ? -1 : 1;

    $seenGenomicSequences{$sequenceSourceId}++;

    unless($seenGenes{$geneSourceId}) {

      next if($soTerm && $gfSoTerm ne $soTerm && !$soExclude); # only rows for this so term
      next if($soTerm && $gfSoTerm eq $soTerm && $soExclude); # only rows which are not this so term

      my $geneLocation = &mapLocation($agpMap, $sequenceSourceId, $geneStart, $geneEnd, $strand, $wantTopLevel);

      $geneModels->{$geneSourceId} = { 'source_id' => $geneSourceId,
                                 'na_feature_id' => $geneNaFeatureId,
                                 'na_sequence_id' => $virtualSequenceIds->{$geneLocation->seq_id()} ? $virtualSequenceIds->{$geneLocation->seq_id()} : $naSequenceId,
                                 'sequence_source_id' => $geneLocation->seq_id,
                                 'start' => $geneLocation->start,
                                 'end' => $geneLocation->end,
                                 'strand' => $geneLocation->strand,
                                 'sequence_ontology_term' => $gfSoTerm,
                                 'sequence_is_piece' => $geneLocation->{_sequence_is_piece},
                                 'ebi_biotype' => $ebiBiotype,
      };


      $seenGenes{$geneSourceId} = 1;
    }


    my $seenTranscript = $seenTranscripts{$transcriptSourceId};
    my $seenExon = $seenExons{$exonSourceId};

    unless($seenTranscript) {
      $geneModels->{$geneSourceId}->{transcripts}->{$transcriptSourceId} = {'source_id' => $transcriptSourceId,
                                                                            'na_feature_id' => $transcriptNaFeatureId,
                                                                            'sequence_ontology_term' => $transcriptSoTerm,
                                                                            'na_sequence_id' => $geneModels->{$geneSourceId}->{na_sequence_id}, # GET FROM THE GENE
                                                                            'sequence_source_id' => $geneModels->{$geneSourceId}->{sequence_source_id}, # GET FROM THE GENE
      };

      $transcriptToGeneMap->{$transcriptSourceId} = $geneSourceId;
      $seenTranscripts{$transcriptSourceId} = 1;
    }



    unless($seenExon) {
      my $exonLocation = &mapLocation($agpMap, $sequenceSourceId, $exonStart, $exonEnd, $strand, $wantTopLevel);

      $geneModels->{$geneSourceId}->{exons}->{$exonSourceId} = {'source_id' => $exonSourceId,
                                                    'na_feature_id' => $exonNaFeatureId,
                                                    'start' => $exonLocation->start,
                                                    'end' => $exonLocation->end,
                                                    'strand' => $exonLocation->strand,
                                                    'sequence_source_id' => $exonLocation->seq_id,
                                                    'na_sequence_id' => $virtualSequenceIds->{$exonLocation->seq_id()} ? $virtualSequenceIds->{$exonLocation->seq_id()} : $naSequenceId,
                                                    'sequence_is_piece' => $exonLocation->{_sequence_is_piece},
      };

      $seenExons{$exonSourceId} = 1;
    }

    my $transcriptExonSourceId = "$transcriptSourceId.$exonSourceId";
    
    unless($seenTranscriptExons{$transcriptExonSourceId}) {
      push @{$geneModels->{$geneSourceId}->{exons}->{$exonSourceId}->{transcripts}}, $transcriptSourceId;
    }

    $seenTranscriptExons{$transcriptExonSourceId} = 1;

    if(defined($proteinSourceId) && !$seenProteins{$proteinSourceId}) {
      $geneModels->{$geneSourceId}->{transcripts}->{$transcriptSourceId}->{proteins}->{$proteinSourceId} = {'source_id' => $proteinSourceId,
                                                                                                            'aa_feature_id' => $proteinAaFeatureId,
                                                                                                            'aa_sequence_id' => $proteinAaSequenceId,
                                                                                                            'translation_start' => $translationStart,
                                                                                                            'translation_end' =>  $translationStop,
                                                                                                            'na_sequence_source_id' => $geneModels->{$geneSourceId}->{sequence_source_id},
                                                                                                            'na_sequence_id' => $geneModels->{$geneSourceId}->{na_sequence_id},
      };

      $proteinToTranscriptMap->{$proteinSourceId} = $transcriptSourceId;
      $seenProteins{$proteinSourceId} = 1;
    }

    push @{$geneModels->{$geneSourceId}->{transcripts}->{$transcriptSourceId}->{exonSourceIds}}, $exonSourceId;

    my $isAllUtr;
    $isAllUtr = 1 unless($codingStart && $codingEnd);

    if(defined $proteinSourceId) {
      my @sortedCds = sort {$a <=> $b} ($codingStart, $codingEnd);

      die join(' ', @sortedCds) if ($sortedCds[0] > $sortedCds[1]);

      my $cdsStart = $sortedCds[0];
      my $cdsEnd = $sortedCds[1];

      my $cdsLocation;
      $cdsLocation= &mapLocation($agpMap, $sequenceSourceId, $cdsStart, $cdsEnd, $strand, $wantTopLevel) unless($isAllUtr);

      my $exon = {source_id => $exonSourceId,
                 na_feature_id => $exonNaFeatureId,
                 strand => $geneModels->{$geneSourceId}->{exons}->{$exonSourceId}->{strand},
                 exon_start => $geneModels->{$geneSourceId}->{exons}->{$exonSourceId}->{start},
                 exon_end => $geneModels->{$geneSourceId}->{exons}->{$exonSourceId}->{end},
                 cds_start => defined($isAllUtr) ? undef : $cdsLocation->start,
                 cds_end => defined($isAllUtr) ? undef : $cdsLocation->end,
                 is_all_utr => $isAllUtr,
                  sequence_is_piece => $cdsLocation->{_sequence_is_piece},
      };


      unless(defined($isAllUtr)) {
        if($geneModels->{$geneSourceId}->{min_cds_start}) {
          $geneModels->{$geneSourceId}->{min_cds_start} = $cdsLocation->start if($cdsLocation->start < $geneModels->{$geneSourceId}->{min_cds_start});
        }
        else {
          $geneModels->{$geneSourceId}->{min_cds_start} = $cdsLocation->start;
        }

        if($geneModels->{$geneSourceId}->{max_cds_end}) {
          $geneModels->{$geneSourceId}->{max_cds_end} = $cdsLocation->end if($cdsLocation->end > $geneModels->{$geneSourceId}->{max_cds_end});
        }
        else {
          $geneModels->{$geneSourceId}->{max_cds_end} = $cdsLocation->end;
        }
      }

      push @{$geneModels->{$geneSourceId}->{transcripts}->{$transcriptSourceId}->{proteins}->{$proteinSourceId}->{exons}}, $exon;
    }
  }


  my @sortedGenes = map { [$_->{sequence_source_id}, $_->{source_id}, $_->{start}, $_->{end}, $_->{min_cds_start}, $_->{max_cds_end}] } sort {$a->{sequence_source_id} cmp $b->{sequence_source_id} || $a->{start} <=> $b->{start} } values %{$geneModels};;

  my @distinctSequenceSourceIds = keys %seenGenomicSequences;
  my $sequenceLengths = $self->getSequenceLengths(\@distinctSequenceSourceIds, $dbh);

  for(my $i = 0; $i < scalar @sortedGenes; $i++) {
    my $geneSourceId = $sortedGenes[$i]->[1];
    my $sequenceSourceId = $sortedGenes[$i]->[0];
    
    my ($lastGeneEndOrSeqStart, $nextGeneStartOrSeqEnd);
    my ($lastCdsEndOrSeqStart, $nextCdsStartOrSeqEnd);

    if($i == 0 || $sortedGenes[$i-1]->[0] ne $sequenceSourceId) {
      $lastGeneEndOrSeqStart = 1;
      $lastCdsEndOrSeqStart = 1;
    }
    else {
      $lastGeneEndOrSeqStart = $sortedGenes[$i-1]->[3];
      $lastCdsEndOrSeqStart = $sortedGenes[$i-1]->[5];
    }

    if($i == scalar @sortedGenes - 1 || $sortedGenes[$i+1]->[0] ne $sequenceSourceId) {
      $nextGeneStartOrSeqEnd = $sequenceLengths->{$sequenceSourceId};
      $nextCdsStartOrSeqEnd = $sequenceLengths->{$sequenceSourceId};
    }
    else {
      $nextGeneStartOrSeqEnd = $sortedGenes[$i+1]->[2];
      $nextCdsStartOrSeqEnd = $sortedGenes[$i+1]->[4];
    }

    $geneModels->{$geneSourceId}->{last_cds_end_or_seq_start} = $lastCdsEndOrSeqStart;
    $geneModels->{$geneSourceId}->{next_cds_start_or_seq_end} = $nextCdsStartOrSeqEnd;

    $geneModels->{$geneSourceId}->{last_gene_end_or_seq_start} = $lastGeneEndOrSeqStart;
    $geneModels->{$geneSourceId}->{next_gene_start_or_seq_end} = $nextGeneStartOrSeqEnd;
  }

  $self->{_all_models_hash} = $geneModels;
  $self->{_transcript_to_gene_map} = $transcriptToGeneMap;
  $self->{_protein_to_transcript_map} = $proteinToTranscriptMap;

}



sub getSequenceLengths {
  my ($self, $sequenceIds, $dbh) = @_;

  my $sql = "select source_id, length from dots.nasequence where source_id = ?";

  my $sh = getSeqSh();
  unless ($sh) {
    $sh = $dbh->prepare($sql);
    setSeqSh($sh);
  }

  my $rv = {};

  foreach my $sequenceId (@$sequenceIds) {
    $sh->execute($sequenceId);

    while(my ($id, $len) = $sh->fetchrow_array()) {
      $rv->{$id} = $len;
    }
  }

  return $rv;
}



#--------------------------------------------------------------------------------
# Static Methods Start
#--------------------------------------------------------------------------------

sub queryForAgpMap {
  my ($self, $dbh) = @_;

  my %agpMap;

  my $sql = "select sp.virtual_na_sequence_id
                                , p.na_sequence_id as piece_na_sequence_id
                               , decode(sp.strand_orientation, '+', '+1', '-', '-1', '+1') as piece_strand
                               , sp.start_position as piece_start
                               , sp.end_position as piece_end
                               , sp.distance_from_left + 1 as virtual_start_min
                               , sp.distance_from_left + sp.end_position - sp.start_position + 1 as virtual_end_max
                               , p.source_id as piece_source_id
                               , vs.source_id as virtual_source_id
                   from dots.sequencepiece sp
                           , dots.nasequence p
                           ,  dots.nasequence vs
                           ,  sres.ontologyterm o
                   where  sp.piece_na_sequence_id = p.na_sequence_id
                    and p.sequence_ontology_id = o.ontology_term_id
                    and sp.virtual_na_sequence_id = vs.na_sequence_id
                    and o.name not in ('gap')";


  my $sh = getVirtualSeqSh();
  unless ($sh) {
    $sh = $dbh->prepare($sql);
    setVirtualSeqSh($sh);
  }
  $sh->execute();

  while(my $hash = $sh->fetchrow_hashref()) {

    my $virtualSequenceSourceId = $hash->{VIRTUAL_SOURCE_ID};
    my $virtualSequenceNaSequenceId = $hash->{VIRTUAL_NA_SEQUENCE_ID};

    # NEED to map the virtual sequence source_ids to na_sequence_ids
    $self->{_virtual_sequence_ids}->{$virtualSequenceSourceId} = $virtualSequenceNaSequenceId;

    my $ctg = Bio::Location::Simple->new( -seq_id => $hash->{PIECE_SOURCE_ID}, 
                                          -start => $hash->{PIECE_START}, 
                                          -end =>  $hash->{PIECE_END}, 
                                          -strand => '+1' );

    my $ctg_on_chr = Bio::Location::Simple->new( -seq_id =>  $hash->{VIRTUAL_SOURCE_ID}, 
                                                 -start => $hash->{VIRTUAL_START_MIN},
                                                 -end =>  $hash->{VIRTUAL_END_MAX} , 
                                                 -strand => $hash->{PIECE_STRAND} );


    my $agp = Bio::Coordinate::Pair->new( -in  => $ctg, -out => $ctg_on_chr );
    my $pieceSourceId = $hash->{PIECE_SOURCE_ID};
 
    push @{$agpMap{$pieceSourceId}},$agp;
  }

  return \%agpMap;
}


sub mapLocation {
  my ($agpMap, $pieceSourceId, $start, $end, $strand, $wantTopLevel) = @_;

  my $match = Bio::Location::Simple->
    new( -seq_id => $pieceSourceId, -start =>   $start, -end =>  $end, -strand => $strand );

  my $agpArray = $agpMap->{$pieceSourceId};

  if($agpArray) {
    $match->{_sequence_is_piece} = 1;
  }
  else {
    $match->{_sequence_is_piece} = 0;
    return $match;
  }

  return $match unless($wantTopLevel);


  my $agp = &findAgpFromPieceLocation($agpArray, $start, $end, $pieceSourceId);

  my $result = $agp->map($match);

  my $resultIsPiece = $agpMap->{$result->seq_id};

  return $result if($match->seq_id eq $result->seq_id || !$resultIsPiece);

  &mapLocation($agpMap, $result->seq_id, $result->start, $result->end, $result->strand,  $wantTopLevel);
}
 

sub findAgpFromPieceLocation {
  my ($agpArray, $start, $end, $sourceId) = @_;

  foreach my $agp (@$agpArray) {
    my $in = $agp->in();
    if($start >= $in->start() && $end <= $in->end()) {
      return $agp;
    }
  }
  
  die "Could not find a virtual sequence for piece $sourceId w/ start=$start and end=$end";
}

# this is a static function
sub getShortFeatureType {
  my ($featureObject) = @_;

  my $className = ref $featureObject;

  my $featureType;
  if($className =~ /GUS::Community::GeneModelLocations::(.+)/) {
    $featureType = $1;
  }
  else {
    die "Expected Classname to start with GUS::Community::GeneModelLocations but found: $className";
  }
  return $featureType;
}


1;

package GUS::Community::GeneModelLocations::Gene;
use base qw(Bio::SeqFeature::Generic);

1;

package GUS::Community::GeneModelLocations::CDS;
use base qw(Bio::SeqFeature::Generic);
1;

package GUS::Community::GeneModelLocations::Transcript;
use base qw(Bio::SeqFeature::Gene::Transcript);
1;

package GUS::Community::GeneModelLocations::UTR;
use base qw(Bio::SeqFeature::Gene::UTR);
1;

package GUS::Community::GeneModelLocations::Exon;
use base qw(Bio::SeqFeature::Gene::Exon);
1;




1;



