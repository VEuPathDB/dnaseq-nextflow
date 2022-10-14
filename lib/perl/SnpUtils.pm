package VEuPath::SnpUtils;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(sampleCacheFileColumnNames snpFileColumnNames isSameSNP );
use strict;


sub sampleCacheFileColumnNames {
  my @columnNames = 
      ('strain',
       'location',
       'sequence_source_id',
       'reference',
       'base',
       'percent',
       'matches_reference',
       'quality',
       'coverage',
      );

  return wantarray ? @columnNames : \@columnNames;
}


sub snpFileColumnNames {
  my @columnNames = 
      ('position_in_protein',
       'is_coding',
       'source_id',
       'minor_allele',
       'minor_product',
       'reference_strain',
       'has_stop_codon',
       'location' => '52768',
       'positions_in_cds',
       'gene_na_feature_id',
       'minor_allele_count',
       'position_in_cds',
       'distinct_allele_count',
       'major_allele' => 'A',
       'total_allele_count',
       'major_product',
       'reference_aa',
       'na_sequence_id',
       'reference_aa_full',
       'major_allele_count',
       'reference_na' => 'A',
       'distinct_strain_count',
       'has_nonsynonymous_allele',
       'positions_in_protein',
       'external_database_release_id'
      );

   return wantarray ? @columnNames : \@columnNames;
}


sub isSameSNP {
  my ($a, $b) = @_;

    my $strain = $a->[0];
    my $peekstrain = $b->[0];

    my $location = $a->[1];
    my $peekLocation = $b->[1];
  
  return $strain eq $peekstrain && $location == $peekLocation;

}
