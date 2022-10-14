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
    my @columnNames =  ("gene_na_feature_id",
             "source_id",
             "na_sequence_id",
             "location",
             "reference_strain",
             "reference_na",
             "reference_aa",
             "ref_position_in_cds",
             "ref_position_in_protein",
             "external_database_release_id",
             "has_nonsynonymous_allele",
             "major_allele",
             "minor_allele",
             "major_allele_count",
             "minor_allele_count",
             "major_product",
             "minor_product",
             "distinct_strain_count",
             "distinct_allele_count",
             "is_coding",
             "reference_aa_full",
             "total_allele_count",
             "has_stop_codon",
             "transcript",
             "product",
             "codon",
             "position_in_codon",
             "ref_codon",
             "snp_position_in_cds",
             "snp_position_in_protein",
             "shifted_location",
             "has_coding_mutation",
             "is_downstream_of_frameshift"
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
