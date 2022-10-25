package VEuPath::SnpUtils;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(sampleCacheFileColumnNames snpFileColumnNames isSameSNP alleleFileColumnNames productFileColumnNames);
use strict;

sub sampleCacheFileColumnNames {
  my @columnNames = 
       ('sequence_source_id',
       'location',
       'strain',
       'reference',
       'base',
       'coverage',
       'percent',
       'quality',
       'pvalue',
       'snp_source_id',
       'external_database_release_id',
       'matches_reference',
       'product',
       'position_in_cds',
       'position_in_protein',
       'na_sequence_id',
       'ref_na_sequence_id',
       'snp_external_database_release_id',
       'protocol_app_node_id',
       'has_nonsynonomous',
       'is_coding'	
      );

  return wantarray ? @columnNames : \@columnNames;
}


sub snpFileColumnNames {
    my @columnNames =  (
	     "location",
	     "gene_na_feature_id",
             "source_id",
	     "na_sequence_id",
             "reference_strain",
             "reference_na",
             "reference_aa",
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
             "has_coding_mutation",
             "total_allele_count",
             "has_stop_codon",
             "ref_codon",
             );

   return wantarray ? @columnNames : \@columnNames;
}


sub alleleFileColumnNames {
    my @columnNames =  (
	"allele",
	"distinct_strain_count",
	"allele_count",
	"average_coverage",
	"average_read_percent"
	);
   return wantarray ? @columnNames : \@columnNames;
}


sub productFileColumnNames {
    my @columnNames =  (
	     "codon",
	     "position_in_codon",
	     "transcript",
             "count",
             "product",
             "ref_location_cds",
             "ref_location_protein"
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
