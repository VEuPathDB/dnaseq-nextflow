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
       'is_coding',
       'position_in_cds',
       'position_in_codon',
       'downstream_of_frameshift',
       'transcript',
       'product',
       'reference_codon',
       'codon',
       'has_nonsynonomous',
       'shifted_location',
       'cds_number',
       'current_shift'
       );

  return wantarray ? @columnNames : \@columnNames;
}


sub snpFileColumnNames {
    my @columnNames =  (
	     "location",
	     "transcript_id",
             "source_id",
	     "reference_strain",
             "reference_na",
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
	     "ref_codon"
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
	     "position_in_protein",
	     "transcript",
             "count",
             "product",
             "ref_location_cds",
	     "ref_location_protein",
	     "downstream_of_frameshift"
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
