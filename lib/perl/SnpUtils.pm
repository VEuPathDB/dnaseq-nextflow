package VEuPath::SnpUtils;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(fileColumnNames isSameSNP);
use strict;

sub fileColumnNames {
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


sub isSameSNP {
  my ($a, $b) = @_;

    my $strain = $a->[0];
    my $peekstrain = $b->[0];

    my $location = $a->[1];
    my $peekLocation = $b->[1];
  
  return $strain eq $peekstrain && $location == $peekLocation;

}
