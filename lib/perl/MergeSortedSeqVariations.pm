package VEuPath::MergeSortedSeqVariations;

use base 'VEuPath::MergeSortedFiles';
use strict;
use locale;  
use VEuPath::SnpUtils;


sub getSkipCount {
  my ($self) = @_;

  return $self->{_skip_count};
}


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $columnNames = &fileColumnNames();
  $self->setDictionaryNames($columnNames);

  return $self;
}


# @OVERRIDE
sub wantFirstLine {
  my ($self) = @_;

  my @a = @{$self->getFirstLineAsArray()};
  my @b = @{$self->getSecondLineAsArray()};

  return $a[0] lt $b[0] || ($a[0] eq $b[0] && $a[1] <= $b[1]);
}


# @OVERRIDE
sub skipLine {
  my ($self, $line, $lineAsArray, $fh) = @_;

  return 1 unless($line);
  return 0 if($self->readingFile1Fh($fh));

  my $filters = $self->getFilters();

  foreach(@$filters) {
    if($lineAsArray->[0] eq $_) {
      $self->{_skip_count}++;
      return 1;
    }
  }
  return 0;
}


sub nextSNP {
  my ($self) = @_;
  $self->readNextGroupOfLines();
}


# @OVERRIDE
sub isSameGroup {
  my ($self, $a, $b) = @_;
  return &isSameSNP($a, $b);
}

1;
