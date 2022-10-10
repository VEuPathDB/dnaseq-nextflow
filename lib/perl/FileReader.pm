package VEuPath::FileReader;

use strict;

sub getFile {$_[0]->{_file}}
sub setFile {$_[0]->{_file} = $_[1]}

sub getFh {$_[0]->{_fh}}
sub setFh {$_[0]->{_fh} = $_[1]}

sub getFilters {$_[0]->{_filters}}
sub setFilters {$_[0]->{_filters} = $_[1]}

sub getPeekLine {$_[0]->{_peek_line}}
sub setPeekLine {$_[0]->{_peek_line} = $_[1]}

sub getPeekLineAsArray {$_[0]->{_peek_line_as_array}}
sub setPeekLineAsArray {$_[0]->{_peek_line_as_array} = $_[1]}

sub getDelimiter {$_[0]->{_delimiter}}
sub setDelimiter {$_[0]->{_delimiter} = $_[1]}

sub getDictionaryNames {$_[0]->{_dictionary_names}}
sub setDictionaryNames {$_[0]->{_dictionary_names} = $_[1]}

sub hasDictionary { 
  my ($self) = @_;

  if($self->getDictionaryNames()) {
    return 1;
  }

  return 0;
}

sub getPeek {
  my ($self) = @_;

  return wantarray ? @{$self->getPeekLineAsArray()} : $self->getPeekLine();
}

sub new {
  my ($class, $file, $filters, $delimiter) = @_;

  my $self = bless {}, $class;

  $self->setFile($file);
  $self->setFilters($filters);

  if($delimiter) {
    $self->setDelimiter($delimiter);
  }
  else {
    $self->setDelimiter(qr//);
  }

  my ($fh);
  open($fh, $file) or die "Cannot open file $file for reading: $!";

  $self->setFh($fh);

  $self->processNext();

  return $self;
}

# this skips empty lines; can override to skip lines as needed (skip based on filters for example)
sub skipLine {
  my ($self, $line, $a, $fh) = @_;

  return !$line;
}

# used for collecting like rows into an array of rows (each row should be made into a dictionary / hash)
#  Default is false;  Subclasses can override as needed
sub isSameGroup {
  return 0;
}

sub readNextLine {
  my ($self, $fh) = @_;

  # handle empty lines or whatever
  while(!eof($fh)) {
    my $line = readline($fh);
    chomp($line);

    my $delimiter = $self->getDelimiter();
    my @a = split(/$delimiter/, $line);

    next if($self->skipLine($line, \@a, $fh));

    return($line, \@a);
  }
  return(undef, []);
}

sub hasNext {
  my ($self) = @_;

  if($self->getPeek()) {
    return 1;
  }

  $self->closeFileHandle();
  return 0;
}


sub closeFileHandle {
  my ($self) = @_;

  my $fh = $self->getFh();

  close $fh;
}


# this is the one which will be called by users
sub nextLine {
  my ($self) = @_;

  return undef unless($self->hasNext());

  my $line = $self->getPeek();
  my @lineAsArray = $self->getPeek();

  $self->processNext();

  return wantarray ? @lineAsArray : $line;
}

sub processNext {
  my ($self) = @_;

  my $fh = $self->getFh();
  my ($line, $lineAsArray) = $self->readNextLine($fh);

  $self->setPeekLineAsArray($lineAsArray);
  $self->setPeekLine($line);
}

sub readNextGroupOfLines {
  my ($self) = @_;

    my @rv;

  my $isSameGroup = 1;

  while($isSameGroup) {
    last unless($self->hasNext());

    my @a = $self->nextLine();
    my @b = $self->getPeek();

    unless($self->isSameGroup(\@a, \@b)) {
      $isSameGroup = 0;
    }

    if($self->hasDictionary()) {
      my $lineAsDictionary = $self->makeDictionary(\@a);
      push @rv, $lineAsDictionary;
    }
    else {
      push @rv, \@a;
    }

  }
  return \@rv;

}


sub makeDictionary {
  my ($self, $lineAsArray) = @_;

  my $dictionaryNames = $self->getDictionaryNames();

  unless($dictionaryNames) {
    die "Cannot call makeDictionary unless you have first setDictionaryNames(\@names)";
  }

  my %dict;

  for(my $i = 0; $i < scalar @$dictionaryNames ; $i++) {
    my $key = $dictionaryNames->[$i];
    my $value = $lineAsArray->[$i];

    $dict{$key} = $value;
  }

  return \%dict;
}


1;

