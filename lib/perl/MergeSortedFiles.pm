package VEuPath::MergeSortedFiles;

use base 'VEuPath::FileReader';
use strict;

#--------------------------------------------------------------------------------

## SUBCLASSES SHOULD override this w/ something interesting
## NOTE:  Pay special attention to "locale" and LC_COLLATE if you are comparing lines from unix sorted files;  Perl's default is not to 'use locale'.  Ie. your subclass must know how the files were sorted originally
sub wantFirstLine {
  my ($self) = @_;

  my $firstLine = $self->getFirstLine();
  my $secondLine = $self->getSecondLine();

  return $firstLine le $secondLine
}

#--------------------------------------------------------------------------------

sub getFirstLine {$_[0]->{_first_line}}
sub setFirstLine {$_[0]->{_first_line} = $_[1]}

sub getSecondLine {$_[0]->{_second_line}}
sub setSecondLine {$_[0]->{_second_line} = $_[1]}

sub getFirstLineAsArray {$_[0]->{_first_line_as_array}}
sub setFirstLineAsArray {$_[0]->{_first_line_as_array} = $_[1]}

sub getSecondLineAsArray {$_[0]->{_second_line_as_array}}
sub setSecondLineAsArray {$_[0]->{_second_line_as_array} = $_[1]}

sub getFirstFh {$_[0]->{_first_fh}}
sub setFirstFh {$_[0]->{_first_fh} = $_[1]}

sub getSecondFh {$_[0]->{_second_fh}}
sub setSecondFh {$_[0]->{_second_fh} = $_[1]}

#--------------------------------------------------------------------------------

sub readingFile1Fh {
  my ($self, $fh) = @_;

  return $self->getFh() eq $fh;
}

sub readingFile2Fh {
  my ($self, $fh) = @_;

  return $self->getFh() ne $fh;
}

#--------------------------------------------------------------------------------
# @OVERRIDE
sub new {
  my ($class, $file1, $file2, $filters, $delimiter) = @_;

  my $self = bless {}, $class;

  if($delimiter) {
    $self->setDelimiter($delimiter);
  }
  else {
    $self->setDelimiter(qr//);
  }

  $self->setFile($file1);
  $self->setFilters($filters);

  my ($file1Fh, $file2Fh);
  open($file1Fh, $file1) or die "Cannot open file $file1 for reading: $!";
  open($file2Fh, $file2) or die "Cannot open file $file2 for reading: $!";

  $self->setFh($file1Fh);

  my ($file1Line, $file1LineAsArray) = $self->readNextLine($file1Fh);
  my ($file2Line, $file2LineAsArray) = $self->readNextLine($file2Fh);

  unless($file1Line || $file2Line) {
    print STDERR "WARN:  Neither of the 2 input files contain any rows.\n";
  }

  $self->setFirstLine($file1Line);
  $self->setFirstLineAsArray($file1LineAsArray);
  $self->setFirstFh($file1Fh);

  $self->setSecondLine($file2Line);
  $self->setSecondLineAsArray($file2LineAsArray);
  $self->setSecondFh($file2Fh);

  $self->processNext();

  return $self;
}


sub merge {
  my ($self) = @_;

  my $firstFh = $self->getFirstFh();
  my $secondFh = $self->getSecondFh();

  my $firstLine = $self->getFirstLine();
  my $secondLine = $self->getSecondLine();

  my $firstLineAsArray = $self->getFirstLineAsArray();
  my $secondLineAsArray = $self->getSecondLineAsArray();

  my $rv;
  my $rvAsArray;

  if($self->wantFirstLine()) {
    $rv = $firstLine;
    $rvAsArray = $firstLineAsArray;
  }

  else {
    $self->setFirstLine($secondLine);
    $self->setFirstLineAsArray($secondLineAsArray);
    $self->setFirstFh($secondFh);

    $self->setSecondLine($firstLine);
    $self->setSecondLineAsArray($firstLineAsArray);
    $self->setSecondFh($firstFh);

    $rv = $secondLine;
    $rvAsArray = $secondLineAsArray;
  }

  
  my $fh = $self->getFirstFh();
  ($firstLine, $firstLineAsArray) = $self->readNextLine($fh);

  $self->setFirstLine($firstLine);
  $self->setFirstLineAsArray($firstLineAsArray);

  return($rv, $rvAsArray);
}


# @OVERRIDE
sub closeFileHandle {
  my ($self) = @_;
  $self->closeFileHandles(); 
}


sub closeFileHandles {
  my ($self) = @_;

  my $fh1 = $self->getFirstFh();
  my $fh2 =  $self->getSecondFh();

  close $fh1;
  close $fh2;
}


# @OVERRIDE
sub processNext {
  my ($self) = @_;

  my $firstFh = $self->getFirstFh();
  my $secondFh = $self->getSecondFh();

  my $firstLine = $self->getFirstLine();
  my $secondLine = $self->getSecondLine();

  my $nextLine;
  my $nextLineAsArray;

  if(!$firstLine) {
    $nextLine = $secondLine;
    $nextLineAsArray = $self->getSecondLineAsArray();

    my ($l, $laa) = $self->readNextLine($secondFh);
    $self->setSecondLine($l);
    $self->setSecondLineAsArray($laa);
  }
  elsif(!$secondLine) {
    $nextLine = $firstLine;
    $nextLineAsArray = $self->getFirstLineAsArray();

    my ($l, $laa) = $self->readNextLine($firstFh);
    $self->setFirstLine($l);
    $self->setFirstLineAsArray($laa);
  }
  else {
    ($nextLine, $nextLineAsArray) = $self->merge();
  }

  $self->setPeekLineAsArray($nextLineAsArray);
  $self->setPeekLine($nextLine);
}


1;
