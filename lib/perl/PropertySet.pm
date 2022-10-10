package VEuPath::PropertySet;

use strict;
use Carp;

# smart parsing of properties file.
# file is of the format:
# name=value
# '#' at beginning of line is a comment
#
# propDeclaration is a reference to an array of (name, default, comment) for
# each property.  if default is not provided, property is required.
#
# if $relax is set, then allow properties in the file that are not in the declaration.
sub new {
    my ($class, $propsFile, $propsDeclaration, $relax) = @_;

    my $self = {};
    bless($self, $class);

    $self->{props} = {};
    $self->{dynamicDefaults} = {};
    $self->{decl} = $propsDeclaration;
    $self->{file} = $propsFile;

    my $fatalError;

    foreach my $decl (@$propsDeclaration) {
        my $name = $decl->[0];
        my $value = $decl->[1];
        $self->{props}->{$name} = $value? $value : "REQD_PROP";
	if (defined($decl->[3])) {
	  $self->{dynamicDefaults}->{$name} = $decl->[3];
	}
      }

    if ($propsFile) {
#      print STDERR "Reading properties from $propsFile\n";

      open(F, $propsFile) || die "Can't open property file $propsFile";

      my %duplicateCheck;
      while (<F>) {
        chomp;
        s/\s+$//;
        next if (!$_ || /^\s*#/);
	if (! /(\S+?)\s*=\s*(.*)\s*/) {
	  print STDERR "In property file '$propsFile', the following line is malformed: '$_'\n";
	  $fatalError = 1;
	}
        my $key = $1;
        my $value = $2;

        if ($duplicateCheck{$key}) {
          print STDERR "Property name '$key' is duplicated in property file '$propsFile'\n";
          $fatalError = 1;
	}
        $duplicateCheck{$key} = 1;

        if (!$relax && !$self->{props}->{$key}) {
          print STDERR "Invalid property name '$key' in property file '$propsFile'\n";
          $fatalError = 1;
	}

        # allow value to include $ENV{} expressions to include environment vars
        $value =~ s/\$ENV\{"?'?(\w+)"?'?\}/$ENV{$1}/g;

        $self->{props}->{$key} = $value;
      }
      close(F);
    }

    foreach my $name (keys %{$self->{props}}) {

      # if we don't have a value for a required key, try to find a
      # dynamic default
      if ($self->{props}->{$name} eq "REQD_PROP") {
	my $dynamicDefaultKey = $self->{dynamicDefaults}->{$name};
	if ($dynamicDefaultKey && defined($self->{props}->{$dynamicDefaultKey})) {
	  $self->{props}->{$name} = $self->{props}->{$dynamicDefaultKey};
	}
      }

      # if still don't then error
      if ($self->{props}->{$name} eq "REQD_PROP")  {
	print STDERR "Required property '$name' must be specified in property file '$propsFile'\n";
	$fatalError = 1;
      }
    }

    die "Fatal PropertySet error(s)" if $fatalError;
    return $self;
}

sub getProp {
    my ($self, $name) = @_;

    my $value = $self->{props}->{$name};
    confess "Cannot find value for property '$name' in file '$self->{file}'\n" unless ($value ne "");
    return $value;
}

sub getPropRelaxed {
    my ($self, $name) = @_;

    return $self->{props}->{$name};
}


sub setProp {
    my ($self, $name, $value) = @_;
    $self->{props}->{$name} = $value;
}

sub addProperty {
    my ($self, $name, $value) = @_;
    $self->{props}->{$name} = $value;
}

sub toString {
  my $self = shift;
  my $ret = "property = value, help\n----------------------\n";
  foreach my $p (sort keys %{$self->{props}}){
    $ret .= "$p = $self->{props}->{$p}".($self->{help}->{$p} ? ", \"$self->{help}->{$p}\"\n" : "\n");
  }
  return $ret;
}


1;
