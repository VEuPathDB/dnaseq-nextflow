#!/usr/bin/perl

package VEuPath::calcPossibleCodons;
use strict;
use Set::CrossProduct;
use Data::Dumper;

#===================================== GENERATING Possible Codons =========================================================================

sub calcPossibleCodons {
    my ($codon) = @_;
    my @codonArray=split(//, $codon);
    my %translate = (A => ['A'],
                     G => ['G'],
                     C => ['C'],
                     T => ['T'],
                     R => ['A','G'],
                     Y => ['C','T'],
                     K => ['G','T'],
                     M => ['A','C'],
                     S => ['G','C'],
                     W => ['A','T'],
                     B => ['G','T','C'],
                     D => ['G','A','T'],
                     H => ['A','C','T'],
                     V => ['G','C','A'],
                     N => ['A','G','C','T']
                     );
    my @expanded = map { $translate{$_} } @codonArray;
    my $iterator = Set::CrossProduct->new(\@expanded);
    my $codonList;
    foreach my $codon ($iterator->combinations) {
	my $string = join(",", @$codon);
	$string = $string =~ s/,//gr;
	push @{ $codonList }, $string;
    }
    my $count = scalar @{ $codonList };
    my $codonString;
    foreach my $codon (0..$count) {
        $codonString = $codonString ."$codonList->[$codon]";
    }
    return $codonString;
}
   
1;
