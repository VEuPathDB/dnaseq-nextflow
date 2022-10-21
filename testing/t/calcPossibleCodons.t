#!/usr/bin/perl

use strict;
use warnings;
use Test2::V0;
use VEuPath::calcPossibleCodons;

# ================ TESTS ==============================================================================================================================

# This function does not return the codons in a string as is shown here. It sets the array of codons to the {product} field of the variation hash. There is additional functionally added to the test calcPossibleCodons module to return a string for ease of testing

# No ambiguities
is( VEuPath::calcPossibleCodons::calcPossibleCodons('AAA'), 'AAA' );
# Ambiguity with two options in last position
is( VEuPath::calcPossibleCodons::calcPossibleCodons('AAR'), 'AAAAAG' );
# Ambiguity with two options in middle position
is( VEuPath::calcPossibleCodons::calcPossibleCodons('AYA'), 'ACAATA' );
# Ambiguity with two options in first position
is( VEuPath::calcPossibleCodons::calcPossibleCodons('KAA'), 'GAATAA' );
# Ambiguity with three options in last position
is( VEuPath::calcPossibleCodons::calcPossibleCodons('AAB'), 'AAGAATAAC' );
# Ambiguity with three options in middle position
is( VEuPath::calcPossibleCodons::calcPossibleCodons('ADA'), 'AGAAAAATA' );
# Ambiguity with three options in first position
is( VEuPath::calcPossibleCodons::calcPossibleCodons('HAA'), 'AAACAATAA' );
# Ambiguity with four options in last position
is( VEuPath::calcPossibleCodons::calcPossibleCodons('AAN'), 'AAAAAGAACAAT' );
# Ambiguity with four options in middle position
is( VEuPath::calcPossibleCodons::calcPossibleCodons('ANA'), 'AAAAGAACAATA' );
# Ambiguity with four options in first position
is( VEuPath::calcPossibleCodons::calcPossibleCodons('NAA'), 'AAAGAACAATAA' );

# 2 Ambiguities with two options
is( VEuPath::calcPossibleCodons::calcPossibleCodons('RAR'), 'AAAAAGGAAGAG' );
# 2 Ambiguities with three options
is( VEuPath::calcPossibleCodons::calcPossibleCodons('ABB'), 'AGGAGTAGCATGATTATCACGACTACC' );
# 2 Ambiguities with four options
is( VEuPath::calcPossibleCodons::calcPossibleCodons('NNA'), 'AAAAGAACAATAGAAGGAGCAGTACAACGACCACTATAATGATCATTA' );

# amBIGuity
is( VEuPath::calcPossibleCodons::calcPossibleCodons('NNN'), 'AAAAAGAACAATAGAAGGAGCAGTACAACGACCACTATAATGATCATTGAAGAGGACGATGGAGGGGGCGGTGCAGCGGCCGCTGTAGTGGTCGTTCAACAGCACCATCGACGGCGCCGTCCACCGCCCCCTCTACTGCTCCTTTAATAGTACTATTGATGGTGCTGTTCATCGTCCTCTTTATTGTTCTTT' );

done_testing();
