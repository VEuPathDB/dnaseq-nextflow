#!/usr/bin/perl

use strict;
use warnings;
use Test2::V0;
use VEuPath::shifted;

# =============== THESE ARRAYS ARE SET WITHIN lib/CalcCoordinates =====================================================================================
# @locationshifts = ([1933,1],[2531,-1],[2800,1],[2900,-1],[3037,-2],[3254,0],[3433,-2],[8334,-1],[13340,-2],[13848,-4],[19255,-5],[20107,-8]);
# @coordinates = ([250,560,0],[1933,2000,1],[2531,2700,0],[2750,2800,0],[2850,2900,1],[3767,4765,1],[5853,7502,1],[9124,11130,1],[12136,12705,1],[15087,17084,1],[18200,18949,1],[20013,21809,1],[22675,23388,1]);


# ================ VALUES I AM PASSING IN =============================================================================================================
# $locationShiftedFrame $value


# ================ TESTS ==============================================================================================================================

# Coordinate Prior to any indels 
is( VEuPath::shifted::locationShifted(0,0), 250 );
is( VEuPath::shifted::locationShifted(0,1), 560 );

# CDS start equal to SNP location with a positive shift
is( VEuPath::shifted::locationShifted(1,0), 1933 );

# CDS start equal to SNP location with a negative shift
is( VEuPath::shifted::locationShifted(2,0), 2532 );

# CDS end equal to SNP location with a positive shift
is( VEuPath::shifted::locationShifted(3,1), 2801 );

# CDS end equal to SNP location with a negative shift
is( VEuPath::shifted::locationShifted(4,1), 2901 );

# Coordinate inside indels
is( VEuPath::shifted::locationShifted(7,0), 9123 );
is( VEuPath::shifted::locationShifted(10,1), 18945 );

# Coordinate inside indels, proving shiftFrame can be earlier, will shift to appropriate value
is( VEuPath::shifted::locationShifted(7,0), 9123 );

# Indel Occuring Between Coordinates as well as Ocurring in the Final Shift Frame 
is( VEuPath::shifted::locationShifted(11,0), 20008 );
is( VEuPath::shifted::locationShifted(11,1), 21801 );

# Coordinate Greater than last indel location
is( VEuPath::shifted::locationShifted(12,0), 22667 );
is( VEuPath::shifted::locationShifted(12,1), 23380 );

done_testing();
