#!/usr/bin/perl

use strict;
use warnings;
use Test2::V0;
use shifted;

# =============== THESE ARRAYS ARE SET WITHIN lib/CalcCoordinates =====================================================================================
# @locationshifts = ([1933,1],[2531,0],[3037,-2],[3254,0],[3433,-2],[8334,-1],[13340,-2],[13848,-4],[19255,-5],[20107,-8]);
# @coordinates = ([250,560,0],[3767,4765,1],[5853,7502,1],[9124,11130,1],[12136,12705,1],[15087,17084,1],[18200,18949,1],[20013,21809,1],[22675,23388,1]);


# ================ VALUES I AM PASSING IN =============================================================================================================
# $locationShiftedFrame $value


# ================ TESTS ==============================================================================================================================

# Coordinate Prior to any indels 
is( shifted::locationShifted(0,0), 250 );
is( shifted::locationShifted(0,1), 560 );

# Coordinate inside indels
is( shifted::locationShifted(1,0), 3765 );
is( shifted::locationShifted(1,1), 4763 );
is( shifted::locationShifted(3,0), 9123 );
is( shifted::locationShifted(4,1), 12704 );
is( shifted::locationShifted(6,1), 18945 );

# Indel Occuring Between Positions, Last Indel Occuring Between Positions
is( shifted::locationShifted(7,0), 20008 );
is( shifted::locationShifted(7,1), 21801 );

# Coordinate Greater than last indel location
is( shifted::locationShifted(8,0), 22667 );
is( shifted::locationShifted(8,1), 23380 );


# Reverse Value Tests
is( shifted::locationShifted(0,2), 0 );
is( shifted::locationShifted(1,2), 1 );
is( shifted::locationShifted(3,2), 1 );
is( shifted::locationShifted(5,2), 1 );
is( shifted::locationShifted(7,2), 1 );


done_testing();
