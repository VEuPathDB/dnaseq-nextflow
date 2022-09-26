#!/usr/bin/perl

use strict;
use warnings;
use Test2::V0;
use CalcCoordinates;

# =============== THESE ARRAYS ARE SET WITHIN lib/CalcCoordinates =====================================================================================
# @locationshifts = ([1933,1],[2531,0],[3037,-2],[3254,0],[3433,-2],[8334,-1],[13340,-2],[13848,-4],[19255,-5],[20107,-8]);
# @coordinates = ([250,560,0],[3767,4765,1],[5853,7502,1],[9124,11130,1],[12136,12705,1],[15087,17084,1],[18200,18949,1],[20013,21809,1],[22675,23388,1]);


# ================ VALUES I AM PASSING IN =============================================================================================================
#$shiftFrame $coordinateFrame $locationshiftsLen $oldShift $i


# ================ TESTS ==============================================================================================================================

# Coordinate Prior to any indels 
is( CalcCoordinates::calcCoordinates(0,0,10,0,0), 250 );
is( CalcCoordinates::calcCoordinates(0,0,10,0,1), 560 );

# Coordinate inside indels
is( CalcCoordinates::calcCoordinates(0,3,10,0,0), 9123 );
is( CalcCoordinates::calcCoordinates(0,6,10,0,1), 18945 );

# Coordinate inside indels, proving shiftFrame can be earlier, will shift to appropriate value
is( CalcCoordinates::calcCoordinates(4,3,10,-2,0), 9123 );

# Indel Occuring Between Coordinates as well as Ocurring in the Final Shift Frame 
is( CalcCoordinates::calcCoordinates(0,7,10,0,0), 20008 );
is( CalcCoordinates::calcCoordinates(0,7,10,0,1), 21801 );

# Coordinate Greater than last indel location
is( CalcCoordinates::calcCoordinates(0,8,10,0,0), 22667 );
is( CalcCoordinates::calcCoordinates(0,8,10,0,1), 23380 );

done_testing();
