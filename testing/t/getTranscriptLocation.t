#!/usr/bin/perl

use getTranscriptLocation;
use strict;
use warnings;
use Test2::V0;

# =============== THESE ARRAYS ARE SET WITHIN lib/calcTranscriptLocation =====================================================================================
# @locationshifts = ([4000,-2],[4334,-1],[5000,0],[9340,-2],[12000,-1],[15848,-4],[18255,-5],[20000,-4],[20107,-8]);
# @coordinatesNonCodingLen = ([3767,4765,3766],[5853,7502,4853],[9124,11130,6474],[12136,12705,7479],[15087,17084,9860],[18200,18949,10975],[20013,21809,12038]);

# ================ VALUES I AM PASSING IN =============================================================================================================
# $snpLocation $coordinatesNonCodingFrame $coordinatesNonCodingFrameLimit $shiftFrame $lastCodingShiftFrame

# ================ TESTS ==============================================================================================================================

is( getTranscriptLocation::getTranscriptLocation(4000,0,6,0,0), 234 );

is( getTranscriptLocation::getTranscriptLocation(4334,0,6,1,0), 566 );

is( getTranscriptLocation::getTranscriptLocation(5000,0,6,2,1), "NC" );

is( getTranscriptLocation::getTranscriptLocation(9340,0,6,3,1), 2865 ); 

is( getTranscriptLocation::getTranscriptLocation(12000,0,6,4,2,3), "NC" );

is( getTranscriptLocation::getTranscriptLocation(18255,0,6,5,5), 7276 );
    
done_testing();
