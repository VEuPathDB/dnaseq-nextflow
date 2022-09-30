#!/usr/bin/perl

use strict;
use warnings;
use Test2::V0;
use calcTranscriptLocation;

# =============== THESE ARRAYS ARE SET WITHIN lib/calcTranscriptLocation =====================================================================================

# @locationshifts = ([4000,-2],[4334,-1],[9340,-2],[15848,-4],[18255,-5],[20107,-8]);
# @coordinatesNonCodingLen = ([3767,4765,3766],[5853,7502,4853],[9124,11130,6474],[12136,12705,7479],[15087,17084,9860],[18200,18949,10975],[20013,21809,12038]);

# ================ VALUES I AM PASSING IN =============================================================================================================

# $snpLocation $shiftFrame $locationsShiftedNonCodingFrame

# ================ TESTS ==============================================================================================================================

is( calcTranscriptLocation::calcTranscriptLocation(4000,0,0), 234 );

is( calcTranscriptLocation::calcTranscriptLocation(4334,1,0), 566 );

is( calcTranscriptLocation::calcTranscriptLocation(9340,2,2), 2865 );

is( calcTranscriptLocation::calcTranscriptLocation(15848,3,4), 5986 );

is( calcTranscriptLocation::calcTranscriptLocation(18225,4,5), 7246 );

is( calcTranscriptLocation::calcTranscriptLocation(20107,5,6), 8064 );
    
done_testing();
