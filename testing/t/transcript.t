#!/usr/bin/perl

use transcript;
use strict;
use warnings;
use Test2::V0;

# =============== THESE ARRAYS ARE SET WITHIN lib/calcTranscriptLocation =====================================================================================
# @locationshifts = ([4000,-2],[4334,-1],[5000,0],[9340,-2],[12000,-1],[15848,-4],[18255,-5],[20000,-4],[20107,-8]);
# @coordinatesNonCodingLen = ([3767,4765,3766],[5853,7502,4853],[9124,11130,6474],[12136,12705,7479],[15087,17084,9860],[18200,18949,10975],[20013,21809,12038]);

# ================ VALUES I AM PASSING IN =============================================================================================================
# $frame

# ================ TESTS ==============================================================================================================================

is( transcript::transcript(0), 234 );

is( transcript::transcript(1), 566 );

is( transcript::transcript(2), "NC" );

is( transcript::transcript(3), 2865 ); 

is( transcript::transcript(4), "NC" );

is( transcript::transcript(5), 5986 );
    
done_testing();
