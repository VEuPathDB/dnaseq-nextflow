#!/usr/bin/perl

package VEuPath::calcTranscriptLocation;
use strict;

my @locationshifts = ([4000,-2],[4334,-1],[9340,-2],[15848,-4],[18255,-5],[20107,-8]);
my @coordinatesNonCodingLen = ([3767,4765,3766],[5853,7502,4853],[9124,11130,6474],[12136,12705,7479],[15087,17084,9860],[18200,18949,10975],[20013,21809,12038]);

sub calcTranscriptLocation {
    my $snpLocation = $_[0];
    my $shiftFrame = $_[1];
    my $coordinatesNonCodingFrame = $_[2];
    my $lastCodingShiftFrame = $_[3];
    my $transcriptLocation;
    if ($shiftFrame == 0) {
        $transcriptLocation = $snpLocation - $coordinatesNonCodingLen[$coordinatesNonCodingFrame][2];
	$lastCodingShiftFrame = $shiftFrame;
    }
    else {
        $transcriptLocation = $snpLocation - $coordinatesNonCodingLen[$coordinatesNonCodingFrame][2] + $locationshifts[$lastCodingShiftFrame][1];
    }
    return $transcriptLocation;
}

1;
