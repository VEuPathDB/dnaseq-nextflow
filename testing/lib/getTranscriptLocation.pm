#!/usr/bin/perl

package getTranscriptLocation;
use strict;

my @locationshifts = ([4000,-2],[4334,-1],[5000,0],[9340,-2],[12000,-1],[15848,-4],[18255,-5],[20000,-4],[20107,-8]);
my @coordinatesNonCodingLen = ([3767,4765,3766],[5853,7502,4853],[9124,11130,6474],[12136,12705,7479],[15087,17084,9860],[18200,18949,10975],[20013,21809,12038]);

sub getTranscriptLocation {
    my $snpLocation = $_[0];
    my $coordinatesNonCodingFrame = $_[1];
    my $coordinatesNonCodingFrameLimit = $_[2];
    my $shiftFrame = $_[3];
    my $lastCodingShiftFrame = $_[4];
    my $transcriptLocation;
    if ($coordinatesNonCodingFrame == $coordinatesNonCodingFrameLimit) {
        if ($snpLocation >= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][0] && $snpLocation <= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][1]) {
            ($transcriptLocation, $lastCodingShiftFrame) = &calcTranscriptLocation($snpLocation, $shiftFrame, $coordinatesNonCodingFrame, $lastCodingShiftFrame);
        }
        else {
           $transcriptLocation = "NC";
        }
    }
    elsif ($snpLocation >= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][0] && $snpLocation <= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][1]) {
         ($transcriptLocation, $lastCodingShiftFrame) = &calcTranscriptLocation($snpLocation, $shiftFrame, $coordinatesNonCodingFrame, $lastCodingShiftFrame);
    }
    elsif ($snpLocation > $coordinatesNonCodingLen[$coordinatesNonCodingFrame][1] || $coordinatesNonCodingFrame == $coordinatesNonCodingFrameLimit) {
        until ($snpLocation <= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][1] || $coordinatesNonCodingFrame == $coordinatesNonCodingFrameLimit) {
            $coordinatesNonCodingFrame++;
        }
        if ($coordinatesNonCodingFrame == $coordinatesNonCodingFrameLimit) {
            if ($snpLocation >= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][0] && $snpLocation <= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][1]) {
               ($transcriptLocation, $lastCodingShiftFrame) = &calcTranscriptLocation($snpLocation, $shiftFrame, $coordinatesNonCodingFrame, $lastCodingShiftFrame);
            }
            else {
                $transcriptLocation = "NC";
            }
        }
        elsif ($snpLocation >= $coordinatesNonCodingLen[$coordinatesNonCodingFrame][0]) {
             ($transcriptLocation, $lastCodingShiftFrame) = &calcTranscriptLocation($snpLocation, $shiftFrame, $coordinatesNonCodingFrame, $lastCodingShiftFrame);
        }
        else {
            $transcriptLocation = "NC";
        }
    }
    else {
        $transcriptLocation = "NC";
    }
    return $transcriptLocation;
}

sub calcTranscriptLocation {
    my ($snpLocation, $shiftFrame, $coordinatesNonCodingFrame, $lastCodingShiftFrame) = @_;
    my $transcriptLocation;
    if ($shiftFrame == 0) {
        $transcriptLocation = $snpLocation - $coordinatesNonCodingLen[$coordinatesNonCodingFrame][2];
	$lastCodingShiftFrame = $shiftFrame;
    }
    else {
        $transcriptLocation = $snpLocation - $coordinatesNonCodingLen[$coordinatesNonCodingFrame][2] + $locationshifts[$lastCodingShiftFrame][1];
	$lastCodingShiftFrame = $shiftFrame;
    }
    return ($transcriptLocation, $lastCodingShiftFrame);
}

1;
