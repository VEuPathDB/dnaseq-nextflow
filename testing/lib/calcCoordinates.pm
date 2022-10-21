#!/usr/bin/perl

package VEuPath::calcCoordinates;
use strict;

my @locationshifts = ([1933,1],[2531,-1],[2800,1],[2900,-1],[3037,-2],[3254,0],[3433,-2],[8334,-1],[13340,-2],[13848,-4],[19255,-5],[20107,-8]);
my @coordinates = ([250,560,0],[1933,2000,1],[2531,2700,0],[2750,2800,0],[2850,2900,1],[3767,4765,1],[5853,7502,1],[9124,11130,1],[12136,12705,1],[15087,17084,1],[18200,18949,1],[20013,21809,1],[22675,23388,1]);

sub calcCoordinates {
    my $shiftFrame = $_[0];
    my $coordinateFrame = $_[1];
    my $locationshiftsLen = $_[2];
    my $oldShift = $_[3];
    my $i = $_[4];
    my $shiftFrameLimit = $locationshiftsLen-1;
    my $coordinate;
    my $oldFrame;
    if ($coordinates[$coordinateFrame][$i] < $locationshifts[$shiftFrame][0]) {
	$coordinate = $oldShift + $coordinates[$coordinateFrame][$i];
    }
    elsif ($locationshifts[$shiftFrame][0] == $coordinates[$coordinateFrame][$i]) {
        my $currentShift = $locationshifts[$shiftFrame][1];
        if ($currentShift == 0) {
          $coordinate = $coordinates[$coordinateFrame][$i];
        }
        elsif ($i == 0) {
          $coordinate = $oldShift + $coordinates[$coordinateFrame][$i];
        }
        elsif ($i == 1 && $currentShift > 0) {
          $coordinate = $currentShift + $coordinates[$coordinateFrame][$i];     
        }
	elsif ($i == 1 && $currentShift < 0) {
          $coordinate = $oldShift + $coordinates[$coordinateFrame][$i];     
        }
    }
    elsif ($coordinates[$coordinateFrame][$i] > $locationshifts[$shiftFrame][0] || $shiftFrame == $shiftFrameLimit) {
	until ($locationshifts[$shiftFrame][0] >= $coordinates[$coordinateFrame][$i] || $shiftFrame == $shiftFrameLimit) {
	    $oldShift = $locationshifts[$shiftFrame][1];
	    $shiftFrame++;
	}
	if ($shiftFrame == $shiftFrameLimit && $coordinates[$coordinateFrame][$i] < $locationshifts[$shiftFrame][0]) {
	    $coordinate = $coordinates[$coordinateFrame][$i] + $locationshifts[$shiftFrame-1][1];
	}
	elsif ($shiftFrame == $shiftFrameLimit && $coordinates[$coordinateFrame][$i] > $locationshifts[$shiftFrame][0]) {
	    $coordinate = $coordinates[$coordinateFrame][$i] + $locationshifts[$shiftFrame][1];
	}
	elsif ($locationshifts[$shiftFrame][0] == $coordinates[$coordinateFrame][$i]) {
	    if ($locationshifts[$shiftFrame][1] == 0) {
                $coordinate = $coordinates[$coordinateFrame][$i];
            }
            elsif ($i == 0) {
	        $oldFrame = $shiftFrame - 1;
                $coordinate = $locationshifts[$oldFrame][1] + $coordinates[$coordinateFrame][$i];
            }
            elsif ($i == 1 && $locationshifts[$shiftFrame][1] > 0) {
                $coordinate = $locationshifts[$shiftFrame][1] + $coordinates[$coordinateFrame][$i];     
            }
	    elsif ($i == 1 && $locationshifts[$shiftFrame][1] < 0) {
                $oldFrame = $shiftFrame - 1;
                $coordinate = $locationshifts[$oldFrame][1] + $coordinates[$coordinateFrame][$i];     
            }
        }
	else {
	    $coordinate = $oldShift + $coordinates[$coordinateFrame][$i];
	}
    }
    return ($coordinate);   
}


1;
