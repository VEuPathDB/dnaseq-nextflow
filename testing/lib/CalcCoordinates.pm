#!/usr/bin/perl

package CalcCoordinates;
use strict;

my @locationshifts = ([1933,1],[2531,0],[3037,-2],[3254,0],[3433,-2],[8334,-1],[13340,-2],[13848,-4],[19255,-5],[20107,-8]);
my @coordinates = ([250,560,0],[3767,4765,1],[5853,7502,1],[9124,11130,1],[12136,12705,1],[15087,17084,1],[18200,18949,1],[20013,21809,1],[22675,23388,1]);
#0,8,10,0,0
sub calcCoordinates {
    my $shiftFrame = $_[0];
    my $coordinateFrame = $_[1];
    my $locationshiftsLen = $_[2];
    my $oldShift = $_[3];
    my $i = $_[4];
    my $coordinate;
    if ($shiftFrame > $locationshiftsLen) {
	$coordinate = $oldShift + $coordinates[$coordinateFrame][$i];
    }
    elsif ($coordinates[$coordinateFrame][$i] < $locationshifts[$shiftFrame][0]) {
	$coordinate = $oldShift + $coordinates[$coordinateFrame][$i];
    }
    elsif ($locationshifts[$shiftFrame][0] == $coordinates[$coordinateFrame][$i]) {
	$coordinate = $locationshifts[$shiftFrame][1] + $coordinates[$coordinateFrame][$i];
    }
    elsif ($locationshifts[$shiftFrame][0] < $coordinates[$coordinateFrame][$i]) {
	until ($locationshifts[$shiftFrame][0] >= $coordinates[$coordinateFrame][$i] || $shiftFrame >= $locationshiftsLen) {
	    $oldShift = $locationshifts[$shiftFrame][1];
	    $shiftFrame++;
	}
	if ($shiftFrame > $locationshiftsLen) {
	    $coordinate = $locationshifts[$shiftFrame][1] + $coordinates[$coordinateFrame][$i];
	}
	elsif ($locationshifts[$shiftFrame][0] == $coordinates[$coordinateFrame][$i]) {
	    $coordinate = $locationshifts[$shiftFrame][1] + $coordinates[$coordinateFrame][$i];
	}
	else {
	    $coordinate = $oldShift + $coordinates[$coordinateFrame][$i];
	}
    }
    return ($coordinate);   
}


1;
