#!/usr/bin/perl

package shifted;
use strict;

my @locationshifts = ([1933,1],[2531,0],[3037,-2],[3254,0],[3433,-2],[8334,-1],[13340,-2],[13848,-4],[19255,-5],[20107,-8]);
my @coordinates = ([250,560,0],[3767,4765,1],[5853,7502,1],[9124,11130,1],[12136,12705,1],[15087,17084,1],[18200,18949,1],[20013,21809,1],[22675,23388,1]);

#===================================== GENERATING LOCATIONSSHIFTED ARRAY =========================================================================

sub locationShifted {
  my $locationShiftedFrame = @_[0];
  my $value = @_[1];
  my $oldShift = 0;  
  my $counter = 0;  
  my @locationsShifted = ();
  my $coordinateFrame=0;
  my $coordinatesLen = scalar @coordinates; 
  my $locationShiftsLen = scalar @locationshifts;
  my $shiftFrame=0;
  my $shiftFrameLimit = $locationShiftsLen-1;
  my $isRev;
  my ($cds_start, $cds_end, $isRev);
  for (my $coordinateFrame;$coordinateFrame<$coordinatesLen;$coordinateFrame++){
      for (my $i=0;$i<3;$i++) {
          if ($i == 0) {
              ($cds_start, $oldShift, $shiftFrame) = &calcCoordinates($shiftFrame, $coordinateFrame, $locationShiftsLen, $oldShift, $i, $shiftFrameLimit);
          }        
          elsif ($i == 1) {
              ($cds_end, $oldShift, $shiftFrame) = &calcCoordinates($shiftFrame, $coordinateFrame, $locationShiftsLen, $oldShift, $i, $shiftFrameLimit);
          }
          elsif ($i == 2) {
              $isRev = $coordinates[$coordinateFrame][$i];
	  }
          $oldShift = $oldShift;
          $shiftFrame = $shiftFrame;
      }
      $oldShift = $oldShift;
      $shiftFrame = $shiftFrame;
      push (@{$locationsShifted[$counter]}, ($cds_start, $cds_end, $isRev));
      $counter++;
  }
  return $locationsShifted[$locationShiftedFrame][$value];
}

#=================================== SUBROUTINES ==========================================================================================

sub calcCoordinates {
    my ($shiftFrame, $coordinateFrame, $locationshiftsLen, $oldShift, $i, $shiftFrameLimit) = @_;
    my $coordinate;
    if ($shiftFrame > $shiftFrameLimit) {
	$coordinate = $oldShift + $coordinates[$coordinateFrame][$i];
    }
    elsif ($coordinates[$coordinateFrame][$i] < $locationshifts[$shiftFrame][0]) {
	$coordinate = $oldShift + $coordinates[$coordinateFrame][$i];
    }
    elsif ($coordinates[$coordinateFrame][$i] == $locationshifts[$shiftFrame][0]) {
	$coordinate = $locationshifts[$shiftFrame][1] + $coordinates[$coordinateFrame][$i];
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
	    $coordinate = $locationshifts[$shiftFrame][1] + $coordinates[$coordinateFrame][$i];
	}
	else {
	    $coordinate = $oldShift + $coordinates[$coordinateFrame][$i];
	}
    }
    return ($coordinate, $oldShift, $shiftFrame);   
}
   
1;
