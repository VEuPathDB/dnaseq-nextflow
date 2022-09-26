#!/usr/bin/perl

#======================================== DEPENDENCIES =========================================================================================
use lib "$ENV{GUS_HOME}/lib/perl";
use strict;
use Data::Dumper;
use Getopt::Long;
use GUS::ObjRelP::DbiDatabase;
use GUS::Supported::GusConfig;
use DBI;
use DBD::Oracle;

#===================================== ARGUMENTS AND SETUP ====================================================================================
my ($testFile, $gusConfigFile, $sequenceId, $regionFile);

&GetOptions("test_file=s"=> \$testFile,
            "sequence_id=s"=> \$sequenceId,
            "region_file=s"=> \$regionFile
           );

$gusConfigFile = $ENV{GUS_HOME} . "/config/gus.config";
my $gusConfig = GUS::Supported::GusConfig->new($gusConfigFile);

open(TEST, ">$testFile") or die "Cannot create a test file: $!";

my $db = GUS::ObjRelP::DbiDatabase->new($gusConfig->getDbiDsn(),
                                         $gusConfig->getDatabaseLogin(),
                                         $gusConfig->getDatabasePassword(),
                                         0, 0, 1,
                                         $gusConfig->getCoreSchemaName()
    );

my $dbh = $db->getQueryHandle();

#======================================= GENERATING LOCATIONSHIFTS ARRAY  =========================================================================

my $INDEL_QUERY = "select i.location, i.shift
                      from apidb.indel i
                      where sample_name = '$sequenceId'";

my $INDEL_QUERY_SH = $dbh->prepare($INDEL_QUERY);

$INDEL_QUERY_SH->execute();

my @locationshifts = ();
my $counter = 0;
my $currentShift = 0;
while (my ($location, $shift) = $INDEL_QUERY_SH->fetchrow_array()) {
    push ( @{$locationshifts[$counter]}, ($location, $shift + $currentShift));
    $counter++;
    $currentShift = $shift + $currentShift;
}
#print TEST Dumper(\@locationshifts), "\n";
#close TEST;

#====================================== GENERATING COORDINATES ARRAY ============================================================================
my $COORDINATE_QUERY = "select listagg(taf.source_id, ',') WITHIN GROUP (ORDER BY taf.aa_feature_id) as transcripts, s.source_id, tf.parent_id as gene_na_feature_id, el.start_min as exon_start, el.end_max as exon_end, decode(el.is_reversed, 1, afe.coding_end, afe.coding_start) as cds_start, decode(el.is_reversed, 1, afe.coding_start, afe.coding_end) as cds_end, el.is_reversed from dots.transcript tf , dots.translatedaafeature taf , dots.aafeatureexon afe , dots.exonfeature ef , dots.nalocation el, dots.nasequence s where tf.na_feature_id = taf.na_feature_id and taf.aa_feature_id = afe.aa_feature_id and afe.exon_feature_id = ef.na_feature_id and ef.na_feature_id = el.na_feature_id  and ef.na_sequence_id = s.na_sequence_id and tf.external_database_release_id = 5981 and s.source_id = 'OU755535' group by s.source_id, tf.parent_id, el.start_min, el.end_max, afe.coding_start, afe.coding_end, el.is_reversed order by s.source_id, el.start_min";
my $COORDINATE_QUERY_SH = $dbh->prepare($COORDINATE_QUERY);
$COORDINATE_QUERY_SH->execute();
my @coordinates = ();
$counter=0;
while (my ($transcript, $source_id, $gene_na_feature_id, $exon_start, $exon_end, $cds_start, $cds_end, $rev) = $COORDINATE_QUERY_SH->fetchrow_array()) {
  push (@{$coordinates[$counter]}, ($cds_start, $cds_end, $rev));
  $counter++;
}
#print TEST Dumper(\@coordinates), "\n";
#close TEST;

#===================================== GENERATING LOCATIONSSHIFTED ARRAY =========================================================================
my @locationsShifted = ();
my $shiftFrame=0;
my $coordinateFrame=0;
my $coordinatesLen = scalar @coordinates;
my $locationShiftsLen = scalar @locationshifts;
my $isRev;
my $oldShift = 0;
$counter = 0;
my ($cds_start, $cds_end, $isRev);
for (my $coordinateFrame;$coordinateFrame<$coordinatesLen;$coordinateFrame++){
    for (my $i=0;$i<3;$i++) {
        if ($i == 0) {
            ($cds_start, $oldShift, $shiftFrame) = &calcCoordinates($shiftFrame, $coordinateFrame, $locationShiftsLen, $oldShift, $i);
        }        
        elsif ($i == 1) {
            ($cds_end, $oldShift, $shiftFrame) = &calcCoordinates($shiftFrame, $coordinateFrame, $locationShiftsLen, $oldShift, $i);
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
#print TEST Dumper(\@locationsShifted), "\n";
#close TEST;

#========================== CREATING REGION FILE FOR SAMTOOLS ===========================================================================

open(REGION, ">$regionFile") or die "Cannot create a region file: $!";
my $locationsShiftedLen = scalar @locationsShifted;
my $adjustedLen = $locationsShiftedLen - 1;
my @lengthArray = (0..$adjustedLen);
for my $i (@lengthArray) {
    print REGION "$sequenceId:$locationsShifted[$i][0]-$locationsShifted[$i][1]\t$locationsShifted[$i][2]\n"; 
}
close REGION;

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


