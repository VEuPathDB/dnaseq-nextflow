#!/usr/bin/env perl

use strict;
use warnings;
use lib "$ENV{GUS_HOME}/lib/perl";
use Getopt::Long;
use GUS::ObjRelP::DbiDatabase;
use GUS::Supported::GusConfig;

# get parameter values
my ($variationFile, $gusConfigFile, $extDbRlsSpec);

&GetOptions("variationFile=s" => \$variationFile,
            "gusConfig=s" => \$gusConfigFile,
            "extdb_spec=s" => \$extDbRlsSpec,
           );

my $gusConfig = GUS::Supported::GusConfig->new($gusConfigFile);

my $db = GUS::ObjRelP::DbiDatabase->new($gusConfig->getDbiDsn(),
                                         $gusConfig->getDatabaseLogin(),
                                         $gusConfig->getDatabasePassword(),
                                         0, 0, 1,
                                         $gusConfig->getCoreSchemaName()
    );
my $dbh = $db->getQueryHandle();

my $thisExtDbRlsId = &queryExtDbRlsIdFromSpec($dbh, $extDbRlsSpec);

my ($outputFile);
open($outputFile, "> ./variationFeature.dat") or die "Cannot open file variationFeature for writing: $!";

open(my $data, '<', $variationFile) || die "Could not open file $variationFile: $!";
while (my $line = <$data>) {
    chomp $line;
    my ($location, $geneNaFeatureId, $sourceId, $referenceStrain, $referenceNa, $referenceAa, $hasNonsynonymousAllele, $majorAllele, $minorAllele, $majorAlleleCount, $minorAlleleCount, $majorProduct, $minorProduct, $distinctStrainCount, $distinctAlleleCount, $hasCodingMutation, $totalAlleleCount, $hasStopCodon, $refCodon, $referenceAAFull) = split(/\t/, $line);
    print $outputFile "$thisExtDbRlsId\t$location\t$geneNaFeatureId\t$sourceId\t$referenceStrain\t$referenceNa\t$referenceAa\t$hasNonsynonymousAllele\t$majorAllele\t$minorAllele\t$majorAlleleCount\t$minorAlleleCount\t$majorProduct\t$minorProduct\t$distinctStrainCount\t$distinctAlleleCount\t$hasCodingMutation\t$totalAlleleCount\t$hasStopCodon\t$refCodon\t$referenceAAFull\n"
}


close $outputFile;
close $data;

sub queryExtDbRlsIdFromSpec {
    my ($dbh, $extDbRlsSpec) = @_;

  my $sql = "select r.external_database_release_id
from sres.externaldatabaserelease r, sres.externaldatabase d
where d.name || '|' || r.version = '$extDbRlsSpec'";


    my $sh = $dbh->prepare($sql);
    $sh->execute();

    my ($extDbRlsId) = $sh->fetchrow_array();

    $sh->finish();
    die "Could not find ext db rls id for spec: $extDbRlsSpec" unless($extDbRlsId);

    return $extDbRlsId;
}
