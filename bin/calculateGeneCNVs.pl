#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use ApiCommonData::Load::CalculationsForCNVs;
use Data::Dumper;
use GUS::Supported::GusConfig;
use DBI;
use DBD::Oracle;

my ($gusConfigFile, $chrPloidyFile, $fpkmFile, $outputDir, $sampleName, $taxonId, $geneFootprintFile);
my $ploidy = 2;
&GetOptions("gusConfigFile=s" => \$gusConfigFile,
            "fpkmFile=s" => \$fpkmFile,
            "ploidy=i" => \$ploidy,
            "outputDir=s" => \$outputDir,
            "sampleName=s" => \$sampleName,
            "taxonId=i"=> \$taxonId,
            "geneFootPrints=s" => \$geneFootprintFile,
            );
 
&usage() unless ($fpkmFile && $outputDir && $sampleName && $taxonId && $geneFootprintFile);

my $sql = "with sequence as (
                select gf.source_id as gene_source_id
                , gf.na_feature_id
                , ns.source_id as contig_source_id
                , ns.source_id as sequence_source_id
                , ns.TAXON_ID
                from dots.genefeature gf
                , DOTS.NASEQUENCE ns
                , SRES.ONTOLOGYTERM ot
                where gf.na_sequence_id = ns.na_sequence_id
                and ot.name = 'chromosome'
                and ns.SEQUENCE_ONTOLOGY_ID = ot.ONTOLOGY_TERM_ID
                and ns.taxon_id = $taxonId
                union
                select gf.SOURCE_ID AS gene_source_id
                , gf.na_feature_id
                , contig.source_id as contig_source_id
                , scaffold.SOURCE_ID as chromosome_source_id
                , scaffold.TAXON_ID
                from dots.genefeature gf
                , dots.nasequence contig
                , dots.nasequence scaffold
                , dots.sequencepiece sp
                , SRES.ONTOLOGYTERM ot
                where gf.na_sequence_id = sp.PIECE_NA_SEQUENCE_ID
                and gf.na_sequence_id = contig.na_sequence_id
                and sp.VIRTUAL_NA_SEQUENCE_ID = scaffold.NA_SEQUENCE_ID
                and ot.name = 'chromosome'
                and scaffold.SEQUENCE_ONTOLOGY_ID = ot.ONTOLOGY_TERM_ID
                and scaffold.taxon_id = $taxonId
            )
  
            , orthologs as (
                select gf.na_feature_id, sg.name
                from dots.genefeature gf
                , dots.SequenceSequenceGroup ssg
                , dots.SequenceGroup sg
                , core.TableInfo ti
                where gf.na_feature_id = ssg.sequence_id
                and ssg.sequence_group_id = sg.sequence_group_id
                and ssg.source_table_id = ti.table_id
                and ti.name = 'GeneFeature'
            )

            select s.gene_source_id
            , o.name
            from sequence s
            , orthologs o
            where s.na_feature_id = o.na_feature_id";

#get hashrefs of chromosome data
my $chrs = ApiCommonData::Load::CalculationsForCNVs::getChrsForCalcs($taxonId);
my $geneData = ApiCommonData::Load::CalculationsForCNVs::getGeneInfo($geneFootprintFile, $chrs);
my $chrValues = ApiCommonData::Load::CalculationsForCNVs::getChrFPKMVals($fpkmFile, $chrs, $geneData);
my $chrMedians = ApiCommonData::Load::CalculationsForCNVs::getChrMedians($chrValues, $chrs);
my $allChrMedian = ApiCommonData::Load::CalculationsForCNVs::getMedianAcrossChrs($chrValues, $chrs);
my $chrPloidies = ApiCommonData::Load::CalculationsForCNVs::getChrPloidies($chrMedians, $allChrMedian, $ploidy, $chrs);

# get OrthoMCL IDs from DB
if (!$gusConfigFile) {
    $gusConfigFile = $ENV{GUS_HOME}."/config/gus.config";
}

die "Config file $gusConfigFile does not exist" unless -e $gusConfigFile;

my $gusConfig = GUS::Supported::GusConfig->new($gusConfigFile);
my ($dbiDsn, $login, $password, $core);
$dbiDsn = $gusConfig->getDbiDsn();
$login = $gusConfig->getDatabaseLogin();
$password = $gusConfig->getDatabasePassword();
$core = $gusConfig->getCoreSchemaName();
my $dbh = DBI->connect($dbiDsn, $login, $password);

my $orthoMclStmt = $dbh->prepare($sql);
$orthoMclStmt->execute();
my $geneOrthoMclIds = {};
while (my @row = $orthoMclStmt->fetchrow_array()){
    $geneOrthoMclIds->{$row[0]} = $row[1];
}
die "ERROR:\tNo orthologs can be found for genes in the organism with taxon id $taxonId. Ortholog data is required to run this script.\n\n
         DATA LOADERS: Please undo this step and re-run after orthologs have been loaded in this database.\n"
         unless scalar keys %$geneOrthoMclIds > 0;
$dbh->disconnect();

my $allGenesData = &geneCopyNumberCalc($fpkmFile, $chrPloidies, $chrMedians, $geneOrthoMclIds, $chrs, $geneData);
&calcCnvForOrthologGroups($allGenesData, $chrPloidies, $outputDir, $sampleName);
 
sub geneCopyNumberCalc {
    my ($fpkmFile, $chrPloidies, $chrMedians, $geneOrthoMclIds, $chrs, $geneInfo) = @_;
    my $allGenesData = [];
    open (IN, $fpkmFile) or die "Cannot read FPKM file $fpkmFile\n$!\n";

    while (<IN>) {
        my $line = $_;
        chomp ($line);

        #ignore header in fpkm file 
        if ($line=~/^PROJECT\t/){
            next;
        }
        else{
            my @data = split(/\t/, $line);
            my ($geneId, $fpkm) = @data;
            die "Error: Gene and FPKM values are missing in line $line\n" unless (defined($geneId) && defined($fpkm));
            my $chr;
            if (exists $geneInfo->{$geneId}) {
                $chr = $geneInfo->{$geneId}->{'chr'};
                my $orthoMclId;
                if (exists $geneOrthoMclIds->{$geneId}) {
                    $orthoMclId = $geneOrthoMclIds->{$geneId};
                } else {
                    $orthoMclId = "UNK_.$geneId";
                }

                if ($chrMedians->{$chr}==0){
                    print STDERR "Error: Bad median for chromosome $chr.  No copy number calculation for $geneId\n";
                    next;
                }
                if ($chrPloidies->{$chr}==0) {
                    print STDERR "Chromosome $chr has a predicted ploidy of 0. The copy number for gene $geneId on chromosome $chr will therefore also be predicted to be 0";
                    my $geneData = [$orthoMclId, $geneId, $chr, 0, $fpkm];
                    push(@{$allGenesData}, $geneData);
                } else {

                    my $rawCN = $fpkm/($chrMedians->{$chr}/$chrPloidies->{$chr});
                    my $geneData = [$orthoMclId, $geneId, $chr, $rawCN, $fpkm];
                    push (@{$allGenesData}, $geneData);
                }
            }
        }
    }
    return $allGenesData;
}

sub calcCnvForOrthologGroups {
    my ($allGenesData, $chrPloidies, $outputDir, $sampleName) = @_;
    my %clusterData;
    # cluster genes that share an orthoMclId
    foreach my $geneData (@{$allGenesData}) {
        my ($orthoMclId, $geneId, $chr, $rawCN, $fpkm) = @{$geneData};
        my $id = $orthoMclId;
        $clusterData{$id}->{"chr"} = $chr;
        push @{$clusterData{$id}->{"genes"}}, $geneId;
        push @{$clusterData{$id}->{"rawCNs"}}, $rawCN;
        $clusterData{$id}->{"rawCNTotal"} += $rawCN;
        push @{$clusterData{$id}->{"fpkms"}}, $fpkm;
    }
    #TSV will contain data as produced by Nick's original script.  TXT contains data in the correct format for loading
    open (TSV, ">$outputDir/$sampleName\_CNVestimations.tsv") or die "Cannot open output file $outputDir/$sampleName\_CNVestimations.tsv for writing\n$!\n";
    open (TXT, ">$outputDir/$sampleName\_geneCNVs.txt") or die "Cannot open outputfile $outputDir/$sampleName\_geneCNVs.txt for writing\n$!\n";

    # Print header lines in both files
    print TSV "OrthoMclId\tChromosome\tPloidy\tGenes In Ref\tTotal CN\tHaploid Number\tGene Dose\tGene List\n";
    print TXT "na_feature_id\thaploid_number\tref_copy_number\n";

    # calculate values for clusters
    foreach my $cluster (keys %clusterData) {
        my $orthoMclId = $cluster;
        my $chr = $clusterData{$cluster}->{"chr"};
        my $chrPloidy = $chrPloidies->{$chr};
        next if $chrPloidy == 0;
        my $haploidCN = $clusterData{$cluster}->{"rawCNTotal"}/$chrPloidy;
        my $geneDose = $haploidCN * $chrPloidy;
        my $genesInRef = scalar(@{$clusterData{$cluster}->{"genes"}});
        my $geneList = join(", ", @{$clusterData{$cluster}->{"genes"}});

        print TSV "$orthoMclId\t$chr\t$chrPloidy\t$genesInRef\t".sprintf("%4.3f", $clusterData{$cluster}->{"rawCNTotal"})."\t".sprintf("%4.2f", $haploidCN)."\t".sprintf("%3.1f", $geneDose)."\t$geneList\n";

        foreach my $gene (@{$clusterData{$cluster}->{"genes"}}) {
            print TXT "$gene\t$haploidCN\t$genesInRef\n";
        }
    }
    close (TSV);
    close (TXT);
}
 
       
sub usage {
    print STDERR "calculateGeneCNVs --fpkmFile <path to genes.fpkm_tracking file generated by Cufflinks> --outputDir <Dir to write output files> --sampleName <sample name> --sql <sql statement to return a list of genes and their OrthoMCL_Ids> --geneFootPrints <gene footprint file> [--gusConfigFile <supply if not using default>] --ploidy <Base ploidy for this organism (i.e., what you expect the majority of chromosomes to be) - default is 2>";
    exit;
}

exit;
