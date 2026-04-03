package CalculationsForTPM;

use strict;
use warnings;

use Exporter qw(import);
our @EXPORT_OK = qw(doTPMCalculation calcTPMForCNVs);

sub _getCounts {
    my %specialCounters = ('__no_feature'=>1, '__ambiguous'=>1, '__too_low_aQual'=>1, '__not_aligned'=>1, '__alignment_not_unique'=>1);
    my ($uniqueCounts, $nonuniqueCounts, $isCNV) = @_;
    my $uniqueCountHash = {};
    my $nonuniqueCountHash = {};

    open (UNIQUE, "<$uniqueCounts") or die "Cannot open file $uniqueCounts. Please check and try again\n$!\n";
    while (<UNIQUE>) {
        my ($geneId, $count) = split /\t/, $_;
        chomp $count;
        next if ($specialCounters{$geneId});
        $uniqueCountHash->{$geneId} = $count;
    }
    close UNIQUE;

    if (defined $isCNV) {
        return $uniqueCountHash;
    }

    open (NONUNIQUE, "<$nonuniqueCounts") or die "Cannot open file $nonuniqueCounts. Please check and try again\n$!\n";
    while (<NONUNIQUE>) {
        my ($geneId, $count) = split /\t/, $_;
        chomp $count;
        next if ($specialCounters{$geneId});

        $count = $count-$uniqueCountHash->{$geneId};
        #NU counts from htseq-count include unique counts too. Subtract unique counts to get only NU
        $nonuniqueCountHash->{$geneId} = $count;
    }
    close NONUNIQUE;
    return($uniqueCountHash, $nonuniqueCountHash);
}

sub _calcRPK {
    my($geneLengths, $countHash, $rpkSum) = @_;
    my $rpkHash;
    while (my ($geneId, $count) = each %{$countHash}) {
        if ($geneLengths->{$geneId}) {
            my $geneLength = $geneLengths->{$geneId}/1000;
            my $rpk = $count/$geneLength;
            $rpkSum += $rpk;
            $rpkHash->{$geneId} = $rpk;
        } else {
            print STDERR "WARNING: Gene $geneId was not found in footprint file. No data will be loaded for this gene.\n";
        }
    }
    return ($rpkSum, $rpkHash);
}

sub _calcTPM {
    my ($rpkHash, $rpkSum) = @_;
    $rpkSum = $rpkSum/1000000;
    my $tpmHash;
    while (my($geneId, $rpk) = each %{$rpkHash}) {
        my $tpm = $rpk/$rpkSum;
        $tpmHash->{$geneId} = $tpm;
    }
    return $tpmHash;
}

#add together for NU to mimic how FPKM profiles look for GraphPackage
sub _nonUniqueTPM {
    my ($uniqueTpmHash, $nonUniqueTpmHash) = @_;
    my $newHash = {};
    foreach my $geneId (keys %{$nonUniqueTpmHash}) {
        if (exists $uniqueTpmHash->{$geneId}) {
            $newHash->{$geneId} = ($nonUniqueTpmHash->{$geneId} + $uniqueTpmHash->{$geneId});
        } else {
            die "Gene $geneId present in non-unique counts hash is not present in unique counts hash\n"
        }
    }
    return $newHash;
}

sub _writeTPM {
    my ($tpmFile, $tpmHash) = @_;
    open (OUT, ">$tpmFile") or die "Cannot open TPM file $tpmFile for writing. Please check and try again.\n$!\n";
    while (my ($geneId, $tpm) = each %{$tpmHash}) {
        print OUT ("$geneId\t$tpm\n") ;
    }
    close OUT;
}

# sub for RNAseq using unique and non-unique with option for sense and antisense
sub doTPMCalculation {
    my ($geneLengths, $senseUniqueCountFile, $senseNUCountFile, $antisenseUniqueCountFile, $antisenseNUCountFile, $senseUniqueTpmFile, $senseNUTpmFile, $antisenseUniqueTpmFile, $antisenseNUTpmFile) = @_;

    # get sense counts
    my ($senseUniqueCountHash, $senseNUCountHash) = &_getCounts($senseUniqueCountFile, $senseNUCountFile);

    #calculate sense RPK
    my ($rpkSum, $senseUniqueRpkHash) = &_calcRPK($geneLengths, $senseUniqueCountHash, 0);
    ($rpkSum, my $senseNURpkHash) = &_calcRPK($geneLengths, $senseNUCountHash, $rpkSum);

    #if we have antisense counts
    if ($antisenseUniqueCountFile && $antisenseNUCountFile) {
        if ($antisenseUniqueTpmFile && $antisenseNUTpmFile) {
            # get antisense counts
            my ($antisenseUniqueCountHash, $antisenseNUCountHash) = &_getCounts($antisenseUniqueCountFile, $antisenseNUCountFile);

            #calculate antisense RPK (must be  done before any TPM calcs to ensure RPK sum includes all values)
            ($rpkSum, my $antisenseUniqueRpkHash) = &_calcRPK($geneLengths, $antisenseUniqueCountHash, $rpkSum);
            ($rpkSum, my $antisenseNURpkHash) = &_calcRPK($geneLengths, $antisenseNUCountHash, $rpkSum);

            #calculate and write out antisense TPM while we are in this loop
            my $antisenseUniqueTpmHash = &_calcTPM($antisenseUniqueRpkHash, $rpkSum);
            my $antisenseNUTpmHash = &_calcTPM($antisenseNURpkHash, $rpkSum);
            $antisenseNUTpmHash = &_nonUniqueTPM($antisenseUniqueTpmHash, $antisenseNUTpmHash);
            &_writeTPM($antisenseUniqueTpmFile, $antisenseUniqueTpmHash);
            &_writeTPM($antisenseNUTpmFile, $antisenseNUTpmHash);
        } else {
            die "Antisense count files $antisenseUniqueCountFile  and $antisenseNUCountFile have been provided, but one or both antisense TPM files have been specified for writing output are missing. Please add these using the --antisenseUniqueTpmFile and --antisenseNUTpmFile flags.\n";
        }
    }

    # calculate and write out sense TPM
    my $senseUniqueTpmHash = &_calcTPM($senseUniqueRpkHash, $rpkSum);
    my $senseNUTpmHash = &_calcTPM($senseNURpkHash, $rpkSum);
    $senseNUTpmHash = &_nonUniqueTPM($senseUniqueTpmHash, $senseNUTpmHash);
    &_writeTPM($senseUniqueTpmFile, $senseUniqueTpmHash);
    &_writeTPM($senseNUTpmFile, $senseNUTpmHash);
}

# sub for CNVs with just one set of counts
sub calcTPMForCNVs {
    my ($geneLengths, $countFile, $tpmFile) = @_;
    my $countHash = &_getCounts($countFile, '', 1);
    my ($rpkSum, $rpkHash) = &_calcRPK($geneLengths, $countHash, 0);
    my $tpmHash = &_calcTPM($rpkHash, $rpkSum);
    &_writeTPM($tpmFile, $tpmHash);
}

1;
