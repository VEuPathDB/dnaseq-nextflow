#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process checkUniqueIds {

  input:
    path params.fastaDir

  output:
    path 'check.txt'

  script:
    template 'checkUniqueIds.bash'
}



process generateRegion {

  input:
    path 'makepositionarraycoding.pl'
    path 'input.fasta'
    
  output:
    path 'input.fasta'
    path 'region.txt'

  script:
    """
    grep ">" input.fasta > temp.fa
    DEFLINE=\$(sed 's/>//' temp.fa)
    /usr/bin/perl makepositionarraycoding.pl --test_file shifted.txt --sequence_id \$DEFLINE --region_file region.txt
    """
}


process runSamtools {

  input:
    path 'genome.fa'
    path 'region.txt'

  output:
    path 'transcriptFinal.fasta'

  script:
    '''
    #!/usr/bin/env perl

use strict;

my $defline;
open(REGION, "<region.txt") or die "Couldn't open regionFile";
while(<REGION>){
    if ($_ =~ /^(.*):.+/) {
         $defline = $1;
	  last;
    }
}
close REGION;

open(REGION, "<region.txt") or die "Couldn't open regionFile";
my $line;
my $seq;
while(<REGION>){
    if ($_ =~ /\t0/) {
    $line = $_;
    $line =~ s/\t0//g;
    open(TEMP, ">temp.txt")  or die "Couldn't open tempFile";
    print TEMP "$line";
    close TEMP;
    $seq = `samtools faidx -r temp.txt genome.fa`;
    }
    else {
    $line = $_;
    $line =~ s/\t1//g;
    open(TEMP, ">temp.txt")  or die "Couldn't open tempFile";
        print TEMP "$line";
	close TEMP;
	$seq = `samtools faidx -r temp.txt -i genome.fa`;
    }
    open(FASTA, ">>transcript.fasta") or die "Couldn't open fastaFile";
    print FASTA $seq;
    close FASTA;
}
close REGION;

open(O,">temp.fasta");

open(I,"<transcript.fasta") || die "Unable to open transcriptFile";

my $line;
while(<I>){
    if (/^>/) {
    print "Defline";
    }
    else {
    $line = $_;
    chomp($line);
    print O $line;
    }
}
close I;
close O;
my $fasta = `fold -w 60 temp.fasta > transcriptFinal.fasta`;
$fasta = `echo '>$defline' | cat - transcriptFinal.fasta > temp && mv temp transcriptFinal.fasta`;
    '''
}


process makeIndex {

  publishDir "$params.outputDir", mode: "copy", pattern: 'combinedFasta.fa.fai'
  publishDir "$params.outputDir", mode: "copy", pattern: 'combinedFasta.fa'

  input:
    path ('combinedFasta.fa')
    path 'check.txt'

  output:
    path('combinedFasta.fa.fai')
    path('combinedFasta.fa') 

  script:
    template 'makeIndex.bash'
}


process mergeVcfs {
  container = 'biocontainers/bcftools:v1.9-1-deb_cv1'

  publishDir "$params.outputDir", mode: "copy", pattern: 'merged.vcf.gz'

  input:
    path '*.vcf.gz'
    path '*.vcf.gz.tbi'

  output:
    path 'merged.vcf.gz'
    path 'toSnpEff.vcf'

  script:
    template 'mergeVcfs.bash'
}


process snpEff {

  publishDir "$params.outputDir", mode: "copy"

  input:
    path 'merged.vcf'
    path 'genes.gtf.gz'
    path 'sequences.fa.gz'

  output:
    path 'merged.ann.vcf'

  script:
    template 'snpEff.bash'    
}


workflow AllExperiments {

  take:
    fastas_qch
    vcfs_qch
    vcfsindex_qch
    
  main:
    
    checkResults = checkUniqueIds(params.fastaDir) 
    generateRegion(params.makepositionarraycoding, fastas_qch) | runSamtools | collectFile(storeDir: params.outputDir, name: 'transcriptFinal.fa', newLine: true)
    combinedFasta = fastas_qch.collectFile(name: 'CombinedFasta.fa')
    makeIndex(combinedFasta, checkResults)
    allvcfs = vcfs_qch.collect()
    allvcfindexes = vcfsindex_qch.collect()
    mergeVcfsResults = mergeVcfs(allvcfs, allvcfindexes)
    snpEff(mergeVcfsResults[1], params.databaseFile, params.sequenceFile)

}
