if(params.fromBAM) {
    samples_ch = Channel.fromPath([params.inputDir + '/**/*.bam']).map { file -> tuple(file.baseName, [file]) }
}
else if(params.isPaired) {
   samples_ch = Channel.fromFilePairs([params.inputDir + '/**/*_{1,2}.fastq', params.inputDir + '/**/*_{1,2}.fastq.gz'])
}
else {
    samples_ch = Channel.fromPath([params.inputDir + '/**/*.fastq', params.inputDir + '/**/*.fastq.gz']).map { file -> tuple(file.baseName, [file]) }
}

process hisat2Index {
     input:
     	path genomeDir from params.genomeDir
     	val genomeFileOrIndexPrefix from params.genomeFileOrIndexPrefix

    output:
	path 'genomeIndex*.ht2' into hisat2_index_file_ch
        val 'genomeIndex' into hisat2_index_ch

    script: 
    if(params.fromBAM)
      """
      touch genomeIndex.1.ht2
      """
    else if( params.createIndex )
      """
      hisat2-build $genomeDir/$genomeFileOrIndexPrefix genomeIndex
      """
    else
      """
      TMP=$genomeDir/$genomeFileOrIndexPrefix
      FILES=\$TMP*
     for f in \$FILES; do cp "\$f" "genomeIndex\${f#\$TMP}" ; done
   
      """
}

process fastqc {
  input:
	tuple val(sampleName), path(sampleFile) from samples_ch

  output:
	tuple val(sampleName), path(sampleFile) into fastqc_samples_ch
        path('fastqc_output', type:'dir') into fastqc_dir_ch

   script:
    if(params.fromBAM)
	"""
        mkdir fastqc_output   
        """

        else 
	"""
        mkdir fastqc_output   
        fastqc -o fastqc_output --extract $sampleFile
        """
}

process fastqc_check {
  input:
	tuple val(sampleName), path(sampleFile) from fastqc_samples_ch
        path 'fastqc_output' from fastqc_dir_ch

  output:
	tuple val(sampleName), path(sampleFile) into fastqc_check_samples_ch
        path 'mateAEncoding' into fastqc_check_encoding_ch

  script:
    if(params.fromBAM)
         '''
           touch mateAEncoding
         '''

        else 
       '''
        #!/usr/bin/perl
        use strict;
        my @dataFiles = glob("fastqc_output/*/fastqc_data.txt");

        my %results;
        foreach(@dataFiles) {
          my $encoding = phred($_);
          $results{$encoding}++;
        }

        my @encodings = keys %results;
        if(scalar @encodings != 1) {
          die "Could not determine mateA encoding";
        }

        open(OUT, ">mateAEncoding") or die "Could not open mateAEncoding for writing: $!";
        print OUT $encodings[0] . "\\n";
        close OUT;

        sub phred {                                                                                                                                                          
            (my $fastqcFile) = @_;
            my %encoding = (
            "sanger" => "phred33",
            "illumina 1.3" => "phred64",
            "illumina 1.4" => "phred64",
            "illumina 1.5" => "phred64",
            "illumina 1.6" => "phred64",
            "illumina 1.7" => "phred64",
            "illumina 1.8" => "phred33",
            "illumina 1.9" => "phred33",
            "illumina 2"   => "phred33",
            "illumina 3"   => "phred33",
            "solexa" => "phred64"
            );
            my $phred;
            open (FH, $fastqcFile) or die "Cannot open $fastqcFile to determine phred encoding: $!";
            while (<FH>) {
                my $line = $_;
                if ($line =~ /encoding/i) {
                    foreach my $format (keys %encoding) {
                        if ($line =~ /$format/i) {
                            if (! defined $phred) {
                                $phred = $encoding{$format};
                            } elsif ($phred ne $encoding{$format}) {
                                $phred = "Error: more than one encoding type";
                            }
                        }
                    }
                    if (! defined $phred) {
                        $phred = "Error: format not recognized on encoding line";
                    }
                }
            }
            if (! defined $phred) {
                $phred = "Error: encoding line not found in file";
            }
            close(FH);
            if ($phred =~ /error/i) {
                die "ERROR: Could not determine phred encoding: $phred";
            }
            return $phred;
        }
       '''
}


process trimmomatic {
  input:
	tuple val(sampleName), path(sampleFile) from fastqc_check_samples_ch
        path('mateAEncoding') from fastqc_check_encoding_ch
        path adaptorsFile from params.trimmomaticAdaptorsFile

  output:
	tuple val(sampleName), path(sampleFile) into trimmomatic_samples_ch
        path 'mateAEncoding' into trimmomatic_encoding_ch
        path 'sample_1P' into trimmomatic_1p_ch
        path 'sample_2P' into trimmomatic_2p_ch

  script:
    if(params.fromBAM)
        """
        touch sample_1P
        touch sample_2P
        """
    else if( params.isPaired )
        """
         mateAEncoding=\$(<mateAEncoding)
         java org.usadellab.trimmomatic.TrimmomaticPE -trimlog trimLog.txt $sampleFile -\$mateAEncoding -baseout sample ILLUMINACLIP:$adaptorsFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
        """
    else 
	"""
        touch sample_2p
        mateAEncoding=\$(<mateAEncoding)
        java org.usadellab.trimmomatic.TrimmomaticSE -trimlog trimLog.txt -$mateAEncoding $sampleFile sample ILLUMINACLIP:$adaptorsFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
        """
}


process hisat2 {
  input:
	tuple val(sampleName), path(sampleFile) from trimmomatic_samples_ch
        path 'mateAEncoding' from trimmomatic_encoding_ch
        path 'sample_1p' from trimmomatic_1p_ch
        path 'sample_2p' from trimmomatic_2p_ch
        val hisat2_index from hisat2_index_ch
        path 'genomeIndex.*.ht2' from hisat2_index_file_ch

    output:
     tuple val(sampleName), path('result_sorted.bam') into hisat2_samples_ch

  script:
    if(params.fromBAM)
        """
        samtools view -bS $sampleFile | samtools sort - > result_sorted.bam
        """

    else if( params.isPaired )
        """
        mateAEncoding=\$(<mateAEncoding)
        hisat2 --no-spliced-alignment -k 1 -p $params.hisat2Threads -q --\$mateAEncoding -x $hisat2_index -1 sample_1p -2 sample_2p  | samtools view -bS - | samtools sort - > result_sorted.bam

        """

    else 
	"""
        mateAEncoding=\$(<mateAEncoding)
        hisat2 --no-spliced-alignment -k 1 -p $params.hisat2Threads -q --\$mateAEncoding -x $hisat2_index -U sample_1p | samtools view -bS - | samtools sort - > result_sorted.bam
        """
}


//fastqc_check_ch.view();
