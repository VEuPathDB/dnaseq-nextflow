if(params.isPaired) {
   samples_ch = Channel.fromFilePairs(['data/fastq_paired_haploid/**/*_{1,2}.fastq','data/fastq_paired_haploid/**/*_{1,2}.fastq.gz'])
}
else {
   samples_ch = Channel.fromPath(['data/fastq_paired_haploid/**/*.fastq', 'data/fastq_paired_haploid/**/*.fastq.gz']).map { file -> tuple(file.baseName, [file]) }
}

process fastqc {
  input:
	tuple val(sampleName), path(readsFn) from samples_ch

  output:
	tuple path('fastqc_output', type:'dir'), val(sampleName), path(readsFn) into fastqc_ch

   script:
	"""
        mkdir fastqc_output   
        fastqc -o fastqc_output --extract $readsFn
        """
}

process fastqc_check {
  input:
	tuple path('fastqc_output'), val(sampleName), path(readsFn) from fastqc_ch

  output:
	tuple val(sampleName), path(readsFn), path('mateAEncoding') into fastqc_check_ch

  script:

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

fastqc_check_ch.view();
