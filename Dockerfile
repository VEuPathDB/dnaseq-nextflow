FROM ubuntu:noble

# set environment variables
ENV TABIX_VERSION=0.2.6

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y git ant build-essential wget unzip bcftools python3 tabix samtools perl default-jre unzip cpanminus bioperl emacs libjson-perl libmodule-install-rdf-perl libxml-parser-perl libdate-manip-perl libtext-csv-perl libstatistics-descriptive-perl libtree-dagnode-perl libxml-simple-perl bwa trimmomatic openjdk-21-jre-headless && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

WORKDIR /gusApp/gus_home/lib/perl

ADD /lib/perl/* /gusApp/gus_home/lib/perl
ENV PERL5LIB=/gusApp/gus_home/lib/perl

# Install Perl modules (no database modules needed)
WORKDIR /usr/local/
RUN cpanm --verbose Bio::Coordinate::GeneMapper Bio::Coordinate::Pair VCF Set::CrossProduct Test2::V0

WORKDIR /usr/bin/

RUN wget https://github.com/freebayes/freebayes/releases/download/v1.3.10/freebayes-1.3.10-linux-amd64-static.gz \
  && gunzip freebayes-1.3.10-linux-amd64-static.gz \
  && chmod +x freebayes-1.3.10-linux-amd64-static \
  && mv freebayes-1.3.10-linux-amd64-static /usr/local/bin/freebayes

RUN wget https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip \
  && unzip snpEff_latest_core.zip \
  && cd snpEff \
  && rm snpEff.config

ADD /bin/snpEff.config /usr/bin/snpEff/snpEff.config

ADD /bin/* /usr/bin/
RUN cd /usr/bin \
  && chmod +x *.pl \
  && chmod +x *.sh

ADD /testing/lib/*.pm /usr/lib/x86_64-linux-gnu/perl5/5.30/VEuPath/

ADD /bin/All_adaptors* /usr/local/bin/

WORKDIR /work
