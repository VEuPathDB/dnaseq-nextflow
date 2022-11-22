FROM ubuntu:focal

# set environment variables
ENV varscan_version 2.3.9
ENV TABIX_VERSION 0.2.6

ENV CLASSPATH /usr/local/VarScan.jar

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y build-essential wget unzip bcftools default-jre python3 python tabix samtools perl default-jre unzip cpanminus bioperl libaio1 emacs && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

WORKDIR /opt/oracle
RUN export INSTANTCLIENT_VER=linux.x64-21.6.0.0.0dbru \
    && wget https://download.oracle.com/otn_software/linux/instantclient/216000/instantclient-basic-$INSTANTCLIENT_VER.zip \
    && wget https://download.oracle.com/otn_software/linux/instantclient/216000/instantclient-tools-$INSTANTCLIENT_VER.zip \
    && wget https://download.oracle.com/otn_software/linux/instantclient/216000/instantclient-sdk-$INSTANTCLIENT_VER.zip \
    && unzip instantclient-basic-$INSTANTCLIENT_VER.zip \
    && unzip instantclient-tools-$INSTANTCLIENT_VER.zip \
    && unzip instantclient-sdk-$INSTANTCLIENT_VER.zip
# Need to change this if we get new version of the instantclient
ENV ORACLE_HOME=/opt/oracle/instantclient_21_6
ENV LD_LIBRARY_PATH=$ORACLE_HOME
ENV PATH=/opt/oracle/instantclient_21_6:$PATH
WORKDIR /tmp
RUN export DBI_VER=1.643 \
    && wget http://www.cpan.org/modules/by-module/DBI/DBI-$DBI_VER.tar.gz \
    && tar xvfz DBI-$DBI_VER.tar.gz \
    && cd DBI-$DBI_VER \
    && perl Makefile.PL \
    && make \
    && make install
RUN export ORACLE_DBD_VER=1.83 \
    && wget http://www.cpan.org/modules/by-module/DBD/DBD-Oracle-$ORACLE_DBD_VER.tar.gz \
    && tar xvfz DBD-Oracle-$ORACLE_DBD_VER.tar.gz \
    && cd DBD-Oracle-$ORACLE_DBD_VER \
    && perl Makefile.PL \
    && make \
    && make install
WORKDIR /work

# download software
WORKDIR /usr/local/
RUN wget -O VarScan.jar https://sourceforge.net/projects/varscan/files/VarScan.v${varscan_version}.jar/download

RUN cpanm Bio::Perl Bio::Seq Bio::Tools::GFF Bio::Coordinate::GeneMapper Bio::Coordinate::Pair Bio::Location::Simple Bio::Tools::CodonTable VCF DBD::Oracle DBI Set::CrossProduct Test2::V0

WORKDIR /usr/bin/

RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
  && unzip snpEff_latest_core.zip \
  && cd snpEff \
  && rm snpEff.config

ADD /bin/snpEff.config /usr/bin/snpEff/snpEff.config

ADD /bin/* /usr/bin/
RUN cd /usr/bin \
  && chmod +x *.pl \
  && chmod +x *.sh

ADD /lib/perl/* /usr/lib/x86_64-linux-gnu/perl5/5.30/VEuPath/
ADD /testing/lib/*.pm /usr/lib/x86_64-linux-gnu/perl5/5.30/VEuPath/

WORKDIR /work