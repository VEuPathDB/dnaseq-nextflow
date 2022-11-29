FROM ubuntu:focal

# set environment variables
ENV varscan_version 2.3.9
ENV TABIX_VERSION 0.2.6

ENV CLASSPATH /usr/local/VarScan.jar

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y git ant build-essential wget unzip bcftools default-jre python3 python tabix samtools perl default-jre unzip cpanminus bioperl libaio1 emacs libjson-perl libmodule-install-rdf-perl libxml-parser-perl openjdk-8-jdk libdate-manip-perl libtext-csv-perl libstatistics-descriptive-perl libtree-dagnode-perl libxml-simple-perl && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

WORKDIR /gusApp
WORKDIR /gusApp/gus_home
WORKDIR /gusApp/project_home

ENV GUS_HOME=/gusApp/gus_home
ENV PROJECT_HOME=/gusApp/project_home
ENV PATH=$PROJECT_HOME/install/bin:$PATH
ENV PATH=$GUS_HOME/bin:$PATH

RUN export INSTALL_GIT_COMMIT_SHA=05197ebc4eb2046cc16e632b0b5852f21727a209 \
    && git clone https://github.com/VEuPathDB/install.git \
    && cd install \
    && git checkout $INSTALL_GIT_COMMIT_SHA

RUN mkdir -p $GUS_HOME/config && cp $PROJECT_HOME/install/gus.config.sample $GUS_HOME/config/gus.config

RUN export CBIL_GIT_COMMIT_SHA=41e17a8c7c61a6ca55fd28bd0f4883c74dcb625c \
    && git clone https://github.com/VEuPathDB/CBIL.git \
    && cd CBIL \
    && git checkout $CBIL_GIT_COMMIT_SHA \
    && bld CBIL

RUN export GUS_GIT_COMMIT_SHA=b11d5a179c5d48af134929c94b68bb908ab53bd6 \
    && git clone https://github.com/VEuPathDB/GusAppFramework.git \
    && mv GusAppFramework GUS \
    && cd GUS \
    && git checkout $GUS_GIT_COMMIT_SHA \
    && bld GUS/PluginMgr \
    && bld GUS/Supported

RUN export APICOMMONDATA_GIT_COMMIT_SHA=cf3e3daf2337a88462060439bd1fcf3a4b714c34 \
    && git clone https://github.com/VEuPathDB/ApiCommonData.git \
    && cd ApiCommonData \
    && git checkout $APICOMMONDATA_GIT_COMMIT_SHA \
    && mkdir -p $GUS_HOME/lib/perl/ApiCommonData/Load/Plugin \
    && cp $PROJECT_HOME/ApiCommonData/Load/plugin/perl/*.pm $GUS_HOME/lib/perl/ApiCommonData/Load/Plugin/ \
    && cp $PROJECT_HOME/ApiCommonData/Load/lib/perl/*.pm $GUS_HOME/lib/perl/ApiCommonData/Load/

RUN mkdir /gusApp/gus_home/lib/perl/GUS/Community \
    mkdir /gusApp/gus_home/lib/perl/VEuPath \
    && cp /gusApp/project_home/GUS/Community/lib/perl/GeneModelLocations.pm /gusApp/gus_home/lib/perl/GUS/Community/
ADD /lib/perl/* /gusApp/gus_home/lib/perl/VEuPath/
ENV PERL5LIB=/gusApp/gus_home/lib/perl

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

ADD /testing/lib/*.pm /usr/lib/x86_64-linux-gnu/perl5/5.30/VEuPath/

WORKDIR /work