FROM ubuntu:focal

# set environment variables
ENV varscan_version 2.3.9
ENV TABIX_VERSION 0.2.6

ENV CLASSPATH /usr/local/VarScan.jar

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y build-essential wget unzip bcftools default-jre python3 python tabix samtools perl default-jre unzip && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

# download software
WORKDIR /usr/local/
RUN wget -O VarScan.jar https://sourceforge.net/projects/varscan/files/VarScan.v${varscan_version}.jar/download

WORKDIR /usr/bin/

RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
  && unzip snpEff_latest_core.zip \
  && cd snpEff \
  && rm snpEff.config

ADD /bin/snpEff.config /usr/bin/snpEff/snpEff.config

ADD /bin/*.pl /usr/bin/
