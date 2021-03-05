FROM ubuntu:20.04

# set environment variables
ENV hisat2_version 2.1.0
ENV trimmomatic_version 0.39
ENV fastqc_version 0.11.9

ENV CLASSPATH /usr/local/Trimmomatic-${trimmomatic_version}/trimmomatic-${trimmomatic_version}.jar

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update -y && apt-get install -y \
    wget \
    unzip \
    build-essential \
    default-jre \
    python3 \
    python 


# download software
WORKDIR /usr/local/
RUN wget https://github.com/DaehwanKimLab/hisat2/archive/v${hisat2_version}.tar.gz
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${fastqc_version}.zip
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${trimmomatic_version}.zip

# unpack
RUN tar -xvzf v${hisat2_version}.tar.gz
RUN unzip fastqc_v${fastqc_version}.zip
RUN unzip Trimmomatic-${trimmomatic_version}.zip

# # install hisat2
WORKDIR /usr/local/hisat2-${hisat2_version}
RUN make
RUN ln -s /usr/local/hisat2-${hisat2_version}/hisat2 /usr/local/bin/hisat2
RUN ln -s /usr/local/hisat2-${hisat2_version}/hisat2-build /usr/local/bin/hisat2-build

# install fastqc
RUN chmod gu+rx /usr/local/FastQC/fastqc
RUN ln -s /usr/local/FastQC/fastqc /usr/local/bin/fastqc

COPY ./data/trimmomatic_adaptors/* /usr/local/etc/