#!/usr/bin/env bash

set -euo pipefail

if [ "$fromBam" = true ]; then

    samtools view -bS $sampleFile | samtools sort - > result_sorted.bam
    
elif [ "$isPaired" = true ]; then

    mateAEncoding=\$(<mateAEncoding)
    hisat2 --no-spliced-alignment \
        -k 1 \
        -p $params.hisat2Threads \
        -q --\$mateAEncoding \
        -x $hisat2_index \
        -1 $sample_1p \
        -2 $sample_2p  \
        | samtools collate -o output.bam -
        samtools fixmate -m output.bam fix.bam
        samtools sort -o sort.bam fix.bam
        samtools markdup -r sort.bam result_sorted.bam
    
else

    mateAEncoding=\$(<mateAEncoding)
    hisat2 --no-spliced-alignment \
        -k 1 \
        -p $params.hisat2Threads \
        -q --\$mateAEncoding \
        -x $hisat2_index \
        -U $sample_1p \
        | samtools collate -o output.bam -
        samtools fixmate -m output.bam fix.bam
        samtools sort -o sort.bam fix.bam
        samtools markdup -r sort.bam result_sorted.bam

fi
