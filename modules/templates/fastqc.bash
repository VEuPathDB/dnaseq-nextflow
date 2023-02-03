#!/usr/bin/env bash

set -euo pipefail

if [ "$fromBam" = true ]; then

    mkdir fastqc_output
    
else

    mkdir fastqc_output   
    fastqc -o fastqc_output --extract $sampleFile
    
fi
