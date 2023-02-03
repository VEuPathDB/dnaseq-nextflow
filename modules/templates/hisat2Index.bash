#!/usr/bin/env bash

set -euo pipefail

if [ "$fromBam" = true ]; then

    touch genomeIndex.1.ht2
    samtools faidx $genomeFasta
    
elif [ "$createIndex" = true ]; then

    hisat2-build $genomeFasta genomeIndex
    samtools faidx $genomeFasta
    
else

    TMP=$params.hisat2Index
    FILES=\$TMP*
    for f in \$FILES; do cp "\$f" "genomeIndex\${f#\$TMP}" ; done
    samtools faidx $genomeFasta

fi
