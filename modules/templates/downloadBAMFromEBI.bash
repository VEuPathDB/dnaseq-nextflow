#!/usr/bin/env bash

set -euo pipefail
wget --ftp-user $params.ebiFtpUser --ftp-password $params.ebiFtpPassword ftp://ftp-private.ebi.ac.uk:/upload/EBI/DNASeq/PlasmoDB_test/$params.organismAbbrev/Broad_HTS_Isolates_QuerySRA/$id/results.bam
mv results.bam ${id}.bam
