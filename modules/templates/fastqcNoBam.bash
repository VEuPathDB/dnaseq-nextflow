#!/usr/bin/env bash

set -euo pipefail
mkdir fastqc_output   
fastqc -o fastqc_output --extract $sampleFile
