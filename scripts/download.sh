#!/bin/bash

accession=$1
prefetch --resume yes $accession -O ./
fastq-dump --split-files ./${accession}
gzip ${accession}*.fastq
rm -rf ${accession}

