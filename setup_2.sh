#!/bin/bash
set -e  # Exit if any command fails

mamba create -y -n ngs -c conda-forge -c bioconda \
  samtools \
  fastqc \
  bwa \
  trimmomatic \
  seqkit

