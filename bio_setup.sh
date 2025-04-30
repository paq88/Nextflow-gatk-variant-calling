#!/bin/bash
set -e  # Exit on error

echo "Creating directories..."
mkdir -p NGS
cd NGS
mkdir -p data/fastq reference adapters

echo "Downloading test FASTQ files (subset)..."
cd data
fastq-dump --split-files -X 1000000 --gzip SRR32281629
fastq-dump --split-files -X 1000000 --gzip SRR32281627
cd ../

echo "Downloading reference genome..."
cd reference
wget -c https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
echo "indexing reference this might take a while ..."
bwa index *

cd ..

echo "Downloading adapter sequences..."
cd adapters
wget -c https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa

echo "Creating placeholder for polyA adapter..."
echo ">polyA
AAAAAAAAAAAAAAAAAAAA" > polyA.fa
cat * > combined_adapters.fa
cd ..

echo "Setup part 3 complete."
