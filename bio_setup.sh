#!/bin/bash
set -e  # Exit on error

echo "Creating directories..."
mkdir -p NGS
cd NGS
mkdir -p data/fastq reference adapters known_sites

echo "Downloading test FASTQ files (subset)..."
cd data
fastq-dump --split-files -X 1000000 --gzip SRR32281629
fastq-dump --split-files -X 1000000 --gzip SRR32281627
cd ../

echo "Downloading reference genome..."
cd reference
wget -c https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
echo "decompresing"
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
echo "indexing reference this might take a while ..."
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
# this is done in pipeline 
#echo "indexing using samtools"
#samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
#gatk CreateSequenceDictionary -R Homo_sapiens.GRCh38.dna.primary_assembly.fa
cd ..  


echo "Downloading known sites"
cd known_sites
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz .

echo "modifications of genome indexes for known sites to match reference"

echo "chr1   1
chr2    2
chr3    3
chr4    4
chr5    5
chr6    6
chr7    7
chr8    8
chr9    9
chr10   10
chr11   11
chr12   12
chr13   13
chr14   14
chr15   15
chr16   16
chr17   17
chr18   18
chr19   19
chr20   20
chr21   21
chr22   22
chrX    X
chrY    Y
chrM    MT"  > chr_rename.txt

echo "fixing known sites"
 
bcftools annotate \
  --rename-chrs chr_rename.txt \
  -Oz -o known_sites/Mills_fixed.vcf.gz \
  known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

echo "indexing fixed known sites"
tabix Mills_fixed.vcf.gz

rm Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

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
