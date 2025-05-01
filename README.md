# Nextflow-gatk-variant-calling


This will install all the neccesary software
- java
- nextflow
- conda
- mamba

as well as crete mamba env with:
- samtools
- fastqc
- trimmomatic
- bwa
- sra-tools

`./setup.sh`

next you have to activate this enviroment:

`mamba activate ngs`
If you'd like you can download test data for this pipeline, this script will download two human samples (only 1 000 000 bp each), human reference genome, index reference using bwa index, illumina adapters and it'll create directory structure. 
`./biosetup.sh`

Next you have to specify workDir and outdir in params section of nextflow pipeline 
And finally run the pipeline 

`nextflow run Variant_calling.nf`

