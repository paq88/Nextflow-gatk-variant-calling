# Nextflow-gatk-variant-calling
setup.sh script should install all the neccesary components, 
In case you'd like to do that manually, you'll need:
mamba and mamba enviroment with:
- samtools
- fastqc
- bwa
- trimmomatic
- sra-tools
- bcftools
- gatk4
- r-base
-  r-ggplot2
-  r-gplots
-  r-gsalib
-  multiqc

`mamba create -y -n ngs -c conda-forge -c bioconda   samtools   fastqc   bwa   trimmomatic   sra-tools   bcftools   gatk4   r-base  r-ggplot2 r-gplots r-gsalib multiqc`

 To run the pipeline:
 - Nextflow
 - Java

also for test data preparation:
- unzip
- zip
- google-cloud-sdk

if you prefere automated setup 
`sudo ./setup.sh`
>>
`source .bashrc`
next you have to activate this enviroment:
 `mamba activate ngs`

If you'd like you can download test data for this pipeline, this script will download two human samples (only 1 000 000 bp each), human reference genome, index reference using bwa index, illumina adapters and it'll create directory structure. 
`./biosetup.sh`

Next you have to specify workDir and outdir in params section of nextflow pipeline 
And finally run the pipeline 

`nextflow run Variant_calling.nf`

