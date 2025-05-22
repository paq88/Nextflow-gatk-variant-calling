#!/usr/bin/env nextflow

// Define the working directory
params.workDir = '/home/paq/NGS/NGS'
params.outdir  = '/home/paq/NGS/NGS/outdir'

// Define input files (paired-end reads)

params.reference = "${params.workDir}/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"


params.adapters = "${params.workDir}/adapters/combined_adapters.fa"
params.reads = "${params.workDir}/data/*_{1,2}.fastq.gz"
params.knownSites = "${params.workDir}/known_sites/Mills_fixed.vcf.gz"

params.trimCPUs = 2
params.bwaCPUs = 4
params.gatkCPUs = 4
params.bwaMemory = '7 GB '

// Create a channel for paired-end reads


// data quality control
process qcBefore {
    publishDir "${params.outdir}/qc/before_trimming/fastqc", mode: 'copy', overwrite: true
    input: 
    tuple val(id), path(reads) 

    output:
    tuple path("*_1_fastqc.zip"), path("*_2_fastqc.zip"), path("*_1_fastqc.html"), path("*_2_fastqc.html")

    script:
    """
    echo "Processing sample: ${id}"
    echo "path to read 1: ${reads[0]}"
    echo "path to read 2: ${reads[1]}"
    fastqc --outdir . ${reads[0]} ${reads[1]}

    """
}

process qcAfter {
    publishDir "${params.outdir}/qc/after_trimming/fastqc", mode: 'copy', overwrite: true
    input: 
    tuple val(id), 
    path(trimmedReads_1),
    path(trimmedReads_2) 

    output:
    tuple path("*_1_trimmed_fastqc.zip"), path("*_2_trimmed_fastqc.zip"), path("*_1_trimmed_fastqc.html"), path("*_2_trimmed_fastqc.html")


    script:
    """
    echo "Processing sample: ${id}"
    echo "path to read 1: ${trimmedReads_1}"
    echo "path to read 2: ${trimmedReads_2}"
    fastqc --outdir . ${trimmedReads_1} ${trimmedReads_2}
    """
}  

process multiqcBefore {
    publishDir "${params.outdir}/qc/before_trimming", mode: 'copy', overwrite: true
    input:
    val allQcBefore


    output:
    path "multiqc_report.html"
    path "multiqc_data.zip"

    script:
    """
    multiqc \\
    ${allQcBefore.join(' ')} \\
    --outdir ./ \\
    --zip-data-dir \\
    --force || true 

    if [ ! -f multiqc_report.html ]; then
        echo "MultiQC failed to produce a report!" >&2
        exit 1
    fi
    """
}

process multiqcAfter {
    publishDir "${params.outdir}/qc/final", mode: 'copy', overwrite: true
    input:
    val(allQcAfter)
    val(allFlagstat)
    val(allCoverage)
    val(allBQSR)
    path(gatkEvalGATK)
    path(gatkEvalBCF)


    output:
    path "multiqc_report.html"
    path "multiqc_data.zip"

    script:
    """
    multiqc \\
    ${allQcAfter.join(' ')} \\
    ${allFlagstat.join(' ')} \\
    ${allCoverage.join(' ')} \\
    ${allBQSR.join(' ')} \\
    ${gatkEvalGATK} \\
    ${gatkEvalBCF} \\
    --outdir ./ \\
    --zip-data-dir \\
    --force || true 

    if [ ! -f multiqc_report.html ]; then
        echo "MultiQC failed to produce a report!" >&2
        exit 1
    fi
    """
}

// trimming reads
process trimReads {
    //publishDir "${params.outdir}/trimmed_reads", mode: 'copy', overwrite: true
    cpus = params.trimCPUs

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), 
    path("${id}_1_trimmed.fastq.gz"), 
    path("${id}_2_trimmed.fastq.gz") 

    script:
    """
    echo "Processing sample: ${id}"
    
    trimmomatic PE -threads ${task.cpus} \
    ${reads[0]} ${reads[1]} \
    ${id}_1_trimmed.fastq.gz \
    ${id}_1_untrimmed.fastq.gz \
    ${id}_2_trimmed.fastq.gz \
    ${id}_2_untrimmed.fastq.gz \
    ILLUMINACLIP:${params.adapters}:2:30:10 \
    SLIDINGWINDOW:4:25 \
    LEADING:25 \
    TRAILING:25 \
    MINLEN:50
    """
}

// optional process for indexing the reference genome
process indexReference{
    input:
    path reference

    output:
    path("${reference}.*")



    script:
    """
    echo "Indexing reference genome: ${reference}"
    # bwa index ${reference}
    """
}

// bwa mapping
process bwaMapping {
    cpus = params.bwaCPUs
    memory = params.bwaMemory// for wsl its minimum and maximum memory
    input:
    tuple val(id),
    path(trimmedReads_1),
    path(trimmedReads_2)

    output: 
    tuple val(id),
    path("${id}.bam")

    script:
        """
        echo "Mapping sample: ${id}"
        echo "Reference: ${params.reference}"
        echo "Trimmed Read 1: ${trimmedReads_1}"
        echo "Trimmed Read 2: ${trimmedReads_2}"

        echo "Reading readgroups"
        header=\$(zcat "${trimmedReads_1}" | head -n 1 )
        sample_id=\$(echo "\$header" | awk '{print \$1}' | tr -d '@' | cut -d '.' -f1)
        rest=\$(echo "\$header" | awk '{print \$2}')
        flowcell_barcode=\$(echo "\$rest" | cut -d ':' -f3)
        lane=\$(echo "\$rest" | cut -d ':' -f4)
        rg="@RG\\tID:\$flowcell_barcode.\$lane\\tSM:\$sample_id\\tPL:ILLUMINA\\tLB:\$sample_id\\tPU:\$flowcell_barcode.\$lane.\$sample_id"
        echo "Read group: \$rg"

        bwa mem -t ${task.cpus} -R "\$rg" \\
        ${params.reference} \\
        ${trimmedReads_1} \\
        ${trimmedReads_2} | samtools view -bS - > ${id}.bam
        """

}

// post aligment steps
process sortByQueryName {
    input:
    tuple val(id), path(mappedBam)

    output:
    tuple val(id), path("${mappedBam}.sorted.bam")

    """
    echo "Sorting BAM file by query name for sample: ${id}"
    samtools sort -n ${mappedBam} -o ${mappedBam}.sorted.bam
    """
}

process fixMate {
    input:
    tuple val(id), path(sortedBam)

    output:
    tuple val(id), path("${sortedBam}.fixmate.bam")

    script:
    """
    echo "Fixing mate information for sample: ${id}"
    samtools fixmate -mc ${sortedBam} ${sortedBam}.fixmate.bam
    """
}

process sortByCoordinate {
    input:
    tuple val(id), path(fixmateBam)

    output:
    tuple val(id), path("${fixmateBam}.sorted.bam")

    script:
    """
    echo "Sorting BAM file by coordinate for sample: ${id}"
    samtools sort ${fixmateBam} -o ${fixmateBam}.sorted.bam
    """
}

process markDuplicates {
    input:
    tuple val(id), path(sortedBam)

    output:
    tuple val(id), path("${id}.dedup.bam")

    script:
    """
    echo "Marking duplicates for sample: ${id}"
    samtools markdup -r ${sortedBam} ${id}.dedup.bam 

    # optional to add --json for statistics file that can be imported to multiqc 
    """
}

process indexBam {
    input:
    tuple val(id), path(dedupBam)

    output:
    tuple val(id), path("${dedupBam}_index.bai")

    script:
    """
    echo "Indexing BAM file for sample: ${id}"
    samtools index ${dedupBam} ${dedupBam}_index.bai
    """
}

// stats for deduplicated bam files
process flagstat {
    publishDir "${params.outdir}/qc/flagstat", mode: 'copy', overwrite: true
    input:
    tuple val(id), path(dedupBam)

    output:
    path ("${id}_flagstat.txt")


    script:
    """
    echo "Generating flagstat for sample: ${id}"
    samtools flagstat ${dedupBam} > ${id}_flagstat.txt
    """
}

process depth {
    publishDir "${params.outdir}/qc/depth", mode: 'copy', overwrite: true
    input:
    tuple val(id), path(dedupBam)

    output:
    path ("${id}_depth.txt")


    script:
    """
    echo "Generating depth for sample: ${id}"
    samtools depth  ${dedupBam} > ${id}_depth.txt # flag -a is for all positions (is correct but Cretes unnecessary large files) we have to calculkate mean from that 
    """

}

process coverage {
    publishDir "${params.outdir}/qc/coverage", mode: 'copy', overwrite: true
    input:
    tuple val(id), path(dedupBam)

    output:
    path ("${id}_coverage.txt")


    script:
    """
    echo "Calculating coverage for sample: ${id}"
    samtools coverage ${dedupBam} -o ${id}_coverage.txt
    """
}

// variant calling with bcftools
process bcfCalling {
    input:
    tuple val(id), path(dedupBam)
    tuple path(reference), path(refDict), path(refFai)
   

    output:
    tuple val(id), path("${id}.bcf")

    script:
    """
    echo "Calling variants using bcftools for sample: ${id}"

    bcftools mpileup -Ou -f ${reference} ${dedupBam} \\
    | bcftools call -mv -Ou -o ${id}.bcf
    """
}

process bcfFiltering {
    input:
    tuple val(id), path(bcf)

    output:
    tuple path("${id}.filtered.bcf"), path("${id}.filtered.csi")


    script:
    """
    echo "Filtering variants for sample: ${id}"
    bcftools filter -i 'QUAL >= 0 && DP >= 0' ${bcf} -Ou -o ${id}.filtered.bcf
    # i for include e for exclude 
    bcftools index -c ${bcf} -o ${id}.filtered.csi
    """
}
// not used indexing in filtering process
process indexBCF {
    input:
    path(bcf)

    output:
    path("${bcf}.csi")

    script:
    """
    bcftools index -c ${bcf}
    """
}

process bcfMerge {
    publishDir "${params.outdir}/vcf", mode: 'copy', overwrite: true
    input:
    val(bcfFiles)
    val(bcfIndexFiles)

    output:
    tuple path("cohort_bcftools.vcf.gz"), path("cohort_bcftools.vcf.gz.tbi")

    script:
    """
    bcftools merge \\
    ${bcfFiles.join(' ')} \\
    -O z \\
    -o cohort_bcftools.vcf.gz

    gatk IndexFeatureFile -I cohort_bcftools.vcf.gz
    """
}

// GATK 
process prepareReferenceGATK {
    input:
    path reference

    output:
    tuple path(reference),
    path("${reference.baseName}.dict"),
    path("${reference}.fai")

    script:
    """
    echo "reference : ${reference.baseName}"
    echo "Preparing reference genome for GATK: ${reference}"
    samtools faidx ${reference}
    gatk CreateSequenceDictionary -R ${reference} 
    """
}

process gatkBQSR {
    input:
    tuple val(id), path(dedupBam)
    tuple path(reference), path(refDict), path(refFai)
    path(knownSites)

    output: 
    tuple val(id), path("${id}.bqsr.bam"), path("${id}.bqsr.bai"), path ("${id}_covariates.pdf"), path("${id}_pre.table"), path("${id}_after.table")



    script:
    """
    echo "Running GATK Base Quality Score Recalibration for sample: ${id}"
    echo "Reference: ${reference}"
    echo "reference dict: ${refDict}"
    echo "Reference fai: ${refFai}"
    echo "Known sites: ${knownSites}"
    echo "Input BAM: ${dedupBam}"
    echo "Output BAM: ${id}.bqsr.bam"
    echo "Output BAI: ${id}.bqsr.bai"

    gatk BaseRecalibrator \\
    -I ${dedupBam} \\
    -R ${reference} \\
    --known-sites ${params.knownSites} \\
    -O ${id}_pre.table

    gatk ApplyBQSR \\
    -I ${dedupBam} \\
    -R ${reference} \\
    --bqsr-recal-file ${id}_pre.table \\
    -O ${id}.bqsr.bam

    gatk BaseRecalibrator \\
    -I ${id}.bqsr.bam \\
    -R ${reference} \\
    --known-sites ${params. knownSites} \\
    -O ${id}_after.table

    gatk AnalyzeCovariates \\
    -before ${id}_pre.table \\
    -after ${id}_after.table \\
    -plots ${id}_covariates.pdf \\
    """

}

process gatkHaplotypeCaller {
    cpus = params.gatkCPUs
    input: 
    tuple val(id), path(recalibratedBam), path(bqsrBai), path(covariatesPdf), path(pre_table), path(after_table)
    tuple path(reference), path(refDict), path(refFai)

    output: 
    tuple val(id), path("${id}.g.vcf.gz")


    script:
    """
    echo "Running GATK HaplotypeCaller for sample: ${id}"
    gatk HaplotypeCaller \\
    -R ${reference} \\
    -I ${recalibratedBam} \\
    -O ${id}.g.vcf.gz \\
    -ERC GVCF \\
    --native-pair-hmm-threads ${task.cpus} \\

    """
}

process indexGVCF {
    tag "$id"
    input:
    tuple val(id), path(gvcf)

    output:
    tuple val(id), path("${gvcf.getName()}"), path("${gvcf.getName()}.tbi")

    script:
    """
    echo "Indexing GVCF for sample: ${id}"
    gatk IndexFeatureFile -I ${gvcf}
    """
}

process combineGVCF {
    input: 
    tuple path(reference), path(refDict), path(refFai)
    path gvcfs 
    path tbis

    output: 
    path("cohort.g.vcf.gz")
    path("cohort.g.vcf.gz.tbi")
    script:
    """
    gatk CombineGVCFs \\
    -R ${reference} \\
    -O cohort.g.vcf.gz \\
    ${gvcfs.collect { "-V ${it}" }.join(" ")} 

    gatk IndexFeatureFile -I cohort.g.vcf.gz


    """
}

process gatkGenotype {
    publishDir "${params.outdir}/vcf", mode: 'copy', overwrite: true
    input: 
    tuple path(reference), path(refDict), path(refFai)
    path cohortGVCF
    path cohortTBIs

    output:
    tuple path("cohort_gatk.vcf.gz"), path("cohort_gatk.vcf.gz.tbi")


    script:
        """
        gatk GenotypeGVCFs \\
        -R ${reference} \\
        -V cohort.g.vcf.gz \\
        -O cohort_gatk.vcf.gz

        gatk IndexFeatureFile -I cohort_gatk.vcf.gz
        """
}

// vcf evaluation  
process gatkEvalGATK {
    publishDir "${params.outdir}/qc/vcf_eval", mode: 'copy', overwrite: true
    input:
    tuple path(vcf), path(vcfTbi)
    tuple path(reference), path(refDict), path(refFai)

    output:
    path("cohort_GATK_eval.txt")

    script:
    """
    echo "Evaluating VCF file: ${vcf}"
    gatk VariantEval -eval ${vcf} -R ${reference} -O cohort_GATK_eval.txt
    """
}

process gatkEvalBCF {
    publishDir "${params.outdir}/qc/vcf_eval", mode: 'copy', overwrite: true
    input:
    tuple path(vcf), path(vcfTbi)
    tuple path(reference), path(refDict), path(refFai)

    output:
    path("cohort_BCF_eval.txt")

    script:
    """
    echo "Evaluating VCF file: ${vcf}"
    gatk VariantEval -eval ${vcf} -R ${reference} -O cohort_BCF_eval.txt
    """
}
// optional bcf stats
process bcfStats {
    publishDir "${params.outdir}/qc/vcf_eval", mode: 'copy', overwrite: true
    input:
    tuple path(vcf), path(vcfTbi)

    output:
    path("cohort_bcf_stats.txt")

    script:
    """
    ls
    echo "Generating BCF stats for VCF file: ${vcf}"
    bcftools stats ${vcf} > cohort_bcf_stats.txt
    """
}



// saving results
process saveBQSRPDF {
    publishDir "${params.outdir}/qc/bqsr_covariates", mode: 'copy', overwrite: true
    input: 
    tuple val(id), path("${id}.bqsr.bam"), path("${id}.bqsr.bai"), path ("${id}_covariates.pdf"), path("${id}_pre.table"), path("${id}_after.table")

    output:
    path("${id}_covariates.pdf")

    script:
    """
    echo "Saving BQSR files for sample: ${id}"
    ls ${id}_covariates.pdf
    """
}

process saveBCFVCF {
    input:
    tuple val(id), path(vcf)

    output:
    stdout
    script:
    """
    mkdir -p ${params.outdir}/vcf/bcf
    echo "Saving VCF for sample: ${id}"
    cp ${vcf} ${params.outdir}/vcf/bcf/${id}.filtered.vcf.gz
    """
}

process saveMappedReads {
    input:
    tuple val(id), path(mappedBam)

    output:
    stdout
    script:
    """
    mkdir -p ${params.outdir}/mapped_reads
    echo "Saving mapped reads for sample: ${id}"
    cp ${mappedBam} ${params.outdir}/mapped_reads/${id}.bam
    """
}

process saveReads {
    input:
    tuple val(id), 
    path(dedupBam)

    output:
    stdout
    script:
    """
    mkdir -p ${params.outdir}/post_aligment/reads
    echo "Saving reads for sample: ${id}"
    cp ${dedupBam} ${params.outdir}/post_aligment/reads/${id}_dedup.bam
    """
}

process saveBQSR {
    input:
    tuple val(id), path(bqsrBam), path(bqsrBai), path(covariates)

    output:
    stdout
    script:
    """
    mkdir -p ${params.outdir}/post_aligment/bqsr
    echo "Saving BQSR files for sample: ${id}"
    cp ${bqsrBam} ${params.outdir}/post_aligment/bqsr/${id}.bqsr.bam
    cp ${bqsrBai} ${params.outdir}/post_aligment/bqsr/${id}.bqsr.bai
    cp ${covariates} ${params.outdir}/post_aligment/bqsr/${id}_covariates.pdf
    """
}

process saveGATKVCF {
    input:
    tuple val(id), path(gatkVcf)

    output:
    stdout
    script:
    """
    mkdir -p ${params.outdir}/vcf/gatk_vcf
    echo "Saving GATK VCF for sample: ${id}"
    cp ${gatkVcf} ${params.outdir}/vcf/gatk_vcf/${id}.g.vcf.gz
    """
}

process saveGenotypedVCF {
    input:
    path genotypedVCF

    output:
    stdout
    script:
    """
    mkdir -p ${params.outdir}/vcf
    cp ${genotypedVCF} ${params.outdir}/vcf/cohort_genotyped.vcf.gz
    """
}


workflow {
    //referenceChannel = Channel.fromPath(params.reference, checkIfExists: true).view()

    readsChannel             = Channel.fromFilePairs(params.reads, checkIfExists: true).view()
    referenceGATKChannel     = prepareReferenceGATK(params.reference)
    knownSitesChannel        = Channel.fromPath(params.knownSites, checkIfExists: true).view()

    //QC before trimming
    qcBeforeChannel   = qcBefore(readsChannel).collect().set { allQcBefore }
    multiqcBefore(allQcBefore)

    //trimming
    trimmedReadsChannel = trimReads(readsChannel).view()

    //fastQC after trimming 
    qcAfterChannel      = qcAfter(trimmedReadsChannel).collect().set { allQcAfter }

    // bwa mapping
    mappedReadsChannel  = bwaMapping(trimmedReadsChannel).view()

    // post aligment steps
    sortedByQueryNameChannel  = sortByQueryName(mappedReadsChannel)
    fixmateChannel            = fixMate(sortedByQueryNameChannel)
    sortedByCoordinateChannel = sortByCoordinate(fixmateChannel)
    dedupChannel              = markDuplicates(sortedByCoordinateChannel)
    indexedBamChannel         = indexBam(dedupChannel)

    // stats for deduplicated bam files
    flagstatChannel = flagstat(dedupChannel).collect().set { allFlagstat }
    depthChannel    = depth(dedupChannel)
    coverageChannel = coverage(dedupChannel).collect().set { allCoverage }


    // variant calling bcf 
    bcfCalledChannel   = bcfCalling(dedupChannel,referenceGATKChannel)
    filteredBcfChannel = bcfFiltering(bcfCalledChannel)
    filteredBcfChannel.map { it[0] }.collect().set{ allBcf }
    filteredBcfChannel.map { it[1] }.collect().set{ allBcfIndex }

    mergedBcfChannel   = bcfMerge(allBcf, allBcfIndex)
    
    // GATK BQSR
    gatkBQSRChannel = gatkBQSR(dedupChannel, referenceGATKChannel, params.knownSites).view() 
    multiqcBQSRChannel = gatkBQSRChannel.map { id, bqsrBam, bqsrBai, covariatesPdf, preTable, afterTable ->
                                                    tuple(covariatesPdf, preTable, afterTable)
                                                        }.collect().set { allBQSR }

    


    // GATK variant calling - haplotype caller
    gatkCalledChannel = gatkHaplotypeCaller(gatkBQSRChannel, referenceGATKChannel).view()
    // index gvcf
    indexedGVCFChannel = indexGVCF(gatkCalledChannel)
    // collect all gvcfs and index files for them wait for all of them to finish
    // collect gvcf files
    indexedGVCFChannel.map { it [1] }.collect().set{ allGVCFs }
    // collect index files 
    indexedGVCFChannel.map { it [2] }.collect().set{ allTBIs }

    //indexedGVCFChannel.map { [it[1], it[2]] }.flatten().collect().set { allGVCFs }

    cohortChannel = combineGVCF(referenceGATKChannel, allGVCFs, allTBIs)

    genotypedChannel = gatkGenotype(referenceGATKChannel, cohortChannel).view()

    // GATK evaluation
    gatkEvalChannel = gatkEvalGATK(genotypedChannel, referenceGATKChannel)
    bcfEvalChannel  = gatkEvalBCF(mergedBcfChannel, referenceGATKChannel)

    // optional bcf stats
    //bcfStatsChannelBCF = bcfStats(mergedBcfChannel)
    //bcfStatsChannelGATK = bcfStats(genotypedChannel)


    // downstream multiqc
    multiqcAfter(allQcAfter, allFlagstat, allCoverage, allBQSR, gatkEvalChannel, bcfEvalChannel)

    //save results
    saveBQSRPDF(gatkBQSRChannel)
}

