#!/usr/bin/env nextflow

// Define the working directory
params.workDir = '/home/paq/NGS/NGS'
params.outdir  = '/home/paq/NGS/NGS/outdir'

// Define input files (paired-end reads)

params.reference = "${params.workDir}/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

params.adapters = "${params.workDir}/adapters/combined_adapters.fa"
params.reads = "${params.workDir}/data/*_{1,2}.fastq.gz"

params.trimCPUs = 2
params.bwaCPUs = 4
params.bwaMemory = '7 GB '

// Create a channel for paired-end reads



process qcBefore {
    input: 
    tuple val(id), path(reads) // reads is a list of two files

    output:
    stdout

    script:
    """
    mkdir -p ${params.outdir}/qc/before_trimming
    echo "Processing sample: ${id}"
    echo "path to read 1: ${reads[0]}"
    echo "path to read 2: ${reads[1]}"
    fastqc --outdir ${params.outdir}/qc/before_trimming ${reads[0]} ${reads[1]}

    """
}

process trimReads {
    cpus = params.trimCPUs

    input:
    tuple val(id), path(reads) // reads is a list of two files

    output:
    tuple val(id), 
    path("${id}_1.trimmed.fastq.gz"), 
    path("${id}_2.trimmed.fastq.gz") // Output files in the working directory

    script:
    """
    mkdir -p trimmed_reads
    echo "Processing sample: ${id}"
    
    trimmomatic PE -threads ${task.cpus} \
    ${reads[0]} ${reads[1]} \
    ${id}_1.trimmed.fastq.gz \
    ${id}_1.untrimmed.fastq.gz \
    ${id}_2.trimmed.fastq.gz \
    ${id}_2.untrimmed.fastq.gz \
    ILLUMINACLIP:${params.adapters}:2:30:10 \
    SLIDINGWINDOW:4:25 \
    LEADING:25 \
    TRAILING:25 \
    MINLEN:50
    """
}

process qcAfter {
    input: 
    tuple val(id), 
    path(trimmedReads_1),
    path(trimmedReads_2) 

    output:
    stdout

    script:
    """
    mkdir -p ${params.outdir}/qc/after_trimming
    echo "Processing sample: ${id}"
    echo "path to read 1: ${trimmedReads_1}"
    echo "path to read 2: ${trimmedReads_2}"
    fastqc --outdir ${params.outdir}/qc/after_trimming ${trimmedReads_1} ${trimmedReads_2}
    """
}  
// optional process for indexing the reference genome
process indexReference{
    input:
    path reference

    output:
    path("${reference}.*") // Declare all index files as output



    script:
    """
    echo "Indexing reference genome: ${reference}"
    # bwa index ${reference}
    """
}


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


process postAligment { // optional process with all post aligment steps in one process
    input:
    tuple val(id), path(mappedBam)

    output:
    tuple val(id),
    path("${mappedBam}.dedup.bam"),
    path("${mappedBam}_index.bai")

    script:
    """
    mkdir -p ${params.outdir}/post_aligment
    echo "Processing sample: ${id}"
    echo "Mapped BAM file: ${mappedBam}"
    # sort by query name 
    samtools sort -n ${mappedBam} -o ${mappedBam}.sorted.bam
    # fixmate
    samtools fixmate -mc ${mappedBam}.sorted.bam ${mappedBam}.fixmate.bam
    # sort by coordinate
    samtools sort ${mappedBam}.fixmate.bam -o ${mappedBam}.sorted.bam
    # mark duplicates
    samtools markdup -r ${mappedBam}.sorted.bam ${mappedBam}.dedup.bam
    # index 
    samtools index ${mappedBam}.dedup.bam ${mappedBam}_index.bai
    # flagstat
    samtools flagstat ${mappedBam}.dedup.bam > ${params.outdir}/post_aligment/${id}_flagstat.txt
    # depth 
    samtools depth -a ${mappedBam}.dedup.bam -o ${params.outdir}/post_aligment/${id}_depth.txt
    """


}


process sortByQueryName {
    input:
    tuple val(id), path(mappedBam)

    output:
    tuple val(id), path("${mappedBam}.sorted.bam")

    script:
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
    tuple val(id), path("${sortedBam}.dedup.bam")

    script:
    """
    echo "Marking duplicates for sample: ${id}"
    samtools markdup -r ${sortedBam} ${sortedBam}.dedup.bam
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

process flagstat {
    input:
    tuple val(id), path(dedupBam)

    output:
    stdout

    script:
    """
    mkdir -p ${params.outdir}/post_aligment
    echo "Generating flagstat for sample: ${id}"
    samtools flagstat ${dedupBam} > ${params.outdir}/post_aligment/${id}_flagstat.txt
    """
}

process depth {
    input:
    tuple val(id), path(dedupBam)

    output:
    stdout

    script:
    """
    mkdir -p ${params.outdir}/post_aligment
    echo "Generating depth for sample: ${id}"
    samtools depth  ${dedupBam} > ${params.outdir}/post_aligment/${id}_depth.txt # flag -a is for all positions (is correct but Cretes unnecessary large files)
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


workflow {
    readsChannel = Channel.fromFilePairs(params.reads, checkIfExists: true).view()
    //referenceChannel = Channel.fromPath(params.reference, checkIfExists: true).view()


    //qc before trimming
    qcBefore(readsChannel)


    //trimming
    trimmedReadsChannel = trimReads(readsChannel).view()

    //QC after trimming 
    qcAfter(trimmedReadsChannel).view()

    mappedReadsChannel = bwaMapping(trimmedReadsChannel).view()

    // post aligment steps
    sortedByQueryNameChannel = sortByQueryName(mappedReadsChannel)
    fixmateChannel = fixMate(sortedByQueryNameChannel)
    sortedByCoordinateChannel = sortByCoordinate(fixmateChannel)
    dedupChannel = markDuplicates(sortedByCoordinateChannel)
    indexedBamChannel = indexBam(dedupChannel)
    flagstat(dedupChannel)
    depth(dedupChannel)
    
    //postAligmentChannel = postAligment(mappedReadsChannel).view() // optional process with all post aligment steps in one process
    saveMappedReads(mappedReadsChannel).view()
    saveReads(dedupChannel).view()


}

