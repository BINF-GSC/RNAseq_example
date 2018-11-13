#!/usr/bin/env nextflow

// params
params.fastqs = "/media/maxh/storage/RNAseq_example/raw_fastqs"
params.out_dir = "/media/maxh/storage/RNAseq_example/out"
params.trim_cpus = 3
params.fastqc_cpus = 3

fastq_pairs = Channel
    .fromFilePairs( "$params.fastqs/*{.read1,.read2}.fastq.gz" )
    .ifEmpty { exit 1, "Fastq file(s) not found at: $params.fastqs" }

process trim {
    container "binfgsc/trimmomatic:latest"
    cpus params.trim_cpus

    input:
    set pair_id, file(reads) from fastq_pairs

    output:
    set file("*.log"), file("*.sum") into trim_log
    set file("*P.fastq.gz"), file("*U.fastq.gz") into trimmed

    """
    java -jar \$TRIM PE \
        -threads $params.trim_cpus \
        -trimlog ${pair_id}.log \
        -summary ${pair_id}.sum \
        $reads \
        -baseout ${pair_id}.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process fastqc {
    container 'binfgsc/fastqc:latest'
    cpus params.fastqc_cpus
    publishDir "$params.out_dir/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set file(paired), file(unpaired) from trimmed

    output:
    file "*" into fastqc

    script:
    """
    fastqc -t $params.fastqc_cpus --no-extract $paired $unpaired
    """
}

process multiqc {
    container 'binfgsc/multiqc:latest'
    publishDir "$params.out_dir/multiqc", mode: 'copy'

    input:
    file fastqc from fastqc.collect()
    file trim_log from trim_log.collect()

    output:
    file "multiqc_report.html"

    script:
    """
    multiqc $fastqc $trim_log -m fastqc -m trimmomatic
    """
}