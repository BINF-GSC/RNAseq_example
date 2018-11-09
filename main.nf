#!/usr/bin/env nextflow

fastqs = Channel
    .fromPath( params.fastq )
    .ifEmpty { exit 1, "Fastq file(s) not found at: ${params.fastq}" }
    .collect()

process fastqc_raw {
    cpus 3
    publishDir "/media/maxh/storage/rnaseq_example/fastqc_out", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    file fastqs

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc -t 3 --no-extract $fastqs
    """
}