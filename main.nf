#!/usr/bin/env nextflow

Channel
    .fromFilePairs( "$params.fastqs/*{.read1,.read2}.fastq.gz" )
    .ifEmpty { exit 1, "Fastq file(s) not found at: $params.fastqs" }
    .into { fq_pairs_fastqc; fq_pairs_trim }

process fastqc_raw {
    input:
    set pairId, file(fq_pair) from fq_pairs_fastqc

    output:
    file "*" into fastqc_raw

    script:
    """
    fastqc -t $task.cpus --no-extract $fq_pair
    """
}

process trim {
    input:
    set pair_id, file(fq_pair) from fq_pairs_trim

    output:
    set file("*.log"), file("*.sum") into trim_log
    set pair_id, file("*P.fastq.gz"), file("*U.fastq.gz") into trimmed_fastqc, trimmed_align

    """
    java -jar \$TRIM PE \
        -threads $task.cpus \
        -summary ${pair_id}.sum \
        $fq_pair \
        -baseout ${pair_id}.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 MINLEN:36 2> ${pair_id}.log
    """
}

process fastqc_trimmed {
    input:
    set pair_id, file(paired), file(unpaired) from trimmed_fastqc

    output:
    file "*" into fastqc

    script:
    """
    fastqc -t $task.cpus --no-extract $paired $unpaired
    """
}

Channel
    .fromPath( "$params.genome/*.fa.gz" )
    .ifEmpty { exit 1, "Fasta file not found at: $params.genome" }
    .set { genome }
Channel
    .fromPath( "$params.gtf/*.gtf.gz")
    .ifEmpty { exit 1, "GTF file not found at: $params.gtf" }
    .into { gtf_star_index; gtf_star_align }

process starIndex {
    input:
    file genome
    file gtf from gtf_star_index

    output:
    file "star" into star_index
    
    script:
    """
    mkdir star
    star \\
        --runMode genomeGenerate \\
        --runThreadN $task.cpus \\
        --genomeDir star/ \\
        --sjdbGTFfile $gtf \\
        --genomeFastaFiles $genome \\
        --limitGenomeGenerateRAM ${task.memory.toBytes()}
    """
}

process starAlign {
    publishDir "$params.out_dir/star", mode: 'copy'

    input:
    set pair_id, file(paired), file(unpaired) from trimmed_align
    file index from star_index
    file gtf from gtf_star_align

    output:
    set file("*Log.final.out"), file('*.bam') into star_aligned
    file "*.out" into alignment_logs
    file "*SJ.out.tab" into star_sj
    file "*Log.out" into star_log

    script:
    """
    star \\
        --genomeDir $index \\
        --sjdbGTFfile $gtf \\
        --readFilesIn $paired \\
        --runThreadN $tast.cpus \\
        --outWigType bedGraph \\
        --outSameType BAM SortedByCoordinate \\
        --readFilesCommand gunzip -c \\
        --outFileNamePrefix $pair_id
    """

}

process multiqc {
    publishDir "$params.out_dir/multiqc", mode: 'copy'

    input:
    file fastqc from fastqc.collect()
    file trim_log from trim_log.collect()
    file alignment_logs.collect()

    output:
    file "multiqc_report.html"

    script:
    """
    multiqc $fastqc $trim_log $alignment_logs
    """
}
