process FASTP {
    tag "$sample"
    conda "${projectDir}/envs/qc.yaml"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample), path(fq1), path(fq2)

    output:
    tuple val(sample), path("${sample}_R1.trim.fastq.gz"), path("${sample}_R2.trim.fastq.gz"), emit: trimmed
    path "${sample}.fastp.json", emit: json
    path "${sample}.fastp.html", emit: html

    script:
    """
    fastp -i ${fq1} -I ${fq2} \
        -o ${sample}_R1.trim.fastq.gz -O ${sample}_R2.trim.fastq.gz \
        -j ${sample}.fastp.json -h ${sample}.fastp.html \
        --thread ${task.cpus}
    """
}
