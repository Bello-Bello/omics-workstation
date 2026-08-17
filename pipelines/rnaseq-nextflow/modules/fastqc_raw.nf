process FASTQC_RAW {
    tag "$sample"
    conda "${projectDir}/envs/qc.yaml"
    publishDir "${params.outdir}/fastqc/raw", mode: 'copy'

    input:
    tuple val(sample), path(fq1), path(fq2)

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip

    script:
    """
    fastqc --outdir . ${fq1} ${fq2}
    """
}
