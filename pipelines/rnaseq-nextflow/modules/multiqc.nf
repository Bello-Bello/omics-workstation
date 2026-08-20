process MULTIQC {
    conda "${projectDir}/envs/qc.yaml"
    container "quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path fastqc_files
    path fastp_files
    path salmon_dirs

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc . --force
    """
}
