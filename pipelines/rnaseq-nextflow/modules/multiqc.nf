process MULTIQC {
    conda "${projectDir}/envs/qc.yaml"
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
