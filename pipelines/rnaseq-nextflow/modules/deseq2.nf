process DESEQ2 {
    conda "${projectDir}/envs/deseq2.yaml"
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    path quant_dirs
    path samples_tsv
    path gtf

    output:
    path "results.csv"
    path "pca.png"
    path "counts_normalized.csv"

    script:
    """
    deseq2.R ${samples_tsv} ${gtf}
    """
}
