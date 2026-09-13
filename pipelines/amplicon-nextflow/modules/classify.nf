// Taxonomic assignment for each ASV, by naive Bayes over k-mer profiles.
process CLASSIFY_TAXONOMY {
    container params.qiime_container
    publishDir "${params.outdir}/qiime", mode: 'copy'

    input:
    path rep_seqs
    path classifier

    output:
    path "taxonomy.qza", emit: taxonomy
    path "taxonomy.qzv", emit: viz

    script:
    """
    qiime feature-classifier classify-sklearn \\
        --i-classifier ${classifier} \\
        --i-reads ${rep_seqs} \\
        --p-n-jobs ${task.cpus} \\
        --o-classification taxonomy.qza

    qiime metadata tabulate \\
        --m-input-file taxonomy.qza \\
        --o-visualization taxonomy.qzv
    """
}
