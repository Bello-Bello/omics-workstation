// DADA2: error-model denoising into amplicon sequence variants.
//
// ASVs, not 97%-identity OTUs. The distinction matters for everything
// downstream: DADA2 fits a per-run error model, so it can tell a real
// single-nucleotide variant from a sequencing error and resolve exact
// sequences. OTU clustering cannot, and papers a 3% divergence band over
// the difference. ASVs are also comparable across studies without
// re-clustering, which OTUs are not.
process DENOISE_DADA2 {
    container params.qiime_container
    publishDir "${params.outdir}/qiime", mode: 'copy'

    input:
    path demux

    output:
    path "table.qza",            emit: table
    path "rep-seqs.qza",         emit: rep_seqs
    path "denoising-stats.qza",  emit: stats
    path "denoising-stats.tsv",  emit: stats_tsv

    script:
    """
    qiime dada2 denoise-paired \\
        --i-demultiplexed-seqs ${demux} \\
        --p-trim-left-f ${params.trim_left_f} \\
        --p-trim-left-r ${params.trim_left_r} \\
        --p-trunc-len-f ${params.trunc_len_f} \\
        --p-trunc-len-r ${params.trunc_len_r} \\
        --p-n-threads ${task.cpus} \\
        --o-table table.qza \\
        --o-representative-sequences rep-seqs.qza \\
        --o-denoising-stats denoising-stats.qza

    # Export the stats immediately: read attrition through filter/merge/
    # chimera removal is the single most useful diagnostic when a run
    # produces a suspiciously empty table, and it should be in the results
    # directory whether or not the rest of the pipeline succeeds.
    qiime tools export --input-path denoising-stats.qza --output-path stats_export
    mv stats_export/stats.tsv denoising-stats.tsv
    """
}
