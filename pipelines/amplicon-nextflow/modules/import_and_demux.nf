// EMP-protocol multiplexed reads -> per-sample demultiplexed reads.
//
// The Earth Microbiome Project layout keeps the barcodes in a *third*
// FASTQ, parallel to forward and reverse, rather than inline in the read
// headers. That is why this cannot use the manifest-based import path.
process IMPORT_EMP {
    container params.qiime_container
    publishDir "${params.outdir}/qiime", mode: 'copy'

    input:
    path seq_dir

    output:
    path "emp-paired-end-sequences.qza"

    script:
    """
    qiime tools import \\
        --type EMPPairedEndSequences \\
        --input-path ${seq_dir} \\
        --output-path emp-paired-end-sequences.qza
    """
}

process DEMUX {
    container params.qiime_container
    publishDir "${params.outdir}/qiime", mode: 'copy'

    input:
    path emp_qza
    path metadata

    output:
    path "demux.qza",             emit: demux
    path "demux-details.qza",     emit: details
    path "demux.qzv",             emit: viz

    script:
    // --p-rev-comp-mapping-barcodes is required for this dataset: the
    // barcodes as sequenced are the reverse complement of the ones written
    // in the metadata sheet. Omitting it does not error — it demultiplexes
    // almost nothing, which is a much worse failure because it looks like
    // a data quality problem rather than a flag problem.
    """
    qiime demux emp-paired \\
        --m-barcodes-file ${metadata} \\
        --m-barcodes-column barcode-sequence \\
        --p-rev-comp-mapping-barcodes \\
        --i-seqs ${emp_qza} \\
        --o-per-sample-sequences demux.qza \\
        --o-error-correction-details demux-details.qza

    qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
    """
}
