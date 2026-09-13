nextflow.enable.dsl = 2

// Amplicon (16S) pipeline: raw multiplexed reads -> ASVs -> taxonomy ->
// compositional statistics, on real public soil data.
//
// Two halves, deliberately separated at EXPORT_TABLE:
//   1. QIIME2, in a pinned container, doing what QIIME2 is good at —
//      demultiplexing, DADA2 denoising, taxonomic classification.
//   2. Ordinary Python doing the statistics, in a small conda env, with
//      the compositional methods implemented from first principles in
//      compositional/ rather than called from a black box.
//
// See README.md for why the second half is not `qiime composition ancombc`.

include { FETCH_DATA; FETCH_CLASSIFIER } from './modules/fetch_data'
include { IMPORT_EMP; DEMUX }            from './modules/import_and_demux'
include { DENOISE_DADA2 }                from './modules/denoise'
include { CLASSIFY_TAXONOMY }            from './modules/classify'
include { EXPORT_TABLE }                 from './modules/export_table'
include { COMPOSITIONAL_ANALYSIS }       from './modules/compositional_analysis'

workflow {
    FETCH_DATA()
    FETCH_CLASSIFIER()

    IMPORT_EMP(FETCH_DATA.out.seq_dir)
    DEMUX(IMPORT_EMP.out, FETCH_DATA.out.metadata)
    DENOISE_DADA2(DEMUX.out.demux)
    CLASSIFY_TAXONOMY(DENOISE_DADA2.out.rep_seqs, FETCH_CLASSIFIER.out)
    EXPORT_TABLE(DENOISE_DADA2.out.table, CLASSIFY_TAXONOMY.out.taxonomy)
    COMPOSITIONAL_ANALYSIS(
        EXPORT_TABLE.out.table_tsv,
        EXPORT_TABLE.out.taxonomy_tsv,
        FETCH_DATA.out.metadata,
    )
}
