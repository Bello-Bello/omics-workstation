nextflow.enable.dsl = 2

// Separate entrypoint from main.nf (STARsolo quantification) rather than
// one combined workflow: this stage starts from a public, already-quantified
// count matrix, not FASTQ. Real drug-screen FASTQ sets are typically
// hundreds of GB (sci-Plex sequences per-well, not per-cell) — reprocessing
// reprocessing them from raw reads is impractical at this scale, and most
// perturbation analysis in practice starts from a delivered count matrix
// rather than re-running alignment per project. See README.md for the full
// reasoning.

include { FETCH_PERTURBATION_DATA } from './modules/fetch_perturbation_data'
include { PERTURBATION_ANALYSIS }   from './modules/perturbation_analysis'

workflow {
    FETCH_PERTURBATION_DATA(params.dataset)
    PERTURBATION_ANALYSIS(FETCH_PERTURBATION_DATA.out.raw_h5ad)
}
