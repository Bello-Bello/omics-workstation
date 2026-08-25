nextflow.enable.dsl = 2

// Separate entrypoint from main.nf (STARsolo quantification) rather than
// one combined workflow: this stage starts from a public, already-quantified
// count matrix, not FASTQ. Real drug-screen FASTQ sets are typically
// hundreds of GB (sci-Plex sequences per-well, not per-cell) — reprocessing
// them from raw reads isn't practical for a portfolio pipeline, and most
// real perturbation-analysis work in industry starts from a delivered
// count matrix anyway, not from re-running alignment per project. See
// README.md for the full reasoning.

include { FETCH_PERTURBATION_DATA } from './modules/fetch_perturbation_data'
include { PERTURBATION_ANALYSIS }   from './modules/perturbation_analysis'

workflow {
    FETCH_PERTURBATION_DATA(params.dataset)
    PERTURBATION_ANALYSIS(FETCH_PERTURBATION_DATA.out.raw_h5ad)
}
