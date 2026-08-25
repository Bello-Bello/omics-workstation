process PERTURBATION_ANALYSIS {
    conda "${projectDir}/envs/scanpy.yaml"
    publishDir "${params.outdir}/perturbation", mode: 'copy'

    input:
    path raw_h5ad

    output:
    path "processed.h5ad", emit: processed_h5ad
    path "differential_expression.csv", emit: de_csv
    path "summary.json", emit: summary
    path "umap_condition_cluster.png", emit: umap_plot

    script:
    // --project-root points at this pipeline's own checkout so the process
    // can `import perturbation` (the local package at ../perturbation/) —
    // see the docstring in run_perturbation_pipeline.py for why this is
    // passed explicitly instead of inferred from the script's own path.
    """
    run_perturbation_pipeline.py \
        --input ${raw_h5ad} \
        --outdir . \
        --project-root ${projectDir} \
        --perturbation-col ${params.perturbation_col} \
        --dose-col ${params.dose_col} \
        --control-label ${params.control_label} \
        --mt-prefix ${params.mt_prefix} \
        --min-genes ${params.min_genes} \
        --max-pct-mt ${params.max_pct_mt} \
        --resolution ${params.leiden_resolution}
    """
}
