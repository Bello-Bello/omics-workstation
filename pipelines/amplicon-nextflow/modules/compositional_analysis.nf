// CLR -> Aitchison PCA / clustering -> PERMANOVA -> differential abundance.
//
// The only process in this pipeline that does not use the QIIME2 container:
// it runs in the small conda env at envs/compositional.yaml, because
// nothing here needs QIIME2 and pinning the analysis to QIIME2's release
// cadence would be a needless constraint. See that file for the reasoning.
process COMPOSITIONAL_ANALYSIS {
    conda "${projectDir}/envs/compositional.yaml"
    publishDir "${params.outdir}/compositional", mode: 'copy'

    input:
    path table_tsv
    path taxonomy_tsv
    path metadata

    output:
    path "*.csv",       emit: tables
    path "*.png",       emit: figures
    path "summary.json", emit: summary
    path "permanova.json", emit: permanova
    path "report.md",   emit: report

    script:
    // --project-root points at this pipeline's checkout so `import
    // compositional` resolves — same pattern as the perturbation pipeline.
    """
    run_compositional_analysis.py \\
        --table ${table_tsv} \\
        --taxonomy ${taxonomy_tsv} \\
        --metadata ${metadata} \\
        --outdir . \\
        --project-root ${projectDir} \\
        --group-column ${params.group_column} \\
        --group-a ${params.group_a} \\
        --group-b ${params.group_b} \\
        --rank ${params.rank} \\
        --min-prevalence ${params.min_prevalence} \\
        --min-total-count ${params.min_total_count} \\
        --min-depth ${params.min_depth} \\
        --permutations ${params.permutations} \\
        --monte-carlo ${params.monte_carlo} \\
        --alpha ${params.alpha} \\
        --seed ${params.seed}
    """
}
