// Leave the QIIME2 artifact world: BIOM -> TSV, taxonomy -> TSV.
//
// Done as an explicit step rather than inside the analysis script so the
// boundary is visible in the DAG. Everything upstream of here is QIIME2;
// everything downstream is ordinary Python that could just as well read a
// table from a different platform.
process EXPORT_TABLE {
    container params.qiime_container
    publishDir "${params.outdir}/exported", mode: 'copy'

    input:
    path table
    path taxonomy

    output:
    path "feature-table.tsv", emit: table_tsv
    path "taxonomy.tsv",      emit: taxonomy_tsv

    script:
    """
    qiime tools export --input-path ${table} --output-path table_export
    biom convert \\
        -i table_export/feature-table.biom \\
        -o feature-table.tsv \\
        --to-tsv

    qiime tools export --input-path ${taxonomy} --output-path taxonomy_export
    mv taxonomy_export/taxonomy.tsv taxonomy.tsv

    # An empty table here means denoising dropped everything; catching it
    # now gives a readable error instead of a pandas traceback.
    n_features=\$(tail -n +3 feature-table.tsv | wc -l)
    echo "Exported \$n_features features"
    test "\$n_features" -gt 0
    """
}
