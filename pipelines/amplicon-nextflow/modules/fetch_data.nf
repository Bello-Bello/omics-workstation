// Downloads the raw multiplexed reads and the sample metadata.
//
// Kept as its own process (rather than folded into the import step) so
// Nextflow's `-resume` caches the ~300 MB download across runs while you
// are still iterating on DADA2 parameters downstream.
process FETCH_DATA {
    tag "${params.subsample ?: 'full'}"
    container params.qiime_container
    storeDir "${projectDir}/.data_cache/${params.subsample ?: 'full'}"

    output:
    path "emp-paired-end-sequences", emit: seq_dir
    path "sample_metadata.tsv",      emit: metadata

    script:
    def sub = params.subsample ? "/${params.subsample}" : ""
    """
    mkdir -p emp-paired-end-sequences
    for f in forward reverse barcodes; do
        curl -sSL --retry 3 --fail \\
            -o emp-paired-end-sequences/\${f}.fastq.gz \\
            "${params.data_base}${sub}/\${f}.fastq.gz"
        gzip -t emp-paired-end-sequences/\${f}.fastq.gz
    done
    curl -sSL --retry 3 --fail -o sample_metadata.tsv "${params.metadata_url}"

    # Fail loudly here rather than letting QIIME2 report a confusing schema
    # error three steps later.
    test -s sample_metadata.tsv
    grep -q 'barcode-sequence' sample_metadata.tsv
    """
}

// The pretrained naive-Bayes classifier, cached the same way (it is ~30 MB
// and never changes for a pinned URL).
process FETCH_CLASSIFIER {
    container params.qiime_container
    storeDir "${projectDir}/.data_cache/classifier"

    output:
    path "classifier.qza"

    script:
    """
    curl -sSL --retry 3 --fail -o classifier.qza "${params.classifier_url}"
    # A .qza is a zip; an HTML error page saved under the same name is not.
    python -c "import zipfile,sys; sys.exit(0 if zipfile.is_zipfile('classifier.qza') else 1)"
    """
}
