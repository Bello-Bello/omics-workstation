process STAR_INDEX {
    conda "${projectDir}/envs/star.yaml"
    container "quay.io/biocontainers/star:2.7.11b--h5ca1c30_8"

    input:
    path genome_fasta
    path gtf

    output:
    path "star_index", emit: index

    script:
    // --genomeSAindexNbases scaled down from STAR's default of 14 —
    // appropriate for a small single-chromosome genome (~61 Mb); the
    // default is tuned for whole genomes and STAR will warn if it's left
    // too high here.
    """
    mkdir star_index
    STAR \
        --runMode genomeGenerate \
        --genomeDir star_index \
        --genomeFastaFiles ${genome_fasta} \
        --sjdbGTFfile ${gtf} \
        --sjdbOverhang 90 \
        --genomeSAindexNbases 11 \
        --runThreadN ${task.cpus}
    """
}
