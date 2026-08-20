process STARSOLO {
    tag "$sample"
    conda "${projectDir}/envs/star.yaml"
    container "quay.io/biocontainers/star:2.7.11b--h5ca1c30_8"
    publishDir "${params.outdir}/starsolo", mode: 'copy'

    input:
    tuple val(sample), path(r1_files), path(r2_files)
    path index
    path whitelist

    output:
    // --outFileNamePrefix below means every output file/dir STAR writes is
    // actually named "${sample}.<whatever STAR normally calls it>" — easy
    // to miss, since the paths below don't look like it at first glance.
    tuple val(sample), path("${sample}.Solo.out/Gene/raw"), emit: raw_matrix
    path "${sample}.Log.final.out", emit: log

    script:
    // 10x reads come as two files per lane: R1 = 16bp cell barcode + 10bp
    // UMI (chemistry v2), R2 = the actual cDNA read. STARsolo wants the
    // cDNA read listed FIRST — the reverse of what the file names suggest.
    // Multiple lanes of the same sample are just comma-joined, no manual
    // concatenation step needed.
    def r1_list = r1_files instanceof List ? r1_files.join(",") : r1_files
    def r2_list = r2_files instanceof List ? r2_files.join(",") : r2_files
    """
    STAR \
        --runMode alignReads \
        --genomeDir ${index} \
        --readFilesIn ${r2_list} ${r1_list} \
        --readFilesCommand zcat \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist ${whitelist} \
        --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 \
        --soloFeatures Gene \
        --outSAMtype None \
        --outFileNamePrefix ${sample}. \
        --runThreadN ${task.cpus}
    """
}
