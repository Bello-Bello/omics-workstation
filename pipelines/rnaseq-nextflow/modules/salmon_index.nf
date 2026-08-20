process SALMON_INDEX {
    conda "${projectDir}/envs/quant.yaml"
    container "quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5"

    input:
    path transcriptome_fasta

    output:
    path "salmon_index", emit: index

    script:
    """
    salmon index -t ${transcriptome_fasta} -i salmon_index -k 31 -p ${task.cpus}
    """
}
