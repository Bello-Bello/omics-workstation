process SALMON_INDEX {
    conda "${projectDir}/envs/quant.yaml"

    input:
    path transcriptome_fasta

    output:
    path "salmon_index", emit: index

    script:
    """
    salmon index -t ${transcriptome_fasta} -i salmon_index -k 31 -p ${task.cpus}
    """
}
