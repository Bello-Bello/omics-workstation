process SALMON_QUANT {
    tag "$sample"
    conda "${projectDir}/envs/quant.yaml"
    publishDir "${params.outdir}/salmon", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)
    path index

    output:
    tuple val(sample), path("${sample}"), emit: quant_dir

    script:
    """
    salmon quant -i ${index} -l A -1 ${r1} -2 ${r2} -p ${task.cpus} --validateMappings -o ${sample}
    """
}
