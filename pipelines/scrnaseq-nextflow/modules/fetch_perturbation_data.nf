process FETCH_PERTURBATION_DATA {
    tag "$dataset_name"
    conda "${projectDir}/envs/scanpy.yaml"
    publishDir "${params.outdir}/perturbation/raw", mode: 'copy'

    input:
    val dataset_name

    output:
    path "${dataset_name}_raw.h5ad", emit: raw_h5ad

    script:
    """
    fetch_perturbation_data.py ${dataset_name} ${dataset_name}_raw.h5ad
    """
}
