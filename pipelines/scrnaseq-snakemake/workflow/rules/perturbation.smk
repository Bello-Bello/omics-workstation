rule fetch_perturbation_data:
    output:
        "resources/perturbation/{dataset}_raw.h5ad",
    params:
        dataset=lambda wc: wc.dataset,
    conda:
        "../../envs/scanpy.yaml"
    log:
        "logs/fetch_perturbation_data_{dataset}.log",
    script:
        "../scripts/fetch_perturbation_data.py"


rule perturbation_analysis:
    input:
        "resources/perturbation/{}_raw.h5ad".format(config["dataset"]),
    output:
        processed_h5ad="results/perturbation/processed.h5ad",
        de_csv="results/perturbation/differential_expression.csv",
        summary="results/perturbation/summary.json",
        umap_plot="results/perturbation/umap_condition_cluster.png",
    params:
        **config["perturbation"],
    conda:
        "../../envs/scanpy.yaml"
    log:
        "logs/perturbation_analysis.log",
    script:
        "../scripts/run_perturbation_pipeline.py"
