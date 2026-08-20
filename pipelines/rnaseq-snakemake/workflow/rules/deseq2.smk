rule deseq2:
    input:
        quants=expand("results/salmon/{sample}/quant.sf", sample=sample_names()),
        samples=config["samples"],
        gtf=config["ref"]["gtf"],
    output:
        results="results/deseq2/results.csv",
        pca="results/deseq2/pca.png",
        counts="results/deseq2/counts_normalized.csv",
    conda:
        "../envs/deseq2.yaml"
    # Same custom image as the Nextflow pipeline (see
    # ../../../rnaseq-nextflow/docker/deseq2.Dockerfile) — build it locally
    # and reference it as "docker-daemon://omics-deseq2:1.0" for Apptainer to
    # use the local Docker cache instead of pulling from a registry.
    container:
        "docker://omics-deseq2:1.0"
    log:
        "logs/deseq2.log"
    script:
        "../scripts/deseq2.R"
