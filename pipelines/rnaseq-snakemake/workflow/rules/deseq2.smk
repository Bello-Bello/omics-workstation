rule deseq2:
    input:
        quants=expand("results/salmon/{sample}/quant.sf", sample=sample_names()),
        samples=config["samples"],
        gtf=config["ref"]["gtf"],
    output:
        results="results/deseq2/results.csv",
        pca="results/deseq2/pca.png",
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2.log"
    script:
        "../scripts/deseq2.R"
