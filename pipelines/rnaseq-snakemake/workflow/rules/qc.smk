rule fastqc_raw:
    input:
        r1=lambda wc: samples.loc[wc.sample, "fq1"],
        r2=lambda wc: samples.loc[wc.sample, "fq2"],
    output:
        html1="results/fastqc/raw/{sample}_R1_fastqc.html",
        html2="results/fastqc/raw/{sample}_R2_fastqc.html",
    conda:
        "../envs/qc.yaml"
    log:
        "logs/fastqc/raw/{sample}.log"
    shell:
        "fastqc --outdir results/fastqc/raw {input.r1} {input.r2} > {log} 2>&1"
