rule fastp:
    input:
        r1=lambda wc: samples.loc[wc.sample, "fq1"],
        r2=lambda wc: samples.loc[wc.sample, "fq2"],
    output:
        r1="results/trimmed/{sample}_R1.trim.fastq.gz",
        r2="results/trimmed/{sample}_R2.trim.fastq.gz",
        json="results/trimmed/{sample}.fastp.json",
        html="results/trimmed/{sample}.fastp.html",
    params:
        extra=config["params"]["fastp"]["extra"],
    threads: config["threads"]["fastp"]
    conda:
        "../envs/qc.yaml"
    container:
        "docker://quay.io/biocontainers/fastp:0.23.4--h125f33a_5"
    log:
        "logs/fastp/{sample}.log"
    shell:
        "fastp -i {input.r1} -I {input.r2} "
        "-o {output.r1} -O {output.r2} "
        "-j {output.json} -h {output.html} "
        "--thread {threads} {params.extra} > {log} 2>&1"
