rule salmon_index:
    input:
        fasta=config["ref"]["transcriptome_fasta"],
    output:
        index=directory("resources/salmon_index"),
    params:
        extra=config["params"]["salmon_index"]["extra"],
    threads: config["threads"]["salmon_index"]
    conda:
        "../envs/quant.yaml"
    container:
        "docker://quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5"
    log:
        "logs/salmon/index.log"
    shell:
        "salmon index -t {input.fasta} -i {output.index} "
        "-p {threads} {params.extra} > {log} 2>&1"

rule salmon_quant:
    input:
        r1="results/trimmed/{sample}_R1.trim.fastq.gz",
        r2="results/trimmed/{sample}_R2.trim.fastq.gz",
        index="resources/salmon_index",
    output:
        quant="results/salmon/{sample}/quant.sf",
    params:
        libtype=config["params"]["salmon_quant"]["libtype"],
        outdir=lambda wc: f"results/salmon/{wc.sample}",
    threads: config["threads"]["salmon_quant"]
    conda:
        "../envs/quant.yaml"
    container:
        "docker://quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5"
    log:
        "logs/salmon/quant/{sample}.log"
    shell:
        "salmon quant -i {input.index} -l {params.libtype} "
        "-1 {input.r1} -2 {input.r2} -p {threads} "
        "--validateMappings -o {params.outdir} > {log} 2>&1"
