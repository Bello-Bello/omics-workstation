rule multiqc:
    input:
        expand("results/fastqc/raw/{sample}_R1_fastqc.html", sample=sample_names()),
        expand("results/trimmed/{sample}.fastp.json", sample=sample_names()),
        expand("results/salmon/{sample}/quant.sf", sample=sample_names()),
    output:
        "results/multiqc/multiqc_report.html",
    conda:
        "../envs/qc.yaml"
    container:
        "docker://quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0"
    log:
        "logs/multiqc.log"
    shell:
        "multiqc results/fastqc results/trimmed results/salmon "
        "--outdir results/multiqc --force > {log} 2>&1"
