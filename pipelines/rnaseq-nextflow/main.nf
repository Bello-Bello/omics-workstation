nextflow.enable.dsl = 2

include { FASTQC_RAW }   from './modules/fastqc_raw'
include { FASTP }        from './modules/fastp'
include { SALMON_INDEX } from './modules/salmon_index'
include { SALMON_QUANT } from './modules/salmon_quant'
include { MULTIQC }      from './modules/multiqc'
include { DESEQ2 }       from './modules/deseq2'

workflow {
    reads_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample, file(row.fq1), file(row.fq2)) }

    FASTQC_RAW(reads_ch)
    FASTP(reads_ch)

    SALMON_INDEX(file(params.transcriptome_fasta))
    SALMON_QUANT(FASTP.out.trimmed, SALMON_INDEX.out.index)

    MULTIQC(
        FASTQC_RAW.out.zip.collect(),
        FASTP.out.json.collect(),
        SALMON_QUANT.out.quant_dir.map { sample, dir -> dir }.collect()
    )

    DESEQ2(
        SALMON_QUANT.out.quant_dir.map { sample, dir -> dir }.collect(),
        file(params.samples),
        file(params.gtf)
    )
}
