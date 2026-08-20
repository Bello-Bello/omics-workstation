nextflow.enable.dsl = 2

include { STAR_INDEX } from './modules/star_index'
include { STARSOLO }   from './modules/starsolo'

workflow {
    // One row per (sample, lane) pair — Sample_Y has two rows (two lanes),
    // Sample_X has one. .groupTuple() below is what merges same-sample rows
    // together; this raw channel is still one item per row.
    reads_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample, file(row.fq1), file(row.fq2)) }

    // Groups by the first tuple element (sample name), collecting every
    // matching row's fq1/fq2 into lists — Sample_X ends up with 1-item
    // lists, Sample_Y with 2-item lists (one per lane). New operator, not
    // needed in the bulk pipeline since every sample there was single-lane.
    grouped_reads_ch = reads_ch.groupTuple()

    STAR_INDEX(file(params.genome_fasta), file(params.gtf))
    STARSOLO(grouped_reads_ch, STAR_INDEX.out.index, file(params.barcode_whitelist))
}
