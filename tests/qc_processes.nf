nextflow.enable.dsl = 2

include { SEQKIT_LENGTH } from '../modules/seqkit_stats.nf'
include { BAM_FLAGSTAT } from '../modules/bam_stats.nf'

params.fastq = null
params.bam = null
params.qc_dir = null
params.outdir = null

// Serialize only the declared emissions, including the work-directory output paths.
workflow LENGTHS {
    SEQKIT_LENGTH(Channel.of(tuple('toy', file(params.fastq, checkIfExists: true))))
    SEQKIT_LENGTH.out.lengths.collectFile(name: 'lengths.tsv', storeDir: params.outdir, newLine: true) { sample, output ->
        "${sample}\t${output}"
    }
}

workflow FLAGSTAT {
    BAM_FLAGSTAT(Channel.of(tuple('toy', file(params.bam, checkIfExists: true))))
    BAM_FLAGSTAT.out.flagstat.collectFile(name: 'flagstat.tsv', storeDir: params.outdir, newLine: true) { sample, output ->
        "${sample}\t${output}"
    }
    BAM_FLAGSTAT.out.stats.collectFile(name: 'stats.tsv', storeDir: params.outdir, newLine: true) { sample, count ->
        "${sample}\t${count}"
    }
}
