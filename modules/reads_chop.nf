process READS_CHOP {
    tag "$sample"
    publishDir "${params.outdir}/clean_fq", mode: 'copy'
    
    input:
    tuple val(sample), path(raw_fq)
    
    output:
    tuple val(sample), path("${sample}.clean.fastq.gz"), emit: clean_fq
    
    script:
    """
    porechop -i ${raw_fq} -o ${sample}.clean.fastq.gz -t ${task.cpus}
    """
    
    stub:
    """
    touch ${sample}.clean.fastq.gz
    """
}
