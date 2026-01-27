process BAM_FLAGSTAT {
    tag "$sample"
    publishDir "${params.qc_dir}/bam_map", mode: 'copy'
    
    input:
    tuple val(sample), path(bam)
    
    output:
    tuple val(sample), path("${sample}.map.stat"), emit: flagstat
    tuple val(sample), env(mapping_reads), emit: stats
    
    script:
    """
    samtools flagstat ${bam} > ${sample}.map.stat
    export mapping_reads=\$(grep "mapped (" ${sample}.map.stat | head -n 1 | awk '{print \$1}')
    echo "${sample} \$mapping_reads" > ${sample}_mapping.txt
    """
    
    stub:
    """
    touch ${sample}.map.stat
    export mapping_reads="1000000"
    echo "${sample} \$mapping_reads" > ${sample}_mapping.txt
    """
}
