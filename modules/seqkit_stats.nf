process SEQKIT_LENGTH {
    tag "$sample"
    publishDir "${params.qc_dir}/pass", mode: 'copy'
    
    input:
    tuple val(sample), path(fastq)
    
    output:
    tuple val(sample), path("${sample}_Length_pass.txt"), emit: lengths
    
    script:
    """
    seqkit fx2tab -j ${task.cpus} -l -n -i -H ${fastq} > ${sample}_Length_pass.txt
    """
    
    stub:
    """
    echo -e "name\\tlength" > ${sample}_Length_pass.txt
    for i in {1..1000}; do
        echo -e "read_\$i\\t\$((RANDOM % 2000 + 500))" >> ${sample}_Length_pass.txt
    done
    """
}