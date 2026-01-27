process GET_FILE_SIZE {
    tag "$sample"
    
    input:
    tuple val(sample), path(fastq)
    
    output:
    tuple val(sample), env(size_gb), emit: sizes
    
    script:
    """
    size=\$(stat -c%s ${fastq})
    export size_gb=\$(echo "scale=2; \$size / 1024 / 1024 / 1024" | bc)
    echo "${sample} \$size_gb" > ${sample}_size.txt
    """
    
    stub:
    """
    export size_gb="10.5"
    echo "${sample} \$size_gb" > ${sample}_size.txt
    """
}
