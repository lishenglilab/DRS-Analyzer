process GET_FILE_SIZE {
    tag "$sample"
    
    input:
    tuple val(sample), path(fastq)
    
    output:
    tuple val(sample), env(size_gb), emit: sizes
    
    script:
    """
    size_bytes=\$(stat -L -c%s "${fastq}")
    export size_gb=\$(echo "\$size_bytes 1073741824" | awk '{printf "%.2f", \$1/\$2}')
    echo "${sample} \${size_gb}" > ${sample}_size.txt
    """
    
    stub:
    """
    export size_gb="10.5"
    echo "${sample} \$size_gb" > ${sample}_size.txt
    """
}
