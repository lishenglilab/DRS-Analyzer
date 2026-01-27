process CREATE_MAPPING_FILE {
    tag "$sample"
    
    input:
    tuple val(sample), val(mapping_reads)
    
    output:
    path("${sample}_mapping.txt")
    
    script:
    """
    echo '${sample} ${mapping_reads}' > ${sample}_mapping.txt
    """
}

process CREATE_SIZE_FILE {
    tag "$sample"
    
    input:
    tuple val(sample), val(size_gb)
    
    output:
    path("${sample}_size.txt")
    
    script:
    """
    echo '${sample} ${size_gb}' > ${sample}_size.txt
    """
}

process COMBINE_QC_STATS {
    publishDir "${params.outdir}/table", mode: 'copy'
    
    input:
    path(length_files)
    path(mapping_files)
    path(size_files)
    path(fail_reads_file)
    
    output:
    path("QC.xls")
    
    script:
    """
    # All files are already staged in the work directory
    # Length files: *_Length_pass.txt
    # Mapping files: *_mapping.txt
    # Size files: *_size.txt
    
    # List files for debugging
    echo "=== Files in work directory ===" >&2
    ls -la >&2
    echo "================================" >&2
    
    # Handle fail reads file
    if [ -f "${fail_reads_file}" ] && [ -s "${fail_reads_file}" ]; then
        cp ${fail_reads_file} read_number_fail.txt
    else
        touch read_number_fail.txt
    fi
    
    # Run the R script with Rscript (doesn't require execute permission)
    Rscript ${projectDir}/bin/combine_qc_stats.R
    """
    
    stub:
    """
    echo -e "sample_id\\tSize_GB\\tTotal_reads\\tpass_percent\\tMedian_pass_read_length\\tMean_pass_read_length\\tMapped_read\\tMapped_ratio" > QC.xls
    echo -e "test_sample\\t10.5\\t1000000\\t90.5\\t800\\t950\\t850000\\t85.0" >> QC.xls
    """
}
