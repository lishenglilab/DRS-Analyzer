process CREATE_CLEAN_MANIFEST {
    tag "create_manifest"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    tuple val(sample), val(condition), val(batch), path(raw_fq)
    tuple val(sample), path(clean_fq)
    
    output:
    path "clean_manifest.txt", emit: clean_manifest
    
    script:
    """
    # 创建新的manifest文件头
    echo -e "Sample\\tCondition\\tBatch\\tClean_FASTQ" > clean_manifest.txt
    
    # 添加每个样本的信息
    while IFS= read -r line; do
        sample_name=\$(echo "\$line" | cut -f1)
        condition=\$(echo "\$line" | cut -f2)
        batch=\$(echo "\$line" | cut -f3)
        clean_path=\$(find . -name "\${sample_name}.clean.fastq.gz" -print -quit)
        echo -e "\${sample_name}\\t\${condition}\\t\${batch}\\t\${clean_path}" >> clean_manifest.txt
    done < "${params.manifest}"
    """
    
    stub:
    """
    echo -e "Sample\\tCondition\\tBatch\\tClean_FASTQ" > clean_manifest.txt
    echo -e "sample1\\tcondition1\\tbatch1\\t/path/to/sample1.clean.fastq.gz" >> clean_manifest.txt
    echo -e "sample2\\tcondition2\\tbatch2\\t/path/to/sample2.clean.fastq.gz" >> clean_manifest.txt
    """
}

