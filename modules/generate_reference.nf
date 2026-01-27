process GENERATE_REFERENCE {
    tag "all_samples"
    publishDir "${params.outdir}/protein_ref", mode: 'copy'
    
    input:
    path filtered_gtf  // 使用过滤后的GTF文件作为输入
    path fasta         // 参考基因组FASTA文件
    val projectDir
    
    output:
    path "ORF_no_redundancy.fa", emit: orf_fa
    
    script:
    def filtered_gtf_file = filtered_gtf.toString()
    def fasta_file = fasta.toString()
    """
    #从转录本预测ORF
    python ${projectDir}/bin/proteinAnalysis.py generateReference \
        -g ${filtered_gtf_file} \
        -f ${fasta_file} \
        -t ${task.cpus}
    
    # 重命名输出文件
    if [ -f "ORF_no_redundancy.fa" ]; then
        echo "ORF_no_redundancy.fa found"
    else
        # 尝试寻找可能的替代输出文件名
        mv \$(dirname ${filtered_gtf})/ORF_no_redundancy.fa . || true
    fi
    """
    
    stub:
    """
    echo ">protein1" > ORF_no_redundancy.fa
    echo "MAHLEQ..." >> ORF_no_redundancy.fa
    """
}

