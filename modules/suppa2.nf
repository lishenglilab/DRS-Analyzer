process SUPPA2_ANALYSIS {
    tag "SUPPA2"
    publishDir "${params.outdir}/suppa", mode: 'copy'
    
    input:
    path gtf_file  // FLAIR_PIP输出的过滤GTF
    path tpm_file  // FLAIR_PIP输出的定量TPM文件
    
    output:
    path "*.ioe", emit: ioe_files
    path "*.psi", emit: psi_files
    
    script:
    def out_prefix = "transcriptEvents"
    def tpm_suppa = "transcript_expression.tsv"
    """
    # 转换TPM文件格式
    cat ${tpm_file} | sed 's/ID\t//g' > ${tpm_suppa}
    
    # 生成可变剪切事件
    suppa.py generateEvents -i ${gtf_file} -o ${out_prefix} -f ioe -e SE SS MX RI FL
    
    # 合并所有事件文件
    awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' ${out_prefix}*.ioe > ${out_prefix}.all.events.ioe
    
    # 计算PSI值
    suppa.py psiPerEvent --ioe-file ${out_prefix}.all.events.ioe -e ${tpm_suppa} -o ${out_prefix}.all.events.psi
    """
    
    stub:
    """
    touch transcriptEvents.SE.ioe
    touch transcriptEvents.SS.ioe
    touch transcriptEvents.MX.ioe
    touch transcriptEvents.RI.ioe
    touch transcriptEvents.FL.ioe
    touch transcriptEvents.all.events.ioe
    touch transcriptEvents.all.events.psi
    """
}

