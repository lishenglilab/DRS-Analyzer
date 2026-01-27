process REMOVE_HOMO {
    tag "all_samples"
    
    publishDir "${params.outdir}/protein_ref", mode: 'copy', overwrite: true, enable: true

    input:
    path orf_fa
    path ref
    val projectDir
    
    output:
    path "ORF_ref.fa", emit: filtered_orf
    
    script:
    """
    echo "DEBUG: params.outdir = ${params.outdir}" >&2
    echo "DEBUG: Publishing to: ${params.outdir}/protein_ref" >&2
    # 运行同源蛋白过滤
    python ${projectDir}/bin/proteinAnalysis.py removeHomo \
        -r ${orf_fa} \
        -f ${ref} \
        -t ${task.cpus}
    
    # 确保输出文件存在
    if [ ! -f "ORF_unmatch_ref.fa" ]; then
        # 尝试从可能的位置复制
        mv \$(dirname ${orf_fa})/ORF_unmatch_ref.fa . || true
    fi
    
    cat ORF_unmatch_ref.fa ${ref} > ORF_ref.fa
    """
    
    stub:
    """
    cp ${orf_fa} ORF_unmatch_ref.fa
    """
}
