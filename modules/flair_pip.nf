process FLAIR_PIP {
    tag "all_samples"
    publishDir "${params.outdir}/flair", mode: 'copy'
    
    input:
    path manifest
    path gtf
    path fasta
    val project_dir  // 添加项目目录输入
    
    output:
    path "isoforms.fa", emit: isoforms_fa
    path "quant.tpm.tsv", emit: quant_tsv
    path "filtered_isoforms.gtf", emit: filtered_gtf
    path "align", emit: bam_files  // 新增：输出BAM文件目录
    path "flair_quantify_filter.tpm.tsv", emit: filtered_quant_tsv
    
    script:
    def manifest_file = manifest.toString()
    def fasta_file = fasta.toString()
    def gtf_file = gtf.toString()
    """
    # 创建临时目录结构
    mkdir -p align
    mkdir -p correct
    mkdir -p collapse
    mkdir -p quant
    
    # 为每个样本运行比对
    cat ${manifest_file} | while read -r sample condition batch fastq fast5; do
        mkdir -p "align/\${sample}"
        flair align -g ${fasta_file} -r "\${fastq}" -o "align/\${sample}/\${sample}" -t ${task.cpus}
    done
    
    # 合并所有BED文件
    cat align/*/*.bed > align/merge.bed
    
    # 运行校正步骤
    flair correct -q align/merge.bed -f ${gtf_file} -g ${fasta_file} --output correct/flair --threads ${task.cpus}
    
    # 获取所有FASTQ路径
    fastq_list=\$(awk -F'\\t' '{print \$4}' ${manifest_file} | tr '\\n' ',')
    fastq_list="\${fastq_list%,}"
    
    # 运行组装步骤
    flair collapse -g ${fasta_file} -q correct/flair_all_corrected.bed --gtf ${gtf_file} -r "\${fastq_list}" --output collapse/isoforms --threads ${task.cpus} --temp_dir collapse
    
    # 运行定量步骤
    flair quantify -r ${manifest_file} -i collapse/isoforms.isoforms.fa --output quant/flair_quantify --threads ${task.cpus} --temp_dir quant --sample_id_only --tpm
    
    # 运行过滤步骤
    Rscript ${project_dir}/bin/transcriptFilter.R -t quant/flair_quantify.tpm.tsv -a collapse/isoforms.isoforms.gtf -q ${params.tpm_cutoff} -o .
    
    # 重命名输出文件
    cp collapse/isoforms.isoforms.fa isoforms.fa
    cp quant/flair_quantify.tpm.tsv quant.tpm.tsv
    cp isoforms_filter.gtf filtered_isoforms.gtf
    """
    
    stub:
    """
    touch isoforms.fa
    touch quant.tpm.tsv
    touch filtered_isoforms.gtf
    mkdir -p align/sample1
    touch align/sample1/sample1.bam
    """
}

