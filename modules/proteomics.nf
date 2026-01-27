// 添加诱饵数据库的进程
process ADD_DECOY {
    tag "add_decoy"
    publishDir "${params.outdir}/protein_ref", mode: 'copy'
    
    input:
    path ref_fa
    
    output:
    path "*-decoys-*.fas", emit: decoy_db
    
    script:
    def work_dir = "decoy_workspace"
    """
    mkdir -p ${work_dir}
    cd ${work_dir}
    cp ../${ref_fa} ./
    ${params.philosopher_path} workspace --init --nocheck --temp ..
    ${params.philosopher_path} database --custom ./${ref_fa}
    ${params.philosopher_path} workspace --clean --nocheck
    
    # 找到生成的数据库文件并复制到上级目录
    decoy_file=\$(find . -name "*-decoys-*.fas")
    cp \$decoy_file ../
    cd ..
    """
}

// 运行FragPipe的进程
process FRAGPIPE {
    tag "fragpipe"
    publishDir "${params.outdir}/protein", mode: 'copy'
    
    cpus 10
    memory '32 GB'
    time '48 h'
    
    input:
    path decoy_db
    path manifest
    path workflow_config
    val philosopher_path
    val fragpipe_path
    val msfragger_jar
    val ionquant_jar
    
    output:
    path "fragpipe_output", emit: fragpipe_output
    
    script:
    // def decoy_db_path = decoy_db.toAbsolutePath().toString()
    def work_dir = "fragpipe_workspace"
    """
    # 创建工作目录
    mkdir -p ${work_dir}
    cd ${work_dir}
    
    # 复制并修改配置文件
    cp ../${workflow_config} workflow.workflow
    cp ../${manifest} manifest.fp-manifest
    DECOY_DB_PATH=\$(readlink -f ../${decoy_db})
    sed -i "s|^database.db-path=.*|database.db-path=\${DECOY_DB_PATH}|" workflow.workflow
    
    # 运行FragPipe
    #export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
    ${fragpipe_path} --headless \\
        --workflow workflow.workflow \\
        --manifest ${manifest} \\
        --workdir . \\
        --config-msfragger ${msfragger_jar} \\
        --config-ionquant ${ionquant_jar} \\
        --config-philosopher ${philosopher_path} \\
        --threads ${task.cpus} > fragpipe.log 2>&1
    
    # 整理输出
    cd ..
    mv ${work_dir} fragpipe_output
    """
}
