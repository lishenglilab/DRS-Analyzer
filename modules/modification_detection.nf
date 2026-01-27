process MODIFICATION_DETECTION {
  tag "$sample"
  publishDir "${params.outdir}/modifications/$sample", mode: 'copy'
  
  cpus 16
  memory '32 GB'
  time '480 h'
  
  input:
  tuple val(sample), path(fast5_dir), path(clean_fq)
  path ref_tr_fa  // 来自FLAIR的转录本参考
  val projectDir  // 项目目录路径
  
  output:
  path "m6a/${sample}/m6A_output.bed", emit: m6a_results
  path "m5c/${sample}/m5c.prediction.txt", emit: m5c_results
  path "psu/${sample}/alignment/prediction.csv", emit: psu_results  // 需确认nanopsu实际输出名
  
  script:
  def modDir = "$projectDir/bin/mod"
  def cpu1 = (task.cpus / 2).toInteger() - 1
  def cpu2 = task.cpus - 1
  """
  # 0. 创建工作目录
  mkdir -p singleFast5/${sample} m6a/${sample} m5c/${sample} psu/${sample}
  
  # 1. FAST5转换
  multi_to_single_fast5 \\
    -i "$fast5_dir" \\
    -s "singleFast5/${sample}" \\
    -t ${task.cpus}
  
  # 2. Tombo重定信号
  tombo resquiggle \\
    --ignore-read-locks \\
    --overwrite "singleFast5/${sample}" \\
    "$ref_tr_fa" \\
    --processes ${task.cpus} \\
    --num-most-common-errors 5
  
  # 3. 检测修饰
  # 3.1 m6A检测
  tombo detect_modifications de_novo \\
    --fast5-basedirs "singleFast5/${sample}" \\
    --statistics-file-basename m6a/${sample}/m6a \\
    --processes ${task.cpus}
    
  tombo text_output browser_files \\
    --fast5-basedirs "singleFast5/${sample}" \\
    --statistics-filename m6a/${sample}/m6a.tombo.stats \\
    --genome-fasta "${ref_tr_fa}" \\
    --file-types coverage fraction \\
    --browser-file-basename m6a/m6a_result
    
  wig2bed < m6a/m6a_result.fraction_modified_reads.plus.wig > m6a/frac.bed
    
  python "${modDir}/cDNA_MINES.py" \\
    --ref "${ref_tr_fa}" \\
    --coverage m6a/m6a_result.coverage.plus.bedgraph \\
    --fraction_modified m6a/frac.bed \\
    --output m6a/${sample}/m6A_output.bed \\
    --kmer_models "$modDir/Final_Models/names.txt"
  
  # 3.2 m5C检测
  tombo detect_modifications alternative_model \\
    --fast5-basedirs "singleFast5/${sample}" \\
    --statistics-file-basename m5c/${sample}/m5c \\
    --alternate-bases 5mC \\
    --processes ${task.cpus}
    
  tombo text_output browser_files \\
    --fast5-basedirs "singleFast5/${sample}" \\
    --statistics-filename m5c/${sample}/m5c.5mC.tombo.stats \\
    --genome-fasta "${ref_tr_fa}" \\
    --file-types coverage statistic dampened_fraction fraction \\
    --browser-file-basename m5c/${sample}/m5c
    
  python "${modDir}/m5c_res_process.py" \\
    --out_dir m5c/${sample} \\
    --plus_file m5c/${sample}/m5c.fraction_modified_reads.plus.wig \\
    --minus_file m5c/${sample}/m5c.fraction_modified_reads.minus.wig \\
    --ref_tr_fa "${ref_tr_fa}"
  
  # 3.3 psU检测
  cd psu/${sample}
  ln -s "../../${clean_fq}" "./${sample}.clean.fastq.gz"
  ln -s "../../${ref_tr_fa}" ./
  pigz -dc "${sample}.clean.fastq.gz" > pass.fastq
  nanopsu alignment -i . -r ./\$(basename "${ref_tr_fa}")
  nanopsu remove_intron
  nanopsu extract_features
  nanopsu prediction
  """
}

