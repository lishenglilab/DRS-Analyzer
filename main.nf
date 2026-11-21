#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

projectDir = file("$baseDir")

// Parameters
params.outdir = "/home/yvzeng/Project/mouseDRS/data/nextflow/testdata/result/"
params.qc_dir = "/home/yvzeng/Project/mouseDRS/data/nextflow/testdata/QC"
params.reference_fasta = "/home/public/reference/fa/mouse/GRCm39.primary_assembly.genomeV29.fa"
params.annotation_gtf = "/home/public/reference/gtf/mouse/gencode.vM29.primary_assembly.annotation.gtf"
params.ref_prot = "/home/yvzeng/Project/mouseDRS/data/refPep/uniprot_sprot_mouse.fasta"
params.manifest = "/home/yvzeng/Project/mouseDRS/data/nextflow/testdata/manifest.txt"
params.tpm_cutoff = 0.1
params.run_proteomics = false
params.proteomics_manifest = null
params.proteom_workflow = null
params.philosopher_path = "/path/to/MSPipe/philosopher/philosopher"
params.fragpipe_path = "/path/to/MSPipe/fragpipe/bin/fragpipe"
params.msfragger_jar = null
params.ionquant_jar = null
params.fragpipe_cli = 'legacy'
params.fragpipe_tools_dir = null
params.fragpipe_python = null
params.fragpipe_dry_run = false
params.fragpipe_ram_gb = null

// Module selection parameter
params.module = null  // Options: reads_chop, flair_pip, quality_report, reference_generate, modification, proteomics, suppa2, or null for full pipeline

// Module-specific input parameters
params.clean_fastq_dir = null  // For modules that need clean fastq
params.bam_dir = null  // For quality_report
params.filtered_gtf = null  // For reference_generate and suppa2
params.isoforms_fa = null  // For modification
params.quant_tsv = null  // For suppa2
params.protein_ref = null  // For proteomics

// Include processes
include { BAM_FLAGSTAT } from './modules/bam_stats.nf'
include { SEQKIT_LENGTH } from './modules/seqkit_stats.nf'
include { GET_FILE_SIZE } from './modules/file_size.nf'
include { CREATE_MAPPING_FILE; CREATE_SIZE_FILE; COMBINE_QC_STATS } from './modules/combine_stats.nf'
include { READS_CHOP } from './modules/reads_chop.nf'
include { FLAIR_PIP } from './modules/flair_pip.nf'
include { GENERATE_REFERENCE } from './modules/generate_reference.nf'
include { REMOVE_HOMO } from './modules/remove_homo.nf'
include { MODIFICATION_DETECTION } from './modules/modification_detection.nf'
include { ADD_DECOY } from './modules/proteomics.nf'
include { FRAGPIPE } from './modules/proteomics.nf'
include { SUPPA2_ANALYSIS } from './modules/suppa2.nf'

// Individual module workflows
workflow reads_chop_module {
    if (!params.manifest) {
        error "ERROR: --manifest is required for reads_chop module"
    }
    
    // Create output directories
    file("${params.outdir}/clean_fq").mkdirs()
    
    samples_raw_fq_ch = Channel
        .fromPath(params.manifest)
        .splitCsv(sep: "\t", header: false)
        .map { row -> tuple(row[0], file(row[3].trim())) }
    
    READS_CHOP(samples_raw_fq_ch)
    
    // Save outputs and create manifest
    READS_CHOP.out.clean_fq
        .collectFile(name: "clean_manifest.txt", storeDir: "${params.outdir}") { sample, fq ->
            "${sample}\t${fq}\n"
        }
    
    emit:
    clean_fq = READS_CHOP.out.clean_fq
}

workflow flair_pip_module {
    take:
    clean_manifest
    
    main:
    if (!params.annotation_gtf || !params.reference_fasta) {
        error "ERROR: --annotation_gtf and --reference_fasta are required for flair_pip module"
    }
    
    // Create output directories
    file("${params.outdir}/flair").mkdirs()
    
    FLAIR_PIP(
        clean_manifest, 
        file(params.annotation_gtf), 
        file(params.reference_fasta), 
        projectDir
    )
    
    // Save key outputs to specified locations
    FLAIR_PIP.out.filtered_gtf
        .collectFile(name: "filtered_isoforms.gtf", storeDir: "${params.outdir}/flair")
    
    FLAIR_PIP.out.quant_tsv
        .collectFile(name: "quant.tpm.tsv", storeDir: "${params.outdir}/flair")
    
    FLAIR_PIP.out.isoforms_fa
        .collectFile(name: "isoforms.fa", storeDir: "${params.outdir}/flair")
    
    FLAIR_PIP.out.filtered_quant_tsv
        .collectFile(name: "flair_quantify_filter.tpm.tsv", storeDir: "${params.outdir}/flair")
    
    emit:
    filtered_gtf = FLAIR_PIP.out.filtered_gtf
    quant_tsv = FLAIR_PIP.out.quant_tsv
    isoforms_fa = FLAIR_PIP.out.isoforms_fa
    bam_files = FLAIR_PIP.out.bam_files
}

workflow quality_report_module {
    take:
    bam_ch
    clean_fq_ch
    
    main:
    // Create output directories
    file("${params.qc_dir}/bam_map").mkdirs()
    file("${params.qc_dir}/pass").mkdirs()
    file("${params.outdir}/table").mkdirs()
    
    // Run the individual stat processes
    BAM_FLAGSTAT(bam_ch)
    SEQKIT_LENGTH(clean_fq_ch)
    GET_FILE_SIZE(clean_fq_ch)
    
    // Get outputs
    bam_stats = BAM_FLAGSTAT.out.stats
        .view { "BAM_STATS: ${it}" }
    length_files = SEQKIT_LENGTH.out.lengths
        .view { "LENGTH_FILES: ${it}" }
        .map { sample, file -> file }  // Extract just the file path
    size_stats = GET_FILE_SIZE.out.sizes
        .view { "SIZE_STATS: ${it}" }
    
    // Create mapping and size text files
    CREATE_MAPPING_FILE(bam_stats)
    CREATE_SIZE_FILE(size_stats)
    
    // Collect all files
    all_length_files = length_files.collect()
    all_mapping_files = CREATE_MAPPING_FILE.out.collect()
    all_size_files = CREATE_SIZE_FILE.out.collect()
    
    // Combine all QC stats
    COMBINE_QC_STATS(
        all_length_files,
        all_mapping_files,
        all_size_files,
        file("${params.qc_dir}/read_number_fail.txt")
    )
    
    emit:
    qc_report = COMBINE_QC_STATS.out
}

workflow reference_generate_module {
    take:
    filtered_gtf
    
    main:
    if (!params.reference_fasta || !params.ref_prot) {
        error "ERROR: --reference_fasta and --ref_prot are required for reference_generate module"
    }
    
    // Create output directories
    file("${params.outdir}/protein_ref").mkdirs()
    
    GENERATE_REFERENCE(filtered_gtf, file(params.reference_fasta), projectDir)
    REMOVE_HOMO(GENERATE_REFERENCE.out.orf_fa, file(params.ref_prot), projectDir)
    
    // Save final protein reference
    REMOVE_HOMO.out.filtered_orf
        .collectFile(name: "final_protein_reference.fa", storeDir: "${params.outdir}/protein_ref")
    
    emit:
    filtered_orf = REMOVE_HOMO.out.filtered_orf
}

workflow modification_module {
    take:
    mod_input_ch
    
    main:
    // Create output directories
    file("${params.outdir}/modifications").mkdirs()
    
    sample_data = mod_input_ch.map { sample, fast5_dir, clean_fq, isoforms_fa -> 
        tuple(sample, fast5_dir, clean_fq) 
    }
    isoforms = mod_input_ch.map { sample, fast5_dir, clean_fq, isoforms_fa -> 
        isoforms_fa 
    }.first()
    
    MODIFICATION_DETECTION(sample_data, isoforms, projectDir)
    
    // Outputs are already published by the MODIFICATION_DETECTION process
    emit:
    m6a = MODIFICATION_DETECTION.out.m6a_results
    m5c = MODIFICATION_DETECTION.out.m5c_results
    psu = MODIFICATION_DETECTION.out.psu_results
}

def proteomicsInput(value, option, executable = false) {
    if (!value) {
        error "ERROR: ${option} is required for proteomics module"
    }
    def path = file(value, glob: false).toAbsolutePath().normalize()
    if (!java.nio.file.Files.isRegularFile(path) || !java.nio.file.Files.isReadable(path)) {
        error "ERROR: ${option} must point to a readable regular file (${path})"
    }
    if (executable && !java.nio.file.Files.isExecutable(path)) {
        error "ERROR: ${option} must point to an executable file (${path})"
    }
    path
}

def proteomicsManifest(manifest) {
    def spectra = []
    def fields = []
    manifest.eachLine('UTF-8') { line, number ->
        if (line.trim()) {
            def columns = line.split('\t', -1)
            if (columns.length != 4 || line.contains('\u0000')) {
                error "ERROR: proteomics manifest line ${number} must contain four tab-separated columns"
            }
            if (!columns[0].trim()) {
                error "ERROR: proteomics manifest line ${number} mzML path is required"
            }
            if (columns[3] != 'DDA') {
                error "ERROR: proteomics manifest line ${number} must have acquisition type DDA"
            }
            def spectrum = manifest.parent.resolve(columns[0]).normalize()
            if (!spectrum.fileName.toString().toLowerCase(java.util.Locale.ROOT).endsWith('.mzml') ||
                !java.nio.file.Files.isRegularFile(spectrum) || !java.nio.file.Files.isReadable(spectrum) ||
                spectrum.toString().find(/[\t\r\n]/)) {
                error "ERROR: proteomics manifest line ${number} must reference a readable mzML file"
            }
            if (spectrum in spectra) {
                error "ERROR: proteomics manifest line ${number} repeats an mzML path"
            }
            spectra << spectrum
            fields << columns[1..3].toList()
        }
    }
    if (!spectra) {
        error "ERROR: proteomics manifest must contain at least one DDA mzML row"
    }
    [spectra, fields]
}

def fragpipeRamOverride() {
    def configured = params.fragpipe_ram_gb
    if (configured == null) {
        return null
    }

    // Nextflow retains only the last value of a repeated pipeline option.
    // Use the real argument boundaries, including for values with embedded quotes.
    def option = '--fragpipe_ram_gb'
    def argv = ProcessHandle.current().info().arguments().orElse(null)
    def cmdline = java.nio.file.Path.of('/proc/self/cmdline')
    if (!argv && java.nio.file.Files.isReadable(cmdline)) {
        argv = new String(java.nio.file.Files.readAllBytes(cmdline), 'UTF-8').split('\u0000')
    }
    if (!argv && workflow.commandLine?.contains(option)) {
        error 'ERROR: Cannot verify --fragpipe_ram_gb CLI arguments on this host; set FRAGPIPE task memory in a Nextflow config instead'
    }
    if (argv?.count { it == option || it.startsWith(option + '=') } > 1) {
        error 'ERROR: --fragpipe_ram_gb may be specified only once'
    }

    def value = configured.toString()
    if (!(value ==~ /[1-9][0-9]*/) || !value.isInteger()) {
        error 'ERROR: --fragpipe_ram_gb must be a positive integer in GiB'
    }
    return value.toInteger()
}


def fragpipeRuntime(launcher, philosopher) {
    def cli = params.fragpipe_cli.toString()
    if (!(cli in ['legacy', '23.1'])) {
        error 'ERROR: --fragpipe_cli must be legacy or 23.1'
    }
    def dryRun = params.fragpipe_dry_run
    if (!(dryRun instanceof Boolean)) {
        error 'ERROR: --fragpipe_dry_run must be true or false'
    }
    def expectedJar = cli == '23.1' ? 'fragpipe-23.1.jar' : 'fragpipe-20.0.jar'
    def lib = launcher.toRealPath().parent.resolve('../lib').normalize().toFile()
    def packaged = lib.listFiles()?.findAll { it.isFile() && it.name ==~ /fragpipe-\d[^\/]*\.jar/ }
    if (packaged && packaged*.name != [expectedJar]) {
        error "ERROR: --fragpipe_cli ${cli} does not match the FragPipe launcher package (${packaged*.name.sort().join(', ')})"
    }

    def executable = [launcher.toString()]
    if (cli == 'legacy') {
        if (params.fragpipe_ram_gb != null) {
            error 'ERROR: --fragpipe_ram_gb requires --fragpipe_cli 23.1'
        }
        if (params.fragpipe_tools_dir || params.fragpipe_python || dryRun) {
            error 'ERROR: --fragpipe_tools_dir, --fragpipe_python and --fragpipe_dry_run require --fragpipe_cli 23.1'
        }
        def msfragger = proteomicsInput(params.msfragger_jar, '--msfragger_jar').toString()
        def ionquant = proteomicsInput(params.ionquant_jar, '--ionquant_jar').toString()
        return [cli: cli, launcher: launcher.toString(), executable: executable + philosopher, readable: [msfragger, ionquant],
                args: ['--config-msfragger', msfragger, '--config-ionquant', ionquant, '--config-philosopher', philosopher],
                dryRun: false]
    }
    def ramGb = fragpipeRamOverride()
    if (params.msfragger_jar || params.ionquant_jar) {
        error 'ERROR: FragPipe 23.1 uses --fragpipe_tools_dir. Remove --msfragger_jar and --ionquant_jar'
    }
    if (!params.fragpipe_tools_dir) {
        error 'ERROR: --fragpipe_tools_dir is required for FragPipe 23.1'
    }
    def tools = file(params.fragpipe_tools_dir, glob: false).toAbsolutePath().normalize()
    if (!java.nio.file.Files.isDirectory(tools) || !java.nio.file.Files.isReadable(tools) ||
        !java.nio.file.Files.isExecutable(tools)) {
        error "ERROR: --fragpipe_tools_dir must point to a readable directory (${tools})"
    }

    // The supported 23.1 bundle matches the pinned engines and native library used in the dry-run.
    def members = ['MSFragger-4.3.jar', 'IonQuant-1.11.11.jar', 'diaTracer-1.3.3.jar',
                   'ext/bruker/libtimsdata-2-21-0-4.so'].collect { name ->
        proteomicsInput(tools.resolve(name), "--fragpipe_tools_dir member ${name}").toString()
    }
    def python = proteomicsInput(params.fragpipe_python, '--fragpipe_python', true).toString()
    def args = ['--config-tools-folder', tools.toString(), '--config-python', python]
    [cli: cli, launcher: launcher.toString(), executable: executable + python, readable: members,
     args: args, dryRun: dryRun, ramGb: ramGb]
}


workflow proteomics_module {
    take:
    protein_ref

    main:
    def manifest = proteomicsInput(params.proteomics_manifest, '--proteomics_manifest')
    def workflowConfig = proteomicsInput(params.proteom_workflow, '--proteom_workflow')
    def philosopher = proteomicsInput(params.philosopher_path, '--philosopher_path', true).toString()
    def fragpipe = proteomicsInput(params.fragpipe_path, '--fragpipe_path', true)
    def runtime = fragpipeRuntime(fragpipe, philosopher)
    def (spectra, manifestFields) = proteomicsManifest(manifest)

    ADD_DECOY(protein_ref, philosopher)
    FRAGPIPE(ADD_DECOY.out.decoy_db, workflowConfig, spectra, manifestFields, runtime)

    emit:
    fragpipe_output = FRAGPIPE.out.fragpipe_output
}

workflow suppa2_module {
    take:
    filtered_gtf
    quant_tsv
    
    main:
    // Create output directories
    file("${params.outdir}/suppa").mkdirs()
    
    SUPPA2_ANALYSIS(filtered_gtf, quant_tsv)
    
    // Outputs are already published by the SUPPA2_ANALYSIS process
    emit:
    ioe_files = SUPPA2_ANALYSIS.out.ioe_files
    psi_files = SUPPA2_ANALYSIS.out.psi_files
}

// Full pipeline workflow
workflow full_pipeline {
    // Read manifest and prepare channels
    samples_raw_fq_ch = Channel
        .fromPath(params.manifest)
        .splitCsv(sep: "\t", header: false)
        .map { row -> 
            def sample = row[0]
            def condition = row[1]
            def batch = row[2]
            def raw_fq = file(row[3].trim())
            def fast5_dir = file(row[4].trim())
            tuple(sample, condition, batch, raw_fq, fast5_dir)
        }
    
    // Step 1: Read chopping
    READS_CHOP(samples_raw_fq_ch.map { sample, condition, batch, raw_fq, fast5_dir -> 
        tuple(sample, raw_fq) 
    })
    clean_fq_ch = READS_CHOP.out.clean_fq
    
    // Create clean manifest
    clean_manifest_ch = samples_raw_fq_ch
        .combine(clean_fq_ch, by: 0)
        .map { sample, condition, batch, raw_fq, fast5_dir, clean_fq -> 
            "${sample}\t${condition}\t${batch}\t${clean_fq}"
        }
        .collectFile(
            name: "clean_manifest.txt",
            storeDir: "${params.outdir}",
            newLine: true
        )
    
    // Step 2: FLAIR pipeline
    flair_pip_module(clean_manifest_ch)
    
    // Step 3: Reference generation
    reference_generate_module(flair_pip_module.out.filtered_gtf)
    
    // Step 4: Quality report
    bam_ch = flair_pip_module.out.bam_files
        .flatMap { align_dir -> 
            def pattern = "${align_dir}/**/*.bam"
            file(pattern).collect { bam_file -> 
              def sample = bam_file.getSimpleName()
              tuple(sample, bam_file)
            }
        }
    quality_report_module(bam_ch, clean_fq_ch)
    
    // Step 5: SUPPA2 analysis
    suppa2_module(flair_pip_module.out.filtered_gtf, flair_pip_module.out.quant_tsv)
    
    // Step 6: Modification detection
    mod_input_ch = clean_fq_ch
        .join(samples_raw_fq_ch.map { s, c, b, fq, f5 -> tuple(s, f5) })
        .map { sample, clean_fq, fast5_dir -> tuple(sample, fast5_dir, clean_fq) }
        .combine(flair_pip_module.out.isoforms_fa.first())
    
    modification_module(mod_input_ch)
    
    // Step 7: Proteomics (optional)
    if (params.run_proteomics) {
        proteomics_module(reference_generate_module.out.filtered_orf.collect().first())
    }
}

// Main entry point
workflow {
    if (params.module == null) {
        // Run full pipeline
        log.info """
        ========================================
        Running FULL PIPELINE
        ========================================
        Output directory: ${params.outdir}
        Manifest: ${params.manifest}
        ========================================
        """
        full_pipeline()
        
    } else if (params.module == "reads_chop") {
        log.info """
        ========================================
        Running MODULE: reads_chop
        ========================================
        Required: --manifest
        Output: Clean FASTQ files
        ========================================
        """
        reads_chop_module()
        
        // Log completion message
        reads_chop_module.out.clean_fq.view { sample, fq ->
            "Processed: ${sample} -> ${fq}"
        }
        
    } else if (params.module == "flair_pip") {
        log.info """
        ========================================
        Running MODULE: flair_pip
        ========================================
        Required: --manifest (clean), --annotation_gtf, --reference_fasta
        Output: Isoforms, GTF, quantification
        ========================================
        """
        
        if (!params.manifest) {
            error "ERROR: --manifest is required for flair_pip module"
        }
        
        clean_manifest_ch = Channel.fromPath(params.manifest)
        flair_pip_module(clean_manifest_ch)
        
        // Log outputs
        flair_pip_module.out.filtered_gtf.view { "GTF: ${it}" }
        flair_pip_module.out.isoforms_fa.view { "Isoforms: ${it}" }
        flair_pip_module.out.quant_tsv.view { "Quantification: ${it}" }
        
    } else if (params.module == "quality_report") {
        log.info """
        ========================================
        Running MODULE: quality_report
        ========================================
        Required: --bam_dir, --clean_fastq_dir
        Output: QC statistics
        ========================================
        """
        
        if (!params.bam_dir || !params.clean_fastq_dir) {
            error "ERROR: --bam_dir and --clean_fastq_dir are required for quality_report module"
        }
        
        // Nextflow 23.10 does not match flat BAMs with the recursive glob.
        def bam_files = (files("${params.bam_dir}/*.bam", type: 'file') +
                         files("${params.bam_dir}/**/*.bam", type: 'file')).unique()
        def clean_fastqs = files("${params.clean_fastq_dir}/*.fastq.gz", type: 'file')
        if (!bam_files) {
            error "ERROR: No BAM files found in --bam_dir (${params.bam_dir})"
        }
        if (!clean_fastqs) {
            error "ERROR: No FASTQ files found in --clean_fastq_dir (${params.clean_fastq_dir})"
        }

        def bam_samples = bam_files.collect { bam -> tuple(bam.getSimpleName(), bam) }
        def fastq_samples = clean_fastqs.collect { fq ->
            tuple(fq.getName().replaceAll(/\.clean\.fastq\.gz$/, ''), fq)
        }
        def bam_ids = bam_samples.collect { it[0] }.toSet()
        def fastq_ids = fastq_samples.collect { it[0] }.toSet()
        if (bam_ids != fastq_ids) {
            error "ERROR: BAM and FASTQ sample IDs do not match. BAM only: ${(bam_ids - fastq_ids).sort()}; FASTQ only: ${(fastq_ids - bam_ids).sort()}"
        }

        bam_ch = Channel.fromList(bam_samples)
        clean_fq_ch = Channel.fromList(fastq_samples)
        
        quality_report_module(bam_ch, clean_fq_ch)
        
        // Log output
        quality_report_module.out.qc_report.view { "QC Report: ${it}" }
        
    } else if (params.module == "reference_generate") {
        log.info """
        ========================================
        Running MODULE: reference_generate
        ========================================
        Required: --filtered_gtf, --reference_fasta, --ref_prot
        Output: Protein reference FASTA
        ========================================
        """
        
        if (!params.filtered_gtf) {
            error "ERROR: --filtered_gtf is required for reference_generate module"
        }
        
        reference_generate_module(Channel.fromPath(params.filtered_gtf))
        
        // Log output
        reference_generate_module.out.filtered_orf.view { "Protein Reference: ${it}" }
        
    } else if (params.module == "modification") {
        log.info """
        ========================================
        Running MODULE: modification
        ========================================
        Required: --manifest, --isoforms_fa
        Output: m6A, m5C, psU modifications
        ========================================
        """
        
        if (!params.manifest || !params.isoforms_fa) {
            error "ERROR: --manifest and --isoforms_fa are required for modification module"
        }
        
        mod_input_ch = Channel
            .fromPath(params.manifest)
            .splitCsv(sep: "\t", header: false)
            .map { row -> 
                def sample = row[0]
                def fast5_dir = file(row[3].trim())
                def clean_fq = file(row[2].trim())
                tuple(sample, fast5_dir, clean_fq)
            }
            .combine(Channel.fromPath(params.isoforms_fa))
        
        modification_module(mod_input_ch)
        
        // Log outputs
        modification_module.out.m6a.view { "m6A results: ${it}" }
        modification_module.out.m5c.view { "m5C results: ${it}" }
        modification_module.out.psu.view { "psU results: ${it}" }
        
    } else if (params.module == "proteomics") {
        log.info """
        ========================================
        Running MODULE: proteomics
        ========================================
        Required: --protein_ref, --proteomics_manifest, --proteom_workflow
        Output: FragPipe results
        ========================================
        """
        
        def proteinRef = proteomicsInput(params.protein_ref, '--protein_ref')
        proteomics_module(proteinRef)

        // Log output
        proteomics_module.out.fragpipe_output.view { "FragPipe results: ${it}" }
        
    } else if (params.module == "suppa2") {
        log.info """
        ========================================
        Running MODULE: suppa2
        ========================================
        Required: --filtered_gtf, --quant_tsv
        Output: Alternative splicing events
        ========================================
        """
        
        if (!params.filtered_gtf || !params.quant_tsv) {
            error "ERROR: --filtered_gtf and --quant_tsv are required for suppa2 module"
        }
        
        suppa2_module(
            Channel.fromPath(params.filtered_gtf),
            Channel.fromPath(params.quant_tsv)
        )
        
        // Log outputs
        suppa2_module.out.ioe_files.view { "IOE files: ${it}" }
        suppa2_module.out.psi_files.view { "PSI files: ${it}" }
        
    } else {
        error """
        ERROR: Unknown module '${params.module}'
        
        Available modules:
        - reads_chop
        - flair_pip
        - quality_report
        - reference_generate
        - modification
        - proteomics
        - suppa2
        
        Or run without --module to execute the full pipeline
        """
    }
}

workflow.onComplete {
    log.info """
    ========================================
    Pipeline completed!
    Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration: ${workflow.duration}
    ========================================
    """
}