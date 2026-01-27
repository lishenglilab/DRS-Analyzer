# DRS-Analyzer

A comprehensive Nextflow pipeline for analyzing Direct RNA Sequencing (DRS) data from Oxford Nanopore Technologies.

## Overview

DRS-Analyzer integrates multiple analysis modules for complete DRS data processing, including:
- ✅ Adapter trimming
- ✅ Transcript identification and quantification
- ✅ Quality control reporting
- ✅ Post-transcriptional modification detection (m6A, m5C, pseudouridine)
- ✅ Alternative splicing analysis
- ✅ ORF prediction and protein reference generation
- ✅ Proteomics integration

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Requirements](#input-requirements)
- [Module Usage](#module-usage)
- [Output Description](#output-description)
- [Test Data](#test-data)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

---

## Installation

### Prerequisites

- Nextflow (version 23.10.1 or higher)
- Docker (recommended) or Conda

### Option 1: Using Docker (Recommended)

```bash
# Pull the pre-configured Docker image
docker pull yvzeng/drs-analyzer-env:latest

# Clone the repository
git clone https://github.com/yourusername/DRS-Analyzer.git
cd DRS-Analyzer
```

### Option 2: Manual Installation

#### Software Dependencies

**Nextflow**
```bash
curl -s https://get.nextflow.io | bash
```

**Python Packages**
```bash
conda create -n drs-analyzer python=3.7
conda activate drs-analyzer
pip install pyfastx==0.8.4 pysam==0.22.0 pandas==1.1.5 scikit-learn==0.22 biopython
```

**R Packages**
```R
install.packages(c("getopt", "tidyverse", "plyr", "ggsci", "Rtsne", "Hmisc", "stringr"))
BiocManager::install(c("diann", "ComplexHeatmap", "Biostrings", "rtracklayer"))
```

**Bioinformatics Tools**
```bash
conda install -c bioconda \
  bedops=2.4.35 \
  tombo=1.5.1 \
  flair=2.0.0 \
  blast=2.17.0 \
  gffread=0.11.8 \
  porechop=0.2.4 \
  samtools=1.20

# Install NanoPsu
git clone https://github.com/sihaohuanguc/Nanopore_psU.git
cd Nanopore_psU
pip install .

# Install SUPPA2
pip install SUPPA
```

**FragPipe (Optional, for proteomics module)**
- Download from: https://fragpipe.nesvilab.org/
- Version: v22.0

---

## Quick Start

### Full Pipeline

```bash
nextflow run main.nf \
  --manifest /path/to/manifest.txt \
  --outdir /path/to/output/ \
  --annotation_gtf /path/to/annotation.gtf \
  --reference_fasta /path/to/genome.fa \
  --ref_prot /path/to/reference_protein.fasta
```

### Run Individual Modules

```bash
# Example: Transcript identification only
nextflow run main.nf \
  --module flair_pip \
  --manifest /path/to/manifest.txt \
  --annotation_gtf /path/to/annotation.gtf \
  --reference_fasta /path/to/genome.fa \
  --outdir /path/to/output/
```

---

## Input Requirements

### 1. Manifest File

Tab-separated file with the following format:

```
Sample1      control      batch1   /path/to/Sample1.fq.gz   /path/to/Sample1/fast5/
Sample2      treated      batch1   /path/to/Sample2.fq.gz   /path/to/Sample2/fast5/
Sample3      treated      batch2   /path/to/Sample3.fq.gz   /path/to/Sample3/fast5/
```

**Column Descriptions:**
- `Sample_ID`: Unique sample identifier
- `Condition`: Experimental condition
- `Batch`: Batch identifier
- `FASTQ_Path`: Path to raw FASTQ file
- `FAST5_Directory`: Path to FAST5 directory (required for modification detection)

### 2. Reference Files

- **Genome FASTA**: Reference genome sequence
- **Annotation GTF**: Gene annotation in GTF format
- **Protein Reference**: Known protein sequences (for proteomics and ORF filtering)

---

## Module Usage

### 1. Read Chopping

Removes adapters from raw FASTQ files.

```bash
nextflow run main.nf \
  --module reads_chop \
  --manifest /path/to/manifest.txt \
  --outdir /path/to/output/
```

**Output:**
```
clean_fq/
├── Sample1.clean.fastq.gz
├── Sample2.clean.fastq.gz
└── Sample3.clean.fastq.gz
```

### 2. FLAIR Pipeline

Identifies and quantifies transcripts.

```bash
nextflow run main.nf \
  --module flair_pip \
  --manifest /path/to/clean_manifest.txt \
  --annotation_gtf /path/to/annotation.gtf \
  --reference_fasta /path/to/genome.fa \
  --outdir /path/to/output/
```

**Output:**
```
flair/
├── filtered_isoforms.gtf           # Filtered transcripts (TPM ≥ 0.1)
├── isoforms.fa                     # Transcript sequences
├── quant.tpm.tsv                   # All transcript quantification
└── flair_quantify_filter.tpm.tsv  # Filtered quantification
```

### 3. Quality Report

Generates QC statistics.

```bash
nextflow run main.nf \
  --module quality_report \
  --bam_dir /path/to/bam_files/ \
  --clean_fastq_dir /path/to/clean_fastq/ \
  --outdir /path/to/output/
```

**Output:**
```
table/
└── QC.xls
```

**QC.xls columns:**
| Column | Description |
|--------|-------------|
| sample_id | Sample identifier |
| Size_GB | FASTQ file size (GB) |
| Total_reads | Total read count |
| pass_percent | QC pass rate (%) |
| Median_pass_read_length | Median read length (bp) |
| Mean_pass_read_length | Mean read length (bp) |
| Mapped_read | Number of mapped reads |
| Mapped_ratio | Mapping rate (%) |

### 4. Alternative Splicing

Identifies alternative splicing events using SUPPA2.

```bash
nextflow run main.nf \
  --module suppa2 \
  --filtered_gtf /path/to/filtered_isoforms.gtf \
  --quant_tsv /path/to/quant.tpm.tsv \
  --outdir /path/to/output/
```

**Output:**
```
suppa/
├── transcriptEvents_SE_strict.ioe      # Skipped exon
├── transcriptEvents_A3_strict.ioe      # Alternative 3' SS
├── transcriptEvents_A5_strict.ioe      # Alternative 5' SS
├── transcriptEvents_MX_strict.ioe      # Mutually exclusive exons
├── transcriptEvents_RI_strict.ioe      # Retained intron
├── transcriptEvents.all.events.ioe     # All events
└── transcriptEvents.all.events.psi.psi # PSI values
```

### 5. Post-transcriptional Modifications

Detects m6A, m5C, and pseudouridine modifications.

```bash
nextflow run main.nf \
  --module modification \
  --manifest /path/to/clean_manifest.txt \
  --isoforms_fa /path/to/isoforms.fa \
  --outdir /path/to/output/
```

**Output:**
```
modifications/
├── Sample1/
│   ├── m6a/Sample1/m6A_output.bed
│   ├── m5c/Sample1/m5c.prediction.txt
│   └── psu/Sample1/alignment/prediction.csv
└── ...
```

#### m6A Results (BED format)

| Column | Description |
|--------|-------------|
| 1 | Transcript ID |
| 2 | Start position |
| 3 | End position |
| 4 | DRACH motif |
| 5 | Modification ID |
| 6 | Strand |
| 7 | Modification fraction (0-1) |
| 8 | Read coverage |

**Example:**
```
mETP-T000002221  272  273  AGACT  mETP-T000002221:272:AGACT:+  +  0.2  5
```

#### m5C Results (Tab-delimited)

| Column | Description |
|--------|-------------|
| transcript | Transcript ID |
| position | Cytosine position |
| base | Base (C) |
| fraction | Modification fraction (0-1) |

**Example:**
```
transcript           position  base  fraction
mETP-T000002206_F    70        C     1.0000
```

#### Pseudouridine Results (CSV)

| Column | Description |
|--------|-------------|
| transcript_ID | Transcript ID |
| position | Uridine position |
| base_type | Base type (T) |
| coverage | Read coverage |
| prob_unmodified | Unmodified probability |
| prob_modified | Modified probability |

**Example:**
```csv
transcript_ID,position,base_type,coverage,prob_unmodified,prob_modified
mETP-T000001375_F,122,T,22,0.955,0.045
```

### 6. Protein Reference Generation

Predicts ORFs and generates protein database.

```bash
nextflow run main.nf \
  --module reference_generate \
  --filtered_gtf /path/to/filtered_isoforms.gtf \
  --reference_fasta /path/to/genome.fa \
  --ref_prot /path/to/uniprot_reference.fasta \
  --outdir /path/to/output/
```

**Output:**
```
protein_ref/
├── ORF_no_redundancy.fa              # Non-redundant ORFs
└── final_protein_reference.fa        # Final protein database
```

### 7. Proteomics Analysis

Analyzes mass spectrometry data using FragPipe.

```bash
nextflow run main.nf \
  --module proteomics \
  --protein_ref /path/to/protein_reference.fasta \
  --proteomics_manifest /path/to/manifest.fp-manifest \
  --proteom_workflow /path/to/workflow.workflow \
  --philosopher_path /path/to/philosopher \
  --fragpipe_path /path/to/fragpipe \
  --msfragger_jar /path/to/MSFragger.jar \
  --ionquant_jar /path/to/IonQuant.jar \
  --outdir /path/to/output/
```

**Output:**
```
protein/fragpipe_output/
├── combined_protein.tsv          # Protein quantification
├── combined_peptide.tsv          # Peptide quantification
├── combined_ion.tsv              # Ion quantification
├── MSstats.csv                   # MSstats format
└── [sample_directories]/         # Per-sample results
```

---

## Output Description

### Directory Structure

```
output/
├── clean_fq/                    # Cleaned FASTQ files
├── flair/                       # Transcript identification results
│   ├── align/                   # BAM alignments
│   ├── filtered_isoforms.gtf   # Main transcript annotation
│   ├── isoforms.fa             # Transcript sequences
│   └── quant.tpm.tsv           # Quantification
├── table/
│   └── QC.xls                  # Quality control summary
├── suppa/                      # Alternative splicing events
├── modifications/              # RNA modifications
│   ├── Sample1/
│   │   ├── m6a/               # m6A modifications
│   │   ├── m5c/               # m5C modifications
│   │   └── psu/               # Pseudouridine
├── protein_ref/                # Protein reference database
└── protein/                    # Proteomics results (optional)
```

---

## Test Data

We provide a small test dataset to help you verify the installation and test individual modules. 

### Download Test Data

The test dataset is available at: **http://www.shenglilabs.com/mETP/**


### Test Data Structure

```
DRS-Analyzer-testdata/
├── fast5/                      # FAST5 files for modification detection
│   ├── F1-1/
│   │   └── PAJ07356_pass_546ffc41_0.fast5
│   ├── F1-2/
│   │   └── PAJ07356_pass_546ffc41_0.fast5
│   └── F1-3/
│       └── PAJ07356_pass_546ffc41_0.fast5
├── fq/                         # Clean FASTQ files
│   ├── F1-1.clean.fastq.gz
│   ├── F1-2.clean.fastq.gz
│   └── F1-3.clean.fastq.gz
└── protein/                    # Proteomics test data
    ├── LFQ-MBR.workflow       # FragPipe workflow configuration
    ├── M2brain.mzML           # Mass spectrometry data
    └── manifest.fp-manifest   # FragPipe manifest file
```

### Important Notes

⚠️ **The test dataset is intentionally small and is designed for testing individual modules only.**

- **Dataset size**: Minimal FAST5/FASTQ files (1 file per sample)
- **Purpose**: Installation verification and module testing
- **Limitation**: NOT suitable for running the full pipeline
- **Recommendation**: Use your own complete dataset for actual analysis

### Testing Individual Modules

#### Example 1: Test Quality Report Module

```bash
# Assuming you have BAM files from FLAIR pipeline
nextflow run main.nf \
  --module quality_report \
  --bam_dir /path/to/testdata/bam_files/ \
  --clean_fastq_dir /path/to/testdata/fq/ \
  --outdir ./test_output/
```

#### Example 2: Test Modification Detection Module

First, create a test manifest file (`test_manifest.txt`):
```
F1-1    control    batch1    /path/to/testdata/fq/F1-1.clean.fastq.gz    /path/to/testdata/fast5/F1-1/
F1-2    control    batch1    /path/to/testdata/fq/F1-2.clean.fastq.gz    /path/to/testdata/fast5/F1-2/
F1-3    treated    batch1    /path/to/testdata/fq/F1-3.clean.fastq.gz    /path/to/testdata/fast5/F1-3/
```

Then run:
```bash
nextflow run main.nf \
  --module modification \
  --manifest test_manifest.txt \
  --isoforms_fa /path/to/isoforms.fa \
  --outdir ./test_output/
```

#### Example 3: Test Proteomics Module

```bash
nextflow run main.nf \
  --module proteomics \
  --protein_ref /path/to/protein_reference.fasta \
  --proteomics_manifest /path/to/testdata/protein/manifest.fp-manifest \
  --proteom_workflow /path/to/testdata/protein/LFQ-MBR.workflow \
  --philosopher_path /path/to/philosopher \
  --fragpipe_path /path/to/fragpipe \
  --msfragger_jar /path/to/MSFragger.jar \
  --ionquant_jar /path/to/IonQuant.jar \
  --outdir ./test_output/
```

### Expected Results

After running the test, you should see:
- ✅ Nextflow pipeline completes without errors
- ✅ Output files are generated in the specified output directory
- ✅ Log files show successful completion

If tests pass, your installation is correct and you can proceed with analyzing your own data.

### Using Your Own Data

For production analysis with complete datasets:

1. Prepare your full manifest file with all samples
2. Ensure you have adequate computational resources (see [Performance Tips](#performance-tips))
3. Run the full pipeline or individual modules as needed


## Version History

- v1.0.0 (2026-01-27): Initial release
  - Complete pipeline implementation
  - All seven modules functional
  - Docker image available
  - Test dataset provided

---

**Developed by [lishenglilab]**

⭐ If you find this tool useful, please give it a star on GitHub!
