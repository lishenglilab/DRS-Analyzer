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
- [Run the synthetic tests](#run-the-synthetic-tests)
- [Test Data](#test-data)

---

## Installation

### Prerequisites

- Nextflow (version 23.10.1 or higher)
- Docker (recommended) or Conda

### Option 1: Using Docker (Recommended)

```bash
# Pull the pre-configured Docker image
docker pull yvzeng/drs-analyzer-env:latest

# You can also create the container by yourself based on the Dockerfile.
git clone https://github.com/lishenglilab/DRS-Analyzer.git
cd DRS-Analyzer
docker build -t drs-analyzer-env .

# Clone the repository
git clone https://github.com/lishenglilab/DRS-Analyzer.git
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
- FragPipe 23.1 is the supported modern launcher. The verified dry-run uses Java 11,
  MSFragger 4.3, IonQuant 1.11.11, diaTracer 1.3.3, Philosopher 5.1.2,
  MSBooster 1.3.17, Percolator 3.7.1 and a Python 3.11 executable.
- The legacy per-tool flags target the FragPipe 20.0 launcher. They have mock
  regression coverage. Use explicit CLI mode and check your installed toolset.

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

The module accepts label-free DDA (regular MS) with mzML. This example verifies
the workflow and tool configuration using a small target-only FASTA and one mzML.
The dry-run builds a decoy database for the supplied FASTA, then asks FragPipe
to write a plan. Use a small FASTA to limit this preflight.

```bash
nextflow -c resources.config run main.nf \
  --module proteomics \
  --protein_ref /path/to/small_target_only.fasta \
  --proteomics_manifest /path/to/single.fp-manifest \
  --proteom_workflow bin/protein/LFQ_DDA_single_run_FragPipe23.1.workflow \
  --fragpipe_cli 23.1 \
  --fragpipe_path /path/to/fragpipe-23.1/bin/fragpipe \
  --philosopher_path /path/to/philosopher-v5.1.2 \
  --fragpipe_tools_dir /path/to/engines \
  --fragpipe_python /path/to/python3.11/bin/python \
  --fragpipe_ram_gb 5 \
  --fragpipe_dry_run true \
  --outdir /path/to/scratch/output/
```

The manifest can have any filename. Each nonblank line contains four
tab-separated fields: mzML path, experiment, biological replicate (may be empty),
and `DDA`. Relative paths resolve from the manifest's directory. Nextflow stages
copies of the mzML files and rebuilds the manifest with task-local paths. Allow
scratch space for those copies. The workflow requires one `database.db-path`
setting and `database.decoy-tag=rev_`. The single-run template keeps IonQuant
LFQ enabled and sets `ionquant.mbr=0`. Use one DDA mzML for this template.

In 23.1 mode, `--fragpipe_tools_dir` must expose `MSFragger-4.3.jar`,
`IonQuant-1.11.11.jar`, `diaTracer-1.3.3.jar` and the matching MSFragger
`ext` tree at its top level. Readable symlinks to the installed files are
supported. `--fragpipe_python` takes an executable file. Bind the tools and
their installation directories at matching paths when running in a container.
The standalone Philosopher executable builds the decoy database. FragPipe
23.1 resolves its own bundled Philosopher for downstream processing. Install
FragPipe 23.1 with Java 11 or later and use a Python version supported by that
release (3.9–3.11). Obtain MSFragger, IonQuant and diaTracer through their
official distribution and follow the terms for your installation and use.
[FragPipe 23.1 release](https://github.com/Nesvilab/FragPipe/releases/tag/23.1) ·
[MSFragger](https://msfragger.nesvilab.org/) ·
[IonQuant](https://ionquant.nesvilab.org/) ·
[Component licenses](https://github.com/Nesvilab/FragPipe/blob/23.1/FRAGPIPE-THIRD-PARTY-LICENSES.txt)

FragPipe 23.1 receives exactly one `--ram` in both modes.
`--fragpipe_ram_gb` accepts a positive whole number of GiB up to the Nextflow
`FRAGPIPE` task memory allocation. With the option omitted, the pipeline
uses the task allocation in whole GiB. Supply the option once. If process
arguments are unavailable on your host, use the task memory setting instead.
Save a site-specific `resources.config` with an allocation appropriate for your
scheduler or container. The following settings and the example RAM value
were used with an enforced 8 GiB job limit and leave room for non-heap memory
and other processes.

```groovy
process {
    withName: FRAGPIPE {
        cpus = 2
        memory = '5 GB'
        time = '45 min'
    }
}
```

Set the scheduler or container memory limit separately.
[FragPipe headless options](https://fragpipe.nesvilab.org/docs/tutorial_headless.html)

Philosopher can launch an embedded ProteinProphet executable from its
temporary directory. If ProteinProphet fails with `Permission denied`,
inspect the mount options of the reported temporary path for `noexec`.
For an isolated container, use a small, owner-only executable temporary
mount for that tool, with `nosuid` and `nodev`, and keep other scratch
mounts non-executable. Run as a non-root user with read-only input and
software mounts and a restricted network. Apply the executable mount
only inside the disposable job environment.

Use `--fragpipe_cli legacy` with `--msfragger_jar` and `--ionquant_jar` for the
20.0 per-tool CLI. The legacy mode is the default for prior invocations.
Use the modern tools directory, Python, RAM and dry-run options with
23.1 mode. The launcher checks a versioned FragPipe JAR in the distribution
when the launcher exposes one. Wrapper executables require an explicit and
accurate mode selection.

The dry-run publishes `fragpipe.log`, `fragpipe.job`, the exported
`fragpipe.workflow` and `fragpipe-files.fp-manifest` under
`protein/fragpipe_output/`. The pipeline checks the log for errors, the 23.1
version and the no-execution confirmation, plus nonempty plan files. Set
`--fragpipe_dry_run false` for an actual run. The real run checks for a
nonempty log and at least one nonempty TSV/CSV.

In a prior private server smoke test using a separately managed, resource-limited
deployment, FragPipe 23.1, MSFragger 4.3, Philosopher 5.1.2 and IonQuant 1.11.11
completed a single DDA acquisition. After final filtering, 18 of 20 peptide
anchors selected in advance from an existing Proteome Discoverer analysis
were recovered. With MBR disabled, 26 of 37 ion rows and 19 of 24 protein rows
had finite positive MS1 intensities. IonQuant skipped MaxLFQ calculation
for the single experiment. This PD-informed test validates the software
workflow on a selected example. Full-database and tdORF-specific FDR,
PD/FragPipe equivalence, cross-run LFQ and MBR are outside its scope.
The repository CI exercises synthetic fixtures and the FragPipe interface
with mock tools. Running the real workflow requires a compatible installation
and local input data.

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

## Run the synthetic tests

On Linux x86-64, create and activate the lightweight test environment with Micromamba, then run the same command used in GitHub Actions.

```bash
micromamba create --strict-channel-priority -f tests/environment.yml
micromamba activate drs-tests
bash tests/run.sh
```

The environment specifies Nextflow 23.10.1 with Java 11, SeqKit 2.12.0, Samtools 1.20, R with tidyverse and data.table, and Python with Biopython. The test driver prints an output directory under `TMPDIR` (or the system temporary directory) with `versions.log`, `results.log`, command logs and results for each test. Each run starts with fresh work directories. CI sets a 20-minute limit for dependency installation and testing, and saves logs on test failure.

The tiny examples in `tests/data/` are synthetic FASTQ, SAM, FASTA and WIG. The driver generates the gzip FASTQ and BAM, verifies SeqKit lengths and Samtools alignment counts, checks the R and Python reports against fixed values, and runs the `quality_report` command on both flat and nested BAM directories. Coverage includes invalid FASTQ, unknown transcript IDs, and missing or mismatched QC inputs. The suite also checks the proteomics manifest, staging, failure handling and publication, plus both FragPipe CLI forms and 23.1 dry-run output checks with mock tools. Adapter trimming, FLAIR, SUPPA2, full modification detection and real proteomics searches are outside its coverage.

For `--module quality_report`, use matching sample names such as `toy.bam` and `toy.clean.fastq.gz`. BAM files can be directly inside `--bam_dir` or in subdirectories. The command returns a nonzero exit when a BAM or FASTQ input set is empty or when their sample IDs differ. It writes a tab-separated report to `<outdir>/table/QC.xls`. An optional `<qc_dir>/read_number_fail.txt` contains sample IDs and failed FASTQ line counts. Omit it or leave it empty for the all-pass case.

---

## Test Data

The eight small fixtures in [tests/data/](tests/data/) are synthetic.
GitHub Actions runs 68 synthetic QC and mock helper/interface tests using
`bash tests/run.sh`.

Two real-data examples are available on this review branch. The
[Real-mini guide](examples/real-mini/) records their scope and file checks.

- **ONT F1-1** is one mouse direct-RNA sample. The paired FASTQ and FAST5
  contain the same 12 reads. A separate 1,500-read QC FASTQ from that sample
  is for FASTQ parsing and read-length summaries and shares zero read IDs
  with the paired set.
- **M2brain DDA** is a separate mass-spectrometry acquisition. The example
  combines one indexed mzML slice (52 MS1 and 1,144 MS2 spectra) with
  57 full-length reviewed mouse Swiss-Prot targets. The curated release
  candidate passed format and configuration validation. A historical
  private server search used the original mzML and original FASTA, with a
  separate functional audit.

The historical Google Drive archive used three copies of a FAST5 file.
Only F1-1 had read IDs in common with its FASTQ in the audited archive.
Use the Real-mini guide for the verified file pairing and defined checks.

These examples are separate acquisitions. The ONT checks cover file
structure and read pairing. The DDA example documents a single-run
configuration and MS1 quantification context. Full modification analysis,
FLAIR, SUPPA2, cross-sample MaxLFQ/MBR and real database searching are
outside the GitHub Actions test suite.

---

**Developed by [lishenglilab]**

⭐ If you find this tool useful, please give it a star on GitHub!
