Nanopore Direct RNA Sequencing (DRS) Analysis Pipeline
=
Overview
-
This Nextflow pipeline provides a comprehensive analysis solution for Nanopore Direct RNA Sequencing (DRS) data. The workflow performs end-to-end processing including quality control, transcript assembly, quantification, modification detection, alternative splicing analysis, and optional proteomics integration.

Requirements
-
*Nextflow<br>
*Python：<br>
`subprocess`
`sys`
`os`
`argparse`
`datetime`
`pyfastx`
`pysam`
`re`
`pandas`<br>
*R:<br>
`getopt`
`tidyverse`
`diann`
`plyr`
`ComplexHeatmap`
`ggsci`
`Rtsne`
`Hmisc`
`Biostrings`
`stringr`
`rtracklayer`<br>
*other tools:<br>
`tombo`
`NanoPsu`
`flair`
`blast`
`gffread`
`suppa2`
`porechop`
`samtools`
`fragpipe`(optional)<br>

Users can also complete the environment configuration by pulling the docker image<br>
`docker pull yvzeng/drs-analyzer-env:latest`<br>

Installation
-
`git clone https://github.com/lishenglilab/mETP.git`<br>
`cd mETP`

Input Requirements
-
Manifest File
A tab-delimited manifest file is required with the following columns:<br>
'Sample ID', 'Condition', 'Batch', 'Path to raw FASTQ file', 'Path to FAST5 directory'

Example manifest:<br>
sample1  control  batch1  /path/to/sample1.fastq.gz  /path/to/fast5/sample1<br>
sample2  treated  batch1  /path/to/sample2.fastq.gz  /path/to/fast5/sample2<br>

Reference Files<br>
*Genome FASTA file<br>
*GTF annotation file<br>
*Reference protein database (FASTA format)<br>

Running the Pipeline
-
Basic Command
-
`nextflow run main.nf \`<br>
`    --manifest samples.manifest \`<br>
`    --reference_fasta GRCm39.genome.fa \`<br>
`    --annotation_gtf gencode.vM29.annotation.gtf \`<br>
`    --ref_prot uniprot_sprot_mouse.fasta \`<br>
`    --outdir results`<br>

Proteomics Mode
-
`nextflow run main.nf \`<br>
`    --manifest samples.manifest \`<br>
`    --reference_fasta GRCm39.genome.fa \`<br>
`    --annotation_gtf gencode.vM29.annotation.gtf \`<br>
`    --ref_prot uniprot_sprot_mouse.fasta \`<br>
`    --outdir results`<br>
`    --run_proteomics true \`<br>
`    --proteomics_manifest proteomics.manifest`<br>
