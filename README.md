Nanopore Direct RNA Sequencing (DRS) Analysis Pipeline
=
Overview
-
This Nextflow pipeline provides a comprehensive analysis solution for Nanopore Direct RNA Sequencing (DRS) data. The workflow performs end-to-end processing including quality control, transcript assembly, quantification, modification detection, alternative splicing analysis, and optional proteomics integration, and this pipeline is executed in the Linux system..

Requirements
-
*Nextflow<br>
*Python：<br>
`subprocess`
`sys`
`os`
`argparse`
`datetime`
`pyfastx==0.8.4`
`pysam==0.22.0`
`re`
`pandas==1.1.5`<br>
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
`tombo==1.5.1`
`NanoPsu`
`flair==v2.0.0`
`blast==2.17.0+`
`gffread==v0.11.8`
`suppa2==v2.4`
`porechop==0.2.4`
`samtools==1.20`
`fragpipe==v22.0`(optional)<br>

Users can also complete the environment configuration by pulling the docker image<br>
`docker pull yvzeng/drs-analyzer-env:latest`<br>

Installation
-
`git clone https://github.com/lishenglilab/DRS-Analyzer.git`<br>
`cd DRS-Analyzer`

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
