# M2brain DDA

This separate mouse proteomics acquisition contains 52 MS1 spectra and
1,144 MS2 spectra in an indexed mzML. Use the accompanying 57-record
target-only FASTA. DRS-Analyzer builds the decoys from that FASTA with
Philosopher. See [data provenance](DATA_PROVENANCE.md) for UniProt
attribution and the reference selection.

The data are included in this review branch. Run from the repository root.
The format checker uses the Python 3 standard library.

```bash
example_dir="$PWD/examples/real-mini/proteomics-M2brain"
(cd "$example_dir" && sha256sum -c SHA256SUMS.txt)
python3 "$example_dir/validate_mzml.py" "$example_dir/M2brain-DDA-mini.mzML" \
  --fasta "$example_dir/mini_target_only.fasta"
bash "$example_dir/make_fragpipe_manifest.sh" "$example_dir" \
  > "$example_dir/single.fp-manifest"
```

The manifest contains one DDA row with the mzML filename, experiment,
empty biological-replicate field and acquisition type. Keep the manifest
beside the mzML because paths are resolved from that directory.

Install and license FragPipe 23.1, MSFragger 4.3, IonQuant 1.11.11,
diaTracer 1.3.3, and the other required components through their official
channels. In the tools directory, place the three versioned JARs and the
matching MSFragger `ext` tree as described in the project README. Set the
tool paths below and point `JAVA_HOME` to Java 11 or later. Supply an
executable Python 3.9–3.11 for FragPipe and Philosopher 5.1.2 for building
the decoy database. [FragPipe 23.1](https://github.com/Nesvilab/FragPipe/releases/tag/23.1) · [MSFragger](https://msfragger.nesvilab.org/) · [IonQuant](https://ionquant.nesvilab.org/)

```groovy
// resources.config
process {
    withName: ADD_DECOY {
        cpus = 1
        memory = '1 GB'
    }
    withName: FRAGPIPE {
        cpus = 2
        memory = '5 GB'
        time = '45 min'
    }
}
```

```bash
export JAVA_HOME=/path/to/java
export PATH="$JAVA_HOME/bin:$PATH"
nextflow -c resources.config run main.nf \
  --module proteomics \
  --protein_ref "$example_dir/mini_target_only.fasta" \
  --proteomics_manifest "$example_dir/single.fp-manifest" \
  --proteom_workflow "$PWD/bin/protein/LFQ_DDA_single_run_FragPipe23.1.workflow" \
  --fragpipe_cli 23.1 \
  --fragpipe_path /path/to/fragpipe-23.1/bin/fragpipe \
  --philosopher_path /path/to/philosopher-v5.1.2 \
  --fragpipe_tools_dir /path/to/engines \
  --fragpipe_python /path/to/python3.11/bin/python \
  --fragpipe_ram_gb 5 \
  --fragpipe_dry_run true \
  --outdir "$PWD/m2brain-dry-run"
```

The dry-run checks staging, builds decoys, and writes the FragPipe plan.
The provided workflow disables match-between-runs (MBR). The 5 GiB engine
setting must fit the Nextflow process allocation. A server preflight
used two CPUs and an 8 GiB job limit, including process overhead.
An actual search needs a compatible, licensed installation and local
execution resources.

A historical private search of the original slice exercised protein
identification and single-run IonQuant MS1 intensity reporting. IonQuant
skips MaxLFQ with one experiment. The de-identified files have format and
configuration validation. This selected demonstration provides no public
PD benchmark, cross-run LFQ/MBR validation, full-database FDR or evidence of
tdORF novelty or equivalence to Proteome Discoverer.
