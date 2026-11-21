# ONT F1-1

The 12-read FASTQ and `fast5/F1-1-paired-mini.fast5` have identical read ID
sets. A second FASTQ has 1,500 different reads from the same sample and is
for FASTQ parsing and read-length summaries. Use the 12-read FASTQ in the
paired manifest.

The files are included in this review branch. Run these commands
from the repository root. Python 3 with h5py is required for the checker.

```bash
example_dir="$PWD/examples/real-mini/ont-F1-1"
(cd "$example_dir" && sha256sum -c SHA256SUMS.txt)
python3 "$example_dir/check_realmini.py" "$example_dir" --include-qc
bash "$example_dir/make_modification_manifest.sh" "$example_dir" \
  > "$example_dir/modification.manifest.tsv"
```

The headerless standalone `--module modification` manifest has four
tab-separated columns in this order. Sample label, unused placeholder, absolute paired
FASTQ path, and absolute FAST5 directory path. The module also requires
`--isoforms_fa` and a complete local installation. The format check covers
FASTQ, HDF5, exact IDs, and readable signal. One sample provides no
biological replicates.
