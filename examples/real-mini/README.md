# Real-mini examples

These examples are available on this review branch under the linked data terms.
The repository's eight synthetic fixtures remain in `tests/data/`, and
GitHub Actions has run 68 synthetic and mock tests using those inputs.

| Example | Contents | Verified use |
|---|---|---|
| [ONT F1-1](ont-F1-1/) | 12 matched FASTQ/FAST5 reads and an independent 1,500-read QC FASTQ | Input structure, exact read pairing and read lengths |
| [M2brain DDA](proteomics-M2brain/) | One indexed mzML with 52 MS1 and 1,144 MS2 spectra, plus 57 target-only UniProtKB/Swiss-Prot sequences | Validated data format and a single-run FragPipe plan |

F1-1 and M2brain are independent acquisitions. Follow each directory's
manifest and checksum instructions. The [data terms](DATA_LICENSE)
and the dataset provenance pages identify the respective data sources.
