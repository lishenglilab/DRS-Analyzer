#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    printf 'Usage: %s EXAMPLE_DIR\n' "$0" >&2
    exit 2
fi
example_dir=$(cd -- "$1" && pwd -P)
fastq="$example_dir/F1-1-paired-mini.fastq.gz"
fast5_dir="$example_dir/fast5"
[[ -r "$fastq" && -f "$fastq" && -d "$fast5_dir" && -r "$fast5_dir/F1-1-paired-mini.fast5" ]]
case "$example_dir" in
    *$'\t'*|*$'\n'*|*$'\r'*|*'"'*)
        printf 'The example directory must not contain tabs, line breaks, or double quotes.\n' >&2
        exit 2
        ;;
esac
printf 'F1-1\tunused\t%s\t%s\n' "$fastq" "$fast5_dir"
