#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    printf 'Usage: %s EXAMPLE_DIR\n' "$0" >&2
    exit 2
fi
example_dir=$(cd -- "$1" && pwd -P)
spectrum=M2brain-DDA-mini.mzML
if [[ ! -r "$example_dir/$spectrum" || ! -s "$example_dir/$spectrum" ]]; then
    printf 'Expected the indexed mzML in the example directory.\n' >&2
    exit 1
fi
# Save stdout beside the mzML. Nextflow resolves this path from the manifest directory.
printf '%s\t%s\t\tDDA\n' "$spectrum" M2brain-mini
