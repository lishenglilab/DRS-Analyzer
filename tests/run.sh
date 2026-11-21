#!/usr/bin/env bash
set -euo pipefail

cd -- "$(dirname -- "$0")/.."
python3 -E -B tests/test_helpers.py
python3 -E -B tests/test_proteomics.py
exec python3 -E -B tests/test_fragpipe231.py
