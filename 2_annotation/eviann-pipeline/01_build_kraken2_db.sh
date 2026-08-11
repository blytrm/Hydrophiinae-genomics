#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

source "${CONDA_PROFILE}"
conda activate "${RNAQC_ENV}"

if [[ -d "${KRAKEN2_DB}/taxonomy" && -s "${KRAKEN2_DB}/hash.k2d" ]]; then
    echo "[$(date)] Kraken2 DB already built at ${KRAKEN2_DB}; nothing to do."
    exit 0
fi

mkdir -p "${KRAKEN2_DB}"

kraken2-build --download-taxonomy --skip-maps --use-ftp --db "${KRAKEN2_DB}"

kraken2-build --download-library UniVec_Core --use-ftp --no-masking --db "${KRAKEN2_DB}"

kraken2-build --build \
    --threads 8 \
    --db "${KRAKEN2_DB}" \
    --no-masking

echo "[$(date)] Kraken2 DB build complete: ${KRAKEN2_DB}"
