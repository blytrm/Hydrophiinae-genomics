#!/usr/bin/env bash
# One-off: install eggnog-mapper into its own env and fetch its database.
set -euo pipefail
DATA=/scratchdata1/users/a1864358/sanders_lab/prog/eggnog_db
mkdir -p "$DATA"
source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
conda env list | grep -q '^eggnog ' || \
    conda create -y -n eggnog -c conda-forge -c bioconda eggnog-mapper
conda activate eggnog
export EGGNOG_DATA_DIR="$DATA"
# -y accepts the download prompts; -M also fetches the MMseqs/diamond DB
[[ -s "$DATA/eggnog.db" ]] || download_eggnog_data.py -y --data_dir "$DATA"
echo "EGGNOG READY: $(du -sh "$DATA" | cut -f1) in $DATA"
