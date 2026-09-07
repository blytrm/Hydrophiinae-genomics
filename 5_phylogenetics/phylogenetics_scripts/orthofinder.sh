#!/usr/bin/env bash
#SBATCH --job-name=orthofinder
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --time=24:00:00
#SBATCH --mem=128GB
#SBATCH -o ./%x_%j.log

set -euo pipefail
source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
set +u
conda activate of3_env
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
PROTEOMES="${HERE}/proteomes"
PRIMARY="${HERE}/primary_transcripts"

mkdir -p "${PRIMARY}"

for f in "${PROTEOMES}"/*.fa; do
    sp="$(basename "${f}" .fa)"
    python3 "${HERE}/scripts/longest_isoform.py" "${f}" "${PRIMARY}/${sp}.fa"
done

orthofinder -t "${THREADS}" -a "${THREADS}" -f "${PRIMARY}"
