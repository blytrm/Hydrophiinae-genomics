#!/bin/bash
#SBATCH --job-name=agpcorrect-hcy
#SBATCH --output=agpcorrect_%j.out
#SBATCH --error=agpcorrect_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00

set -euo pipefail

cd /hpcfs/users/a1864358/sanders_lab/asm/files/hcy/files/hcy.backup0.backup0

source /hpcfs/users/a1864358/miniconda/etc/profile.d/conda.sh
conda activate general

python3 AGPCorrect \
    hcy-autocur.fa \
    hcy-curated-pretext.agp \
    > hcy-curated-corrected.agp

echo "Done. Output: hcy-curated-corrected.agp"
