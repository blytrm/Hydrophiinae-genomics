#!/usr/bin/env bash
#SBATCH -J postEviann
#SBATCH -p newicelake
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --time=2-00:00:00
#SBATCH -o /hpcfs/users/a1864358/sanders_lab/asm/files/annotation/nf/logs/%x_%j.out
#SBATCH -e /hpcfs/users/a1864358/sanders_lab/asm/files/annotation/nf/logs/%x_%j.err

set -euo pipefail
export PATH="/hpcfs/users/a1864358/miniconda/miniconda3/bin:$PATH"
source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh

cd /hpcfs/users/a1864358/sanders_lab/asm/files/annotation/nf
bash run_hmaj.sh
