#!/usr/bin/env bash
#SBATCH --job-name=of_prot_lift
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --time=48:00:00
#SBATCH --mem=192GB
#SBATCH -o ./%x_%j.log
# Proteome OrthoFinder, 12 taxa, NO cadam, lifted nscut/ptext. Inputs pre-collapsed.
set -euo pipefail
source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
set +u; conda activate of3_env; set -u
ROOT=/hpcfs/users/a1864358/sanders_lab/phylogenomics
T=${SLURM_CPUS_PER_TASK:-32}
orthofinder -t "$T" -a "$T" -f "$ROOT/orthofinder/primary_transcripts_lifted" -n lifted
