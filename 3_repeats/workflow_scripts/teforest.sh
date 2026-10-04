#!/usr/bin/env bash
#SBATCH --job-name=teforest
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --time=3-00:00:00
#SBATCH --mem=200G
#SBATCH --chdir=/hpcfs/users/a1864358/sanders_lab/asm/files/repeat-anno/ra2
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail

BASE=/hpcfs/users/a1864358/sanders_lab/asm/files/repeat-anno/ra2
SAMPLE=SRR16961054
OUTDIR="$BASE/results/03_teforest"

mkdir -p logs "$OUTDIR"

source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
conda run --no-capture-output -n TEforest python /hpcfs/users/a1864358/sanders_lab/prog/TEforest/TEforest.py \
    --workflow_dir /hpcfs/users/a1864358/sanders_lab/prog/TEforest/workflow \
    --workdir "$OUTDIR" \
    --threads "${SLURM_CPUS_PER_TASK:-72}" \
    --samples "$SAMPLE" \
    --fq_base_path "$BASE/data" \
    --consensusTEs "$BASE/results/hybridTEs.fasta" \
    --ref_genome "$BASE/data/10-horn-final.renamed.fa" \
    --ref_te_locations "$BASE/results/reference_TEs_teforest.bed" \
    --euchromatin "$BASE/results/euchromatin_wholegenome.bed" \
    --model "$BASE/models/TEforest_Nonreference50X.pkl" \
    --ref_model "$BASE/models/TEforest_reference50X.pkl" \
    --cleanup_intermediates
