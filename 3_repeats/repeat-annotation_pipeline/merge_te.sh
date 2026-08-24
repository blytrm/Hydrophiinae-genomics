#!/usr/bin/env bash
#SBATCH --job-name=merge_te
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --chdir=/hpcfs/users/a1864358/sanders_lab/asm/files/repeat-anno/ra2
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail

BASE=/hpcfs/users/a1864358/sanders_lab/asm/files/repeat-anno/ra2
SAMPLE=SRR16961054
OUTDIR="$BASE/results/04_merged"

mkdir -p logs "$OUTDIR"

source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
conda run --no-capture-output -n TEforest python "$BASE/merge_te_annotations.py" \
    --teforest "$BASE/results/03_teforest/output/${SAMPLE}_TEforest_bps_nonredundant.bed" \
    --ref_bed "$BASE/results/reference_TEs_teforest.bed" \
    --window 100 \
    --out_gff3 "$OUTDIR/TE_annotation_unified.gff3" \
    --out_bed "$OUTDIR/TE_annotation_unified.bed"
