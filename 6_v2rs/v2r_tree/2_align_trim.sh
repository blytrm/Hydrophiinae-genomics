#!/bin/bash
#SBATCH -J v2r_aln
#SBATCH -p icelake
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o /scratchdata1/users/a1864358/sanders_lab/phylogenomics/v2r_tree/workspace/logs/20_aln_%j.out
set -euo pipefail
PROJ=/scratchdata1/users/a1864358/sanders_lab/phylogenomics/v2r_tree
D=$PROJ/workspace/data
source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
set +u; conda activate general; set -u
T=${SLURM_CPUS_PER_TASK:-32}

echo "[$(date)] MAFFT align $(grep -c '^>' "$D/backbone.faa") seqs"

mafft --auto --anysymbol --thread "$T" "$D/backbone.faa" > "$D/backbone.aln"
echo "[$(date)] aligned: $(grep -c '^>' "$D/backbone.aln") seqs, cols=$(awk '/^>/{if(s)exit;s=1;next}{w+=length($0)}END{print w}' "$D/backbone.aln")"

echo "[$(date)] trimal -automated1"
trimal -in "$D/backbone.aln" -out "$D/backbone.trim.aln" -automated1
COLS=$(awk '/^>/{if(s)exit;s=1;next}{w+=length($0)}END{print w}' "$D/backbone.trim.aln")
echo "[$(date)] trimmed cols=$COLS"

# guard: if -automated1 over-trims below the 200-col floor, fall back to gappy-column removal
if [ "${COLS:-0}" -lt 200 ]; then
  echo "[$(date)] WARN trimmed to $COLS cols (<200); refalling back to -gt 0.1"
  trimal -in "$D/backbone.aln" -out "$D/backbone.trim.aln" -gt 0.1
  echo "[$(date)] fallback trimmed cols=$(awk '/^>/{if(s)exit;s=1;next}{w+=length($0)}END{print w}' "$D/backbone.trim.aln")"
fi
echo "[$(date)] DONE"
