#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem 32G
#SBATCH --time 2-00:00:00
#SBATCH -J generax_v2r
#SBATCH -o /hpcfs/users/a1864358/sanders_lab/phylogenomics/cafe/generax/generax.%j.log

# GeneRax 2.1.3 returns a zero phylogenetic likelihood on this family.
# 2.0.4 reconciles the intact V2R ML tree (2436 tips, 948 columns) under UndatedDL.
# Radius-1 SPR [edit of gene tree ... 1: those edits are limited a neighbouring branch] --> optimises duplication/loss rates and the reconciliation root.

set -euo pipefail
export OMP_NUM_THREADS=12
export OMPI_MCA_btl=^openib
source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
conda activate generax2
cd /hpcfs/users/a1864358/sanders_lab/phylogenomics/cafe/generax

generax \
  -f families.txt \
  -s /hpcfs/users/a1864358/sanders_lab/phylogenomics/cafe/trees/species_tree_ultrametric.nwk \
  -r UndatedDL \
  --strategy SPR \
  --max-spr-radius 1 \
  --per-species-rates \
  -p out_v2r \
  --seed 123
  
echo DONE $(date)
