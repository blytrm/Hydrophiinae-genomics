#!/bin/bash
#SBATCH -J v2r_tree
#SBATCH -p icelake
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH -t 3-00:00:00
#SBATCH -o /scratchdata1/users/a1864358/sanders_lab/phylogenomics/v2r_tree/workspace/logs/30_tree_%j.out
set -euo pipefail
PROJ=/scratchdata1/users/a1864358/sanders_lab/phylogenomics/v2r_tree
D=$PROJ/workspace/data
R=$PROJ/results
mkdir -p "$R"
source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
T=${SLURM_CPUS_PER_TASK:-32}

# 1. exploratory tree (minutes) — FastTree
set +u; conda activate general; set -u
echo "[$(date)] FastTree exploratory"
FastTreeMP -lg -gamma < "$D/backbone.trim.aln" > "$R/mega.fasttree.nwk"
echo "[$(date)] fasttree tips: $(grep -o ',' "$R/mega.fasttree.nwk" | wc -l)"

# 2. definitive ML — IQ-TREE (LG+F+G4 + UFBoot; prototype-proven)
set +u; conda activate of3_env; set -u
echo "[$(date)] IQ-TREE ML (LG+F+G4, UFBoot 1000, --bnni)"
iqtree -s "$D/backbone.trim.aln" --prefix "$R/mega.ml" \
    -m LG+F+G4 -B 1000 --bnni -T AUTO --ntmax "$T" --redo
echo "[$(date)] ML tree: $R/mega.ml.treefile"
echo "[$(date)] ALL DONE — tree ready"
