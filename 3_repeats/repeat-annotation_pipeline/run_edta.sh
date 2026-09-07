#!/bin/bash
#SBATCH --job-name=edta_10horn
#SBATCH --partition=dm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180G
#SBATCH --time=8-00:00:00
#SBATCH --output=edta_%j.log
#SBATCH --error=edta_%j.err

set -euo pipefail

source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
conda activate EDTA

PROJ=/scratchdata1/users/a1864358/sanders_lab/repeatdensity-atv2rs
GENOME=$PROJ/10-horn-final.renamed.fa
CURATED=$PROJ/repeats/hydrophis_ornatus-garvin.fa.mod.EDTA.TElib.fa
THREADS=36
RUNDIR=$PROJ/edta_run

mkdir -p "$RUNDIR"
cd "$RUNDIR"
ln -sf "$GENOME" .
GBASE=$(basename "$GENOME")

EDTA.pl \
  --genome "$GBASE" \
  --species others \
  --step all \
  --curatedlib "$CURATED" \
  --sensitive 1 \
  --anno 1 \
  --threads "$THREADS" \
  --overwrite 0

