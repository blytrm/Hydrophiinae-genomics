#!/bin/bash
#SBATCH --job-name=edta_array
#SBATCH --partition=dm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180G
#SBATCH --time=2-00:00:00
#SBATCH --array=1-4
#SBATCH --output=edta_%A_%a.log
#SBATCH --error=edta_%A_%a.err

set -euo pipefail

source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
conda activate EDTA

ASM_DIR=/hpcfs/users/a1864358/sanders_lab/asm/files/final-asms
CURATED=/scratchdata1/users/a1864358/sanders_lab/repeatdensity-atv2rs/resources-repeats/hydrophis_ornatus-garvin.fa.mod.EDTA.TElib.fa
SCRATCH=/scratchdata1/users/a1864358/sanders_lab/asm/files/repeat-anno
THREADS=36

SAMPLES=(dummy 10-hcure 10-hcurw 10-hcy 10-hmaj)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
GENOME=$ASM_DIR/${SAMPLE}-final.renamed.fa
RUNDIR=$SCRATCH/edta_${SAMPLE}

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
