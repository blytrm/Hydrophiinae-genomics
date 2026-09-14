#!/bin/bash
#SBATCH -J asm_qc
#SBATCH -a 0-3                 
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH -o /hpcfs/users/a1864358/sanders_lab/resources/terstrl-elapids_liftover-anno/logs/asm_qc_%A_%a.out

# 0=nscut_old 1=nscut_new 2=ptext_old 3=ptext_new
#   => array job

# Assembly QC on OLD (2018 RefSeq) and NEW (2025) assemblies, both species.
#   compleasm  
#   quast      

set -euo pipefail

RES=/hpcfs/users/a1864358/sanders_lab/resources
OUT=$RES/terstrl-elapids_liftover-anno
DB=/hpcfs/users/a1864358/dbs/compleasm/mb_downloads     
T=${SLURM_CPUS_PER_TASK:-16}

SP=(nscut nscut ptext ptext)
ERA=(old   new   old   new)
ASM=(
  $RES/nscut/nscut.fa
  $RES/nscut/2025/nscut_2025.fa
  $RES/ptext/ptext.fa
  $RES/ptext/2025/ptext_2025.fa
)

i=$SLURM_ARRAY_TASK_ID
tag=${SP[$i]}_${ERA[$i]}
asm=${ASM[$i]}
o=$OUT/qc/$tag
mkdir -p "$o"

echo "[$(date)] QC $tag -> $asm"

source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh

# compleasm
set +u; conda activate compleasm; set -u
compleasm run -a "$asm" -o "$o/compleasm" -l squamata -L "$DB" -t "$T"

# squamata
set +u; conda activate lifton; set -u
quast.py "$asm" -o "$o/quast" --threads "$T" --eukaryote --plots-format png --no-icarus

echo "[$(date)] done $tag"
