#!/usr/bin/env bash
#SBATCH --job-name=hic-pretext
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=48
#SBATCH --time=24:00:00
#SBATCH --mem=200GB
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com


# AutoHiC: Align Hi-C reads to the AutoHiC-corrected assembly
# and generate a PretextMap file for manual curation in PretextView.

cd /hpcfs/users/a1864358/sanders_lab/asm/files/horn

FA="horn-morehic.fa"  
HIC_R1="SRR16961052_1.fastq.gz"
HIC_R2="SRR16961052_2.fastq.gz"
OUT="/hpcfs/users/a1864358/sanders_lab/asm/files/horn/pre-manCur"
PRETEXT_PREFIX="horn-premancur"                      

set -euo pipefail
module purge

export PATH="/scratchdata1/users/a1864358/miniconda/bin:$PATH"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /hpcfs/users/a1864358/miniconda/envs/general

THREADS=${SLURM_CPUS_PER_TASK}
mkdir -p "${OUT}"

bwa index "${FA}"

bwa mem -5SP -t "${THREADS}" "${FA}" "${HIC_R1}" "${HIC_R2}" \
    | samblaster --ignoreUnmated \
    | samtools view -@ "${THREADS}" -h -F 2316 \
    | PretextMap \
        -o "${OUT}/${PRETEXT_PREFIX}.pretext" \
        --sortby length \
        --sortorder descend \
        --mapq 10
