#!/usr/bin/env bash
#SBATCH --job-name=hic-realign
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=40
#SBATCH --time=20:00:00
#SBATCH --mem=210GB
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com

ASM="/hpcfs/users/a1864358/sanders_lab/asm/files/horn/horn-haphic.fa"             
HIC_R1="/hpcfs/users/a1864358/sanders_lab/asm/files/horn/SRR16961052_1.fastq.gz"       
HIC_R2="/hpcfs/users/a1864358/sanders_lab/asm/files/horn/SRR16961052_2.fastq.gz"       
OUT_BAM="/hpcfs/users/a1864358/sanders_lab/asm/files/horn/horn-realign.bam"

set -euo pipefail
module purge

export PATH="/scratchdata1/users/a1864358/miniconda/bin:$PATH"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /hpcfs/users/a1864358/miniconda/envs/general

THREADS=${SLURM_CPUS_PER_TASK}
mkdir -p "$(dirname "${OUT_BAM}")"

bwa index "${ASM}"

bwa mem -5SP -t "${THREADS}" "${ASM}" "${HIC_R1}" "${HIC_R2}" | \
    samblaster | \
    samtools view - -@ "${THREADS}" -b -F 3340 \
        -o "${OUT_BAM}"
