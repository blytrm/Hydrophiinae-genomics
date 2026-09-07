#!/usr/bin/env bash
#SBATCH --job-name=haphic-hcy
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=40  
#SBATCH --time=16:00:00
#SBATCH --mem=220GB
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com

set -euo pipefail

module purge

export PATH="/scratchdata1/users/a1864358/miniconda/bin:$PATH"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /hpcfs/users/a1864358/miniconda/envs/haphic

cd /hpcfs/users/a1864358/sanders_lab/asm/files/hcy/files/hcy.backup0.backup0

FA="/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/files/hcy.backup0.backup0/3-dnqc/hcy-polish1.fasta"
BAM="/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/files/hcy.backup0.backup0/HiC-hcy.bam"
OUT="/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/files/hcy.backup0.backup0/haphic-2"
SW="/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/HapHiC"
HIC_R1="/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/files/hcy.backup0.backup0/SRR10692938_1.fastq.gz"
HIC_R2="/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/files/hcy.backup0.backup0/SRR10692938_2.fastq.gz"

bwa index "${FA}"
bwa mem -5SP -t "${SLURM_CPUS_PER_TASK}" "${FA}" "${HIC_R1}" "${HIC_R2}" | samblaster | samtools view - -@ "${SLURM_CPUS_PER_TASK}" -S -h -b -F 3340 -o hic-int.bam

/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/HapHiC/utils/filter_bam hic-int.bam 1 --nm 3 --threads "${SLURM_CPUS_PER_TASK}" | samtools view - -b -@ "${SLURM_CPUS_PER_TASK}" -o "${BAM}"

${SW}/haphic pipeline \
    "${FA}" "${BAM}" 14 \
    --RE "GATC" \
    --outdir "${OUT}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --processes "${SLURM_CPUS_PER_TASK}" \
    --correct_nrounds 2 \
    --bin_size 200 \
    --topN 25 \
    --mutprob 0.15 \
    --npop 200 \
    --min_inflation 1.8 \
    --max_inflation 10.0 
