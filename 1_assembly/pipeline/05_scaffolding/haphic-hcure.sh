#!/usr/bin/env bash
#SBATCH --job-name=haphic
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=48
#SBATCH --time=24:00:00
#SBATCH --mem=240GB
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com


ASM="/hpcfs/users/a1864358/sanders_lab/asm/files/hcur-E/hcure-polish1.fa"              
HIC_R1="/hpcfs/users/a1864358/sanders_lab/asm/files/hcur-E/hcurE-hic_1.fastq.gz"       
HIC_R2="/hpcfs/users/a1864358/sanders_lab/asm/files/hcur-E/hcurE-hic_2.fastq.gz"       
HIC_BAM_PREFIX="HiC-hcure"           
N_CHROM=14                              
RE="GATC,GANTC,CTNAG,TTAAC"           
CORRECT_NROUNDS=3                       
NPOP=200
MUTPROB=0.15
OUT="/hpcfs/users/a1864358/sanders_lab/asm/files/hcur-E/haphic_out"
WORKDIR="/hpcfs/users/a1864358/sanders_lab/asm/files/hcur-E"
HAPHIC="/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/HapHiC"


set -euo pipefail
module purge

export PATH="/scratchdata1/users/a1864358/miniconda/bin:$PATH"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /hpcfs/users/a1864358/miniconda/envs/haphic

cd "${WORKDIR}"

bwa index "${ASM}"

bwa mem -5SP -t "${SLURM_CPUS_PER_TASK}" "${ASM}" "${HIC_R1}" "${HIC_R2}" | \
    samblaster | \
    samtools view - -@ "${SLURM_CPUS_PER_TASK}" -S -h -b -F 3340 \
        -o "${HIC_BAM_PREFIX}.bam"

"${HAPHIC}/utils/filter_bam" "${HIC_BAM_PREFIX}.bam" 1 \
    --nm 3 \
    --threads "${SLURM_CPUS_PER_TASK}" | \
    samtools view - -b -@ "${SLURM_CPUS_PER_TASK}" \
        -o "${HIC_BAM_PREFIX}.filtered.bam"

"${HAPHIC}/haphic" pipeline \
    "${ASM}" "${HIC_BAM_PREFIX}.filtered.bam" "${N_CHROM}" \
    --RE "${RE}" \
    --outdir "${OUT}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --processes "${SLURM_CPUS_PER_TASK}" \
    --correct_nrounds "${CORRECT_NROUNDS}" \
    --npop "${NPOP}" \
    --mutprob "${MUTPROB}"
