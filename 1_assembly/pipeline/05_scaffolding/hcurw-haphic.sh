#!/usr/bin/env bash
#SBATCH --job-name=haphic-hcurw
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=48
#SBATCH --time=24:00:00
#SBATCH --mem=240GB
#SBATCH -o %x_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com

module purge

export PATH="/scratchdata1/users/a1864358/miniconda/bin:$PATH"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /hpcfs/users/a1864358/miniconda/envs/haphic

cd /hpcfs/users/a1864358/sanders_lab/asm/files/hcur-w

ASM='/hpcfs/users/a1864358/sanders_lab/asm/files/hcur-w/hcurw-craq.fa'
HIC1='/hpcfs/users/a1864358/sanders_lab/asm/files/hcur-w/SRR16961053_1.fastq.gz'
HIC2='/hpcfs/users/a1864358/sanders_lab/asm/files/hcur-w/SRR16961053_2.fastq.gz'

bwa index "${ASM}"
bwa mem -5SP -t "${SLURM_CPUS_PER_TASK}" "${ASM}" "${HIC1}" "${HIC2}" | samblaster | samtools view - -@ "${SLURM_CPUS_PER_TASK}" -S -h -b -F 3340 -o HiC-hcurw.bam

/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/HapHiC/utils/filter_bam HiC-hcurw.bam 1 --nm 3 --threads "${SLURM_CPUS_PER_TASK}" | samtools view - -b -@ "${SLURM_CPUS_PER_TASK}" -o HiC-hcurw.filtered.bam

BAM="/hpcfs/users/a1864358/sanders_lab/asm/files/hcur-w/HiC-hcurw.filtered.bam"
OUT="/hpcfs/users/a1864358/sanders_lab/asm/files/hcur-w/haphic"
SW="/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/HapHiC"

${SW}/haphic pipeline \
    "${ASM}" "${BAM}" 14 \
    --RE "GATC,GANTC,CTNAG,TTAAC" \
    --outdir ${OUT} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --processes ${SLURM_CPUS_PER_TASK} \
    --correct_nrounds 4 \
    --npop 200 \
    --mutprob 0.15
