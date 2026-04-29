#!/usr/bin/env bash
#SBATCH --job-name=haphic-hcy
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

cd /hpcfs/users/a1864358/sanders_lab/asm/files/hcy

bwa index nd.asm.fasta
bwa mem -5SP -t ${SLURM_CPUS_PER_TASK} nd.asm.fasta ./SRR10692938_1.fastq.gz ./SRR10692938_1.fastq.gz | samblaster | samtools view - -@ 14 -S -h -b -F 3340 -o HiC-hmaj.bam

/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/HapHiC/utils/filter_bam HiC-hmaj.bam 1 --nm 3 --threads 20 | samtools view - -b -@ 14 -o HiC-hmaj.filtered.bam

FA="/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/nd.asm.fasta"
BAM="/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/HiC-hcy.filtered.bam"
OUT="/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/haphic"
SW="/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/HapHiC"

${SW}/haphic pipeline \
    ${FA} ${BAM} 14 \
    --RE "GATC,GANTC,CTNAG,TTAAC" \
    --outdir ${OUT} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --processes ${SLURM_CPUS_PER_TASK} \
    --correct_nrounds 3
