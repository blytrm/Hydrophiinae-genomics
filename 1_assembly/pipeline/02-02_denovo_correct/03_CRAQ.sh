#!/usr/bin/env bash
#SBATCH --job-name=craq-hcy
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
conda activate /hpcfs/users/a1864358/miniconda/envs/general

cd /hpcfs/users/a1864358/sanders_lab/asm/files/hcy

CRAQ='/hpcfs/users/a1864358/miniconda/envs/general/bin/craq'
ASM='/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/ins_corrected.fa'
SHORT_READ_1='/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/hcy-sr-1.fastq.gz'
SHORT_READ_2='/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/hcy-sr-2.fastq.gz'
OUT='/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/craq'
CLR='/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/all_clr.fastq.gz'

# short-read alignment
bwa index "${ASM}"
bwa mem -5SP -t ${SLURM_CPUS_PER_TASK} "${ASM}" \
    "${SHORT_READ_1}" "${SHORT_READ_2}" | \
    samblaster | \
    samtools view - -@ ${SLURM_CPUS_PER_TASK} -S -h -b -F 3340 -o sr-hcy.unsorted.bam

samtools sort -@ ${SLURM_CPUS_PER_TASK} -m 4G sr-hcy.unsorted.bam -o sr-hcy.sorted.bam
samtools index sr-hcy.sorted.bam
rm sr-hcy.unsorted.bam

# long-read alignment
minimap2 -ax map-pb -t ${SLURM_CPUS_PER_TASK} "${ASM}" "${CLR}" | \
    samtools sort -@ ${SLURM_CPUS_PER_TASK} -m 4G -o clr-hcy.sorted.bam
samtools index clr-hcy.sorted.bam

"${CRAQ}" \
    --genome "${ASM}" \
    --sms_input clr-hcy.sorted.bam \
    --ngs_input sr-hcy.sorted.bam \
    --map map-pb \
    --break T \
    --plot T \
    --output_dir "${OUT}" \
    --thread "${SLURM_CPUS_PER_TASK}"


