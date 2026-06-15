#!/usr/bin/env bash
#SBATCH --job-name=mosdepth
#SBATCH -p batch
#SBATCH -N 1
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00
#SBATCH --mem=80GB
#SBATCH -o ./joblogs/%x_%j.log

# mosdepth v0.3.11
# minmap2 v2.31-r1302
# samtools v1.20

set -euo pipefail

DIR="/home/${USER}/al-biohub/projects/billy/honours"
WDIR="${DIR}/01_genome_assembly/<species>/final"
ASM="${WDIR}/<species>.fa"
READS="${WDIR}/corrected_nanopore.fasta"
READS="${DIR}/data/sra-data/genome/<species>/longread_wgs/<species>.fastq.gz"

cd "${WDIR}" || exit 1

# Align long-reads (corrected reads where possible)
minimap2 \
    -a \
    -x map-ont \
    -L \
    -t 18 \
    -2 \
    -K 2G \
    --secondary=no \
    "${ASM}" "${READS}" |
samtools sort \
    -@ 2 \
    -m 2G \
    -O BAM \
    -o "<species>.longread.bam" -

samtools index -@ 16 "<species>.longread.bam"

# Depth
mkdir mosdepth
mosdepth \
    -t ${SLURM_CPUS_PER_TASK} \
    --mapq 20 \
    "mosdepth/<species>.longread" \
    "<species>.longread.bam"
