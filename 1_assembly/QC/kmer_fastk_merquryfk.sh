#!/usr/bin/env bash
#SBATCH --job-name=FastK-MerquryFK
#SBATCH -p batch
#SBATCH -N 1
#SBATCH --cpus-per-task=20
#SBATCH --time=04:00:00
#SBATCH --mem=60GB
#SBATCH -o ./joblogs/%x_%j.log

set -euo pipefail

KMER=31

WDIR="/home/${USER}/al-biohub/projects/billy/honours/01_genome_assembly/<species>"
WGS_DIR="/home/${USER}/al-biohub/projects/billy/honours/data/sra-data/genome/<species>/illumina_wgs"

READ_R1="${WGS_DIR}/<R1>.fastq.gz"
READ_R2="${WGS_DIR}/<R2>.fastq.gz"

DB_PREFIX="${WDIR}/merqury/reads/<species>_short-read"
mkdir -p "$(dirname "${DB_PREFIX}")"

TMP_DIR="${WDIR}/tmp"
mkdir -p "${TMP_DIR}"

echo -e "[$(date)] FastK start\n\treads: ${READ_R1} ${READ_R2}\n\tk=${KMER}\n\tthreads=${SLURM_CPUS_PER_TASK}\n\ttmp=${TMP_DIR}"

# fastk v1.2
FastK \
    -v \
    -k${KMER} \
    -t1 \
    -T"${SLURM_CPUS_PER_TASK}" \
    -M50 \
    -P"${TMP_DIR}" \
    -N"${DB_PREFIX}" \
    "${READ_R1}" "${READ_R2}"

# merquryFK v1.2
MerquryFK \
    -v \
    -T${SLURM_CPUS_PER_TASK} \
    -P"${TMP_DIR}" \
    "${DB_PREFIX}" "${WDIR}/<species>.fa" "${WDIR}/<species>"
