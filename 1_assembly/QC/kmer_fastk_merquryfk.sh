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

WDIR="/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hmaj"

READ_R1="${WDIR}/hydmaj.fastq.gz"

DB_PREFIX="${WDIR}/merqury/reads/maj"
mkdir -p "$(dirname "${DB_PREFIX}")"

TMP_DIR="${WDIR}/tmp"
mkdir -p "${TMP_DIR}"

echo -e "[$(date)] FastK start\n\treads: ${READ_R1} ${READ_R2}\n\tk=${KMER}\n\tthreads=${SLURM_CPUS_PER_TASK}\n\ttmp=${TMP_DIR}"

FastK \
    -v \
    -k${KMER} \
    -t1 \
    -T"${SLURM_CPUS_PER_TASK}" \
    -M50 \
    -P"${TMP_DIR}" \
    -N"${DB_PREFIX}" \
    "${READ_R1}"

MerquryFK \
    -v \
    -T${SLURM_CPUS_PER_TASK} \
    -P"${TMP_DIR}" \
    "${DB_PREFIX}" "${WDIR}/hmaj_filt_denovo.fa" "${WDIR}/hmaj"
