#!/usr/bin/env bash
#SBATCH --job-name=hcurw-gapfill-polish
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --time=3-00:00:00
#SBATCH --mem=750GB
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com

set -euo pipefail

WORKDIR="/hpcfs/users/a1864358/sanders_lab/asm/files/hcur-w"
ASM="${WORKDIR}/hcurw-mancur.fa"
LONGREADS="${WORKDIR}/hcu_west.fastq.gz"
SR1="${WORKDIR}/SRR16961055_1.fastq.gz"
SR2="${WORKDIR}/SRR16961055_2.fastq.gz"

THREADS="${SLURM_CPUS_PER_TASK:-32}"
MIN_READ_LEN=30000
GAPFILL_FA="${WORKDIR}/hcurw-gapfill.fa"
GAPFILL_DIR="${WORKDIR}/gapfill/${SLURM_JOB_ID:-manual_$(date '+%Y%m%d_%H%M%S')}"
POLISH_DIR="${WORKDIR}/nextpolish-final"

module purge
export PATH="/hpcfs/users/a1864358/miniconda/envs/general/bin:$PATH"
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u; conda activate /hpcfs/users/a1864358/miniconda/envs/general; set -u

mkdir -p "${GAPFILL_DIR}" && cd "${GAPFILL_DIR}"

# gapfill w longreads
echo "[$(date)] === Step 1: Gapfilling ==="

# Length-filter ONT reads 
LONG_FILT="${GAPFILL_DIR}/reads.${MIN_READ_LEN}.fa"
seqkit seq -j "${THREADS}" -m "${MIN_READ_LEN}" -o "${LONG_FILT}" "${LONGREADS}"
seqkit stats -T "${LONG_FILT}"

tgsgapcloser \
    --scaff   "${ASM}" \
    --reads   "${LONG_FILT}" \
    --output  "${GAPFILL_DIR}/hcurw-tgs" \
    --tgstype ont \
    --ne \
    --thread  "${THREADS}" \
    --minmap_arg '-ax map-ont' \
    >"${GAPFILL_DIR}/tgsgapcloser.out" \
    2>"${GAPFILL_DIR}/tgsgapcloser.err"

cp "${GAPFILL_DIR}/hcurw-tgs.scaff_seqs" "${GAPFILL_FA}"
seqkit stats -T "${ASM}" "${GAPFILL_FA}"
echo "[$(date)] Step 1 done: ${GAPFILL_FA}"



# nextpolish w shortreads
echo "[$(date)] === Step 2: Polishing ==="

mkdir -p "${POLISH_DIR}" && cd "${POLISH_DIR}"
printf '%s\n%s\n' "${SR1}" "${SR2}" > "${POLISH_DIR}/sgs.fofn"

cat > "${POLISH_DIR}/run.cfg" <<EOF
[General]
job_type = local
job_prefix = hcurw_nextpolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 4
multithread_jobs = 8
genome = ${GAPFILL_FA}
genome_size = auto
workdir = ${POLISH_DIR}/rundir
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ${POLISH_DIR}/sgs.fofn
sgs_options = -max_depth 100 -bwa
EOF

nextPolish "${POLISH_DIR}/run.cfg"

FINAL_FA="${POLISH_DIR}/rundir/genome.nextpolish.fasta"
seqkit stats -T "${GAPFILL_FA}" "${FINAL_FA}"
echo "[$(date)] Step 2 done: ${FINAL_FA}"
