#!/usr/bin/env bash
#SBATCH --job-name=nextpolish
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH --cpus-per-task=40
#SBATCH --time=3-00:00:00
#SBATCH --mem=400GB
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com

ASM="/scratchdata1/users/a1864358/sanders_lab/asm/files/hcur-E/hcure-ins.fa"
CLR="/scratchdata1/users/a1864358/sanders_lab/asm/files/hcur-E/corrected.fasta.gz"
GENOME_SIZE="1900000000"
LGS_MAP="-x map-pb"                     
WORKDIR="/scratchdata1/users/a1864358/sanders_lab/asm/files/hcur-E"
OUTDIR="${WORKDIR}/polish"

set -euo pipefail
module purge

export PATH="/scratchdata1/users/a1864358/miniconda/bin:$PATH"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /hpcfs/users/a1864358/miniconda/envs/general

export TMPDIR="${TMPDIR:-/tmp}/${SLURM_JOB_ID:-$$}"
mkdir -p "${TMPDIR}" "${OUTDIR}"

cd "${OUTDIR}"

LGS_FOFN="${OUTDIR}/lgs.fofn"
echo "${CLR}" > "${LGS_FOFN}"

RUNDIR="${OUTDIR}/rundir.${SLURM_JOB_ID:-$$}"
CONFIG="${OUTDIR}/nextpolish-run.${SLURM_JOB_ID:-$$}.cfg"

cat > "${CONFIG}" <<EOF
[General]
job_type = local
job_prefix = nextPolish
task = 5566
rewrite = yes
rerun = 3
parallel_jobs = 4
multithread_jobs = 10
genome = ${ASM}
genome_size = ${GENOME_SIZE}
workdir = ${RUNDIR}
polish_options = -p 10

[lgs_option]
lgs_fofn = ${LGS_FOFN}
lgs_options = -min_read_len 1k -max_depth 100
lgs_minimap2_options = ${LGS_MAP} -t 10
EOF

nextPolish "${CONFIG}"
