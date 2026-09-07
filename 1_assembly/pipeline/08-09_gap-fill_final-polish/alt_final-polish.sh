#!/usr/bin/env bash
#SBATCH --job-name=gap-pol
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --time=2-00:00:00
#SBATCH --mem=750GB
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com


set -euo pipefail

ASM=""
SR1="/hpcfs/users/a1864358/sanders_lab/asm/files/horn/SRR10617472_1.fastq.gz"
SR2="/hpcfs/users/a1864358/sanders_lab/asm/files/horn/SRR10617472_2.fastq.gz"

THREADS="${SLURM_CPUS_PER_TASK:-32}"

mkdir -p "${OUT}"/{logs,gapfill,polish}
cd "${OUT}"

# 2. POLISH (NextPolish, short reads) -- skip if no SR provided
FINAL="${GAPFILLED}"
if [[ -n "${SR1}" && -n "${SR2}" ]]; then
  PDIR="${OUT}/polish"
  printf '%s\n%s\n' "${SR1}" "${SR2}" > "${PDIR}/sgs.fofn"
  cat > "${PDIR}/run.cfg" <<EOF
[General]
job_type = local
job_prefix = nextpolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 2
multithread_jobs = 16
genome = ${GAPFILLED}
genome_size = auto
workdir = ${PDIR}/rundir
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ${PDIR}/sgs.fofn
sgs_options = -max_depth 100 -bwa
EOF
  nextPolish "${PDIR}/run.cfg"
  FINAL="${PDIR}/rundir/genome.nextpolish.fasta"
fi

seqkit stats -T "${ASM}" "${GAPFILLED}" "${FINAL}" | tee "${OUT}/logs/final.stats.tsv"
echo "[$(date)] DONE -> ${FINAL}"
