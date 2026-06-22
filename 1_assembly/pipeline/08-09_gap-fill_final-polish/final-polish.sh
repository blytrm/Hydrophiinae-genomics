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

module purge
export PATH="/scratchdata1/users/a1864358/miniconda/bin:$PATH"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /hpcfs/users/a1864358/miniconda/envs/general

ASM="/scratchdata1/users/a1864358/sanders_lab/asm/files/hcur-E/hcure-curated.fa"
OUT="/scratchdata1/users/a1864358/sanders_lab/asm/files/hcur-E/polish-final"
SR1="/scratchdata1/users/a1864358/sanders_lab/asm/files/hcur-E/SRR11461820_1.fastq.gz"
SR2="/scratchdata1/users/a1864358/sanders_lab/asm/files/hcur-E/SRR11461820_2.fastq.gz"

TARGET_COV=50          # short-read coverage to subsample to (rasusa)
THREADS="${SLURM_CPUS_PER_TASK:-32}"

mkdir -p "${OUT}"/{logs,polish,tmp}
cd "${OUT}"

# Use a large gpfs-backed TMPDIR (node /tmp is only ~8G).
export TMPDIR="${OUT}/tmp/${SLURM_JOB_ID:-manual}"
mkdir -p "${TMPDIR}"

# Clean any stale NextPolish rundir from a prior failed run.
rm -rf "${OUT}/polish/rundir"

for f in "${ASM}" "${SR1}" "${SR2}"; do
  [[ -f "${f}" ]] || { echo "ERROR: missing input ${f}" >&2; exit 1; }
done

# downsample
SUB1="${OUT}/polish/SRR11461820_sub_1.fastq.gz"
SUB2="${OUT}/polish/SRR11461820_sub_2.fastq.gz"
echo "[$(date)] rasusa: subsampling short reads to ${TARGET_COV}x (genome=${ASM})"
rasusa reads -c "${TARGET_COV}" -g "${ASM}.fai" -v \
  -o "${SUB1}" -o "${SUB2}" "${SR1}" "${SR2}" 2>&1 | tee "${OUT}/logs/rasusa.log" \
  || { [[ -f "${ASM}.fai" ]] || samtools faidx "${ASM}"; \
       rasusa reads -c "${TARGET_COV}" -g "${ASM}.fai" -v -o "${SUB1}" -o "${SUB2}" "${SR1}" "${SR2}" 2>&1 | tee "${OUT}/logs/rasusa.log"; }

# POLISH (NextPolish, short reads -> downsampled)
GAPFILLED="${ASM}"
FINAL="${GAPFILLED}"
PDIR="${OUT}/polish"
printf '%s\n%s\n' "${SUB1}" "${SUB2}" > "${PDIR}/sgs.fofn"
cat > "${PDIR}/run.cfg" <<EOF
[General]
job_type = local
job_prefix = nextpolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 1
multithread_jobs = 24
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

seqkit stats -T "${ASM}" "${GAPFILLED}" "${FINAL}" | tee "${OUT}/logs/final.stats.tsv"
echo "[$(date)] DONE -> ${FINAL}"
