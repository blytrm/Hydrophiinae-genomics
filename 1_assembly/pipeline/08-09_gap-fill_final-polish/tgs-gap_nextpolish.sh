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

# gapfill (TGSGapCloser) + polish (NextPolish) for assemblies


set -euo pipefail

ASM="${1:?asm.fa}"
RAW="${2:?reads}"
TYPE="${3:?ont|hifi|pb}"
OUT="${4:?out_dir}"
SR1="${5:-}"
SR2="${6:-}"

THREADS="${SLURM_CPUS_PER_TASK:-32}"
MIN_ONT_LEN=30000

case "${TYPE}" in
  ont)  TGSTYPE=ont; MINIMAP='-x map-ont' ;;
  hifi) TGSTYPE=pb;  MINIMAP='-x asm20' ;;
  pb)   TGSTYPE=pb;  MINIMAP='-x map-pb' ;;
  *) echo "TYPE must be ont|pb|hifi" >&2; exit 2 ;;
esac

mkdir -p "${OUT}"/{logs,gapfill,polish}
cd "${OUT}"

# PREP READS (length-filter ONT only)
if [[ "${TYPE}" == "ont" ]]; then
  READS="${OUT}/reads.${MIN_ONT_LEN}.fa"
  [[ -s "${READS}" ]] || seqkit seq -j "${THREADS}" -m "${MIN_ONT_LEN}" -o "${READS}" "${RAW}"
else
  READS="${RAW}"
fi
seqkit stats -T "${READS}" | tee "${OUT}/logs/reads.stats.tsv"

# GAPFILL (TGSGapCloser)
NGAPS=$(seqtk gap "${ASM}" | wc -l)
echo "[$(date)] ${NGAPS} gaps in ${ASM}"

GAP_PREFIX="${OUT}/gapfill/gapclosed"
GAPFILLED="${ASM}"
if (( NGAPS > 0 )); then
  tgsgapcloser \
    --scaff   "${ASM}" \
    --reads   "${READS}" \
    --output  "${GAP_PREFIX}" \
    --tgstype "${TGSTYPE}" \
    --minmap_arg "${MINIMAP}" \
    --thread  "${THREADS}" \
    --ne \
    >"${OUT}/logs/tgsgapcloser.out" \
    2>"${OUT}/logs/tgsgapcloser.err"
  for c in "${GAP_PREFIX}.scaff_seqs" "${GAP_PREFIX}.scaff_seq"; do
    [[ -s "$c" ]] && GAPFILLED="$c" && break
  done
fi
seqkit stats -T "${ASM}" "${GAPFILLED}" | tee "${OUT}/logs/gapfill.stats.tsv"

# POLISH (NextPolish, short reads) --> skip if no SR provided
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
