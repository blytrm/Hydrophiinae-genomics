#!/usr/bin/env bash
#SBATCH --job-name=05_eviann
#SBATCH -p highmem
#SBATCH -A honours
#SBATCH -N 1
#SBATCH --cpus-per-task=40
#SBATCH --time=10:00:00
#SBATCH --mem=260GB
#SBATCH -o /hpcfs/users/a1864358/sanders_lab/asm/files/annotation/pipeline/joblogs/%x_%j.out
#SBATCH -e /hpcfs/users/a1864358/sanders_lab/asm/files/annotation/pipeline/joblogs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com
# Step 05 - Run EviAnn.
# Inputs:
#   - genome FASTA            (${GENOME})
#   - sorted BAMs from step 4 (auto-collected into ${FOFN_DIR}/eviann.fofn)
#   - cleaned protein FASTA   (${OUT_PROTQC}/proteins.clean.faa  from step 3)
#   - optional UniProt        (${UNIPROT_REF})
# Submit: sbatch 05_eviann.sh

set -euo pipefail
module purge

SCRIPT_DIR="/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/pipeline"
source "${SCRIPT_DIR}/00_config.sh"

source "${CONDA_PROFILE}"
conda activate "${EVIANN_ENV}"

WORK="${OUT_EVIANN}/${GENOME_LABEL}"
mkdir -p "${WORK}"
cd "${WORK}"

BAM_FOFN="${FOFN_DIR}/eviann.fofn"
PROTEINS="${OUT_PROTQC}/proteins.clean.faa"

# Build BAM fofn from step 04 outputs
: > "${BAM_FOFN}"
shopt -s nullglob
for bam in "${OUT_ALIGN}/bam"/*.sorted.bam; do
    base="$(basename "${bam}")"
    case "${base}" in
        *isoseq*|*pacbio*) tag="bam_isoseq" ;;
        *)                 tag="bam" ;;
    esac
    printf '%s %s\n' "${bam}" "${tag}" >> "${BAM_FOFN}"
done
shopt -u nullglob

if [[ ! -s "${BAM_FOFN}" ]]; then
    echo "ERROR: no BAMs found in ${OUT_ALIGN}/bam (run step 04 first)" >&2
    exit 1
fi

if [[ ! -s "${PROTEINS}" ]]; then
    echo "ERROR: missing protein FASTA ${PROTEINS} (run step 03 first)" >&2
    exit 1
fi

echo "[$(date)] BAM fofn (${BAM_FOFN}):"
cat "${BAM_FOFN}"

# Run EviAnn
EVIANN_ARGS=(
    -t "${SLURM_CPUS_PER_TASK:-16}"
    -g "${GENOME}"
    -r "${BAM_FOFN}"
    -p "${PROTEINS}"
)
if [[ -n "${UNIPROT_REF}" && -s "${UNIPROT_REF}" ]]; then
    EVIANN_ARGS+=(-s "${UNIPROT_REF}")
fi

echo "[$(date)] eviann.sh ${EVIANN_ARGS[*]}"
eviann.sh "${EVIANN_ARGS[@]}"
echo "[$(date)] EviAnn finished. Outputs in ${WORK}"
