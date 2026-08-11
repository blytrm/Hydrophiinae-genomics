#!/usr/bin/env bash
#SBATCH --job-name=04_align_rna
#SBATCH -p icelake
#SBATCH -A honours
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=40GB
#SBATCH -o /hpcfs/users/a1864358/sanders_lab/asm/files/annotation/pipeline/joblogs/%x_%A_%a.out
#SBATCH -e /hpcfs/users/a1864358/sanders_lab/asm/files/annotation/pipeline/joblogs/%x_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=13billy.trim13@gmail.com


# reads ${FOFN_DIR}/rna-reads.tsv with columns:
#     <sample>\t<R1>\t<R2_or_->\t<type>
# where type is one of: fastq | isoseq


set -euo pipefail
module purge

SCRIPT_DIR="/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/pipeline"
source "${SCRIPT_DIR}/00_config.sh"

if [[ -n "${HISAT2_MODULE:-}" ]]; then
    module load "${HISAT2_MODULE}"
fi

source "${CONDA_PROFILE}"
conda activate "${ALIGN_ENV}"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

CLEANED_DIR="${OUT_RNAQC}/cleaned"
INDEX_DIR="$(dirname "${HISAT2_IDX}")"
SAMPLES_TSV="${FOFN_DIR}/rna-reads.tsv"
BAM_DIR="${OUT_ALIGN}/bam"
TMP_DIR="${OUT_ALIGN}/tmp"
mkdir -p "${INDEX_DIR}" "${BAM_DIR}" "${TMP_DIR}"

generate_samples_tsv() {
    echo "[$(date)] Generating ${SAMPLES_TSV} from ${CLEANED_DIR}"
    : > "${SAMPLES_TSV}"
    for r1 in "${CLEANED_DIR}"/*_R1.clean.fastq.gz; do
        [[ -f "${r1}" ]] || continue
        sample="$(basename "${r1}" _R1.clean.fastq.gz)"
        r2="${CLEANED_DIR}/${sample}_R2.clean.fastq.gz"
        if [[ -f "${r2}" ]]; then
            printf '%s\t%s\t%s\tfastq\n' "${sample}" "${r1}" "${r2}" >> "${SAMPLES_TSV}"
        fi
    done
    # IsoSeq long reads (optional): drop *.iso.fastq.gz files into ${RNA_DIR}/isoseq/
    # and uncomment below to include them. They go through minimap2 instead of HISAT2.
    # for lr in "${RNA_DIR}"/isoseq/*.fastq.gz; do
    #     [[ -f "${lr}" ]] || continue
    #     sample="$(basename "${lr}" .fastq.gz)"
    #     printf '%s\t%s\t-\tisoseq\n' "${sample}" "${lr}" >> "${SAMPLES_TSV}"
    # done
    echo "Wrote $(wc -l < "${SAMPLES_TSV}") samples to ${SAMPLES_TSV}"
}

prepare_inputs() {
    generate_samples_tsv
    if [[ ! -s "${SAMPLES_TSV}" ]]; then
        echo "ERROR: no cleaned read pairs found in ${CLEANED_DIR} (run step 02 first)" >&2
        exit 1
    fi

    if [[ ! -s "${HISAT2_IDX}.1.ht2" ]]; then
        echo "[$(date)] Building HISAT2 index for ${GENOME} -> ${HISAT2_IDX}"
        hisat2-build -p "${THREADS}" "${GENOME}" "${HISAT2_IDX}"
    fi
}

if [[ "${1:-}" == "--prepare" ]]; then
    prepare_inputs
    exit 0
fi

if [[ ! -s "${SAMPLES_TSV}" ]]; then
    echo "ERROR: missing ${SAMPLES_TSV}" >&2
    echo "Run first: bash ${BASH_SOURCE[0]} --prepare" >&2
    exit 1
fi

if [[ ! -s "${HISAT2_IDX}.1.ht2" ]]; then
    echo "ERROR: missing HISAT2 index ${HISAT2_IDX}.1.ht2" >&2
    echo "Run first: bash ${BASH_SOURCE[0]} --prepare" >&2
    exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: this script must be submitted as a SLURM array job." >&2
    echo "  sbatch --array=1-$(wc -l < "${SAMPLES_TSV}") 04_align_rna.sh" >&2
    exit 1
fi

LINE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES_TSV}")"
read -r SAMPLE R1 R2 TYPE <<< "${LINE}"

OUT_BAM="${BAM_DIR}/${SAMPLE}.sorted.bam"
echo "[$(date)] task=${SLURM_ARRAY_TASK_ID} sample=${SAMPLE} type=${TYPE}"

if [[ "${TYPE}" == "fastq" ]]; then
    hisat2 \
        -x "${HISAT2_IDX}" \
        -p "${THREADS}" \
        --min-intronlen "${MIN_INTRON}" \
        --max-intronlen "${MAX_INTRON}" \
        -1 "${R1}" \
        -2 "${R2}" |
    samtools sort -@ 4 -m 2G -T "${TMP_DIR}/${SAMPLE}" -O BAM -o "${OUT_BAM}"

elif [[ "${TYPE}" == "isoseq" ]]; then
    cat \
        <(cat "${GENOME}" <(echo "") | ufasta sizes -H | awk '{print "@SQ\tSN:"$1"\tLN:"$2}') \
        <(minimap2 -a -u f -x splice -t "${THREADS}" -G "${MAX_INTRON}" "${GENOME}" "${R1}" | grep -v '^@') |
    samtools sort -@ 4 -m 2G -T "${TMP_DIR}/${SAMPLE}" -O BAM -o "${OUT_BAM}"

else
    echo "ERROR: unknown TYPE='${TYPE}' for sample ${SAMPLE}" >&2
    exit 1
fi

samtools index "${OUT_BAM}"
echo "[$(date)] done: ${OUT_BAM}"
