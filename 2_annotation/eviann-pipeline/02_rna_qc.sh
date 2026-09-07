#!/usr/bin/env bash
#SBATCH --job-name=02_rna_qc
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=40
#SBATCH --time=40:00:00
#SBATCH --mem=150GB
#SBATCH -o /hpcfs/users/a1864358/sanders_lab/asm/files/annotation/pipeline/joblogs/%x_%j.out
#SBATCH -e /hpcfs/users/a1864358/sanders_lab/asm/files/annotation/pipeline/joblogs/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com


# RNA-seq QC for EviAnn input.
# FastQC -> bbduk (adapters + rRNA trim) -> Kraken2 (host filter) -> fastp -> FastQC -> MultiQC.
# Cleaned => ${OUT_RNAQC}/cleaned/<sample>_R{1,2}.clean.fastq.gz

set -euo pipefail
module purge

SCRIPT_DIR="/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/pipeline"
source "${SCRIPT_DIR}/00_config.sh"

source "${CONDA_PROFILE}"
conda activate "${RNAQC_ENV}"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

CLEANED_DIR="${OUT_RNAQC}/cleaned"
OUT_LOGS="${OUT_RNAQC}/logs"
OUT_FASTQC_INIT="${OUT_RNAQC}/fastqc-init"
OUT_FASTQC_POST="${OUT_RNAQC}/fastqc-post"
OUT_BBDUK="${OUT_RNAQC}/bbduk"
OUT_KRAKEN2="${OUT_RNAQC}/kraken2"
OUT_FASTP="${OUT_RNAQC}/fastp"
OUT_RRNA="${OUT_RNAQC}/rrna"
OUT_FINAL="${OUT_RNAQC}/final"
OUT_MULTIQC="${OUT_RNAQC}/multiqc"

mkdir -p "${OUT_RNAQC}" "${CLEANED_DIR}" "${OUT_LOGS}" \
    "${OUT_FASTQC_INIT}" "${OUT_FASTQC_POST}" \
    "${OUT_BBDUK}" "${OUT_KRAKEN2}" "${OUT_FASTP}" \
    "${OUT_RRNA}" "${OUT_FINAL}" "${OUT_MULTIQC}"

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing required command: $1" >&2; exit 127; }
}

require_file() {
    [[ -f "$1" ]] || { echo "ERROR: missing required file: $1" >&2; exit 2; }
}

require_dir() {
    [[ -d "$1" ]] || { echo "ERROR: missing required directory: $1" >&2; exit 2; }
}

preflight() {
    require_dir "${RNA_DIR}"
    require_file "${ADAPTER_REF}"
    require_file "${RRNA_REF}"
    require_dir "${KRAKEN2_DB}"
    require_cmd fastqc
    require_cmd bbduk.sh
    require_cmd kraken2
    require_cmd fastp
    require_cmd multiqc
    require_cmd realpath
}

# Sample discovery (paired-end auto-detect)
find_pair() {
    local r1="$1"
    local candidate
    for candidate in \
        "${r1/_R1_/_R2_}" \
        "${r1/_R1./_R2.}" \
        "${r1/_1_/_2_}" \
        "${r1/_1./_2.}"; do
        if [[ "${candidate}" != "${r1}" ]] && [[ -f "${candidate}" ]]; then
            printf '%s' "${candidate}"
            return
        fi
    done
    printf '%s' "${r1/_R1_/_R2_}"
}

sample_id() {
    basename "${1}" | sed -E 's/[._]R?1([._][0-9]+)?\.(fastq|fq)(\.gz)?$//'
}

declare -a SAMPLES=()
declare -A SAMPLE_R1=()
declare -A SAMPLE_R2=()
declare -A SEEN_SAMPLE=()

populate_samples_array() {
    while IFS= read -r -d '' file; do
        local r1="${file}"
        local r2; r2="$(find_pair "${r1}")"
        if [[ -f "${r2}" ]]; then
            local sample; sample="$(sample_id "${r1}")"
            if [[ -n "${SEEN_SAMPLE[${sample}]:-}" ]]; then
                local base="${sample}"
                local n=2
                while [[ -n "${SEEN_SAMPLE[${base}__${n}]:-}" ]]; do
                    n=$((n + 1))
                done
                sample="${base}__${n}"
            fi
            SEEN_SAMPLE["${sample}"]=1
            SAMPLES+=("${sample}")
            SAMPLE_R1["${sample}"]="${r1}"
            SAMPLE_R2["${sample}"]="${r2}"
            echo "paired sample: ${sample} (${r1} + ${r2})"
        else
            echo "WARN: unpaired R1 (no matching R2): ${r1}"
        fi
    done < <(
        find "${RNA_DIR}" -type f \( \
            -name '*_R1.fastq' -o -name '*_R1.fastq.gz' -o \
            -name '*_R1_*.fastq' -o -name '*_R1_*.fastq.gz' -o \
            -name '*_R1.fq' -o -name '*_R1.fq.gz' -o \
            -name '*_R1_*.fq' -o -name '*_R1_*.fq.gz' -o \
            -name '*_1.fastq' -o -name '*_1.fastq.gz' -o \
            -name '*_1_*.fastq' -o -name '*_1_*.fastq.gz' -o \
            -name '*_1.fq' -o -name '*_1.fq.gz' -o \
            -name '*_1_*.fq' -o -name '*_1_*.fq.gz' \
        \) -print0 | sort -z
    )
}

# Per-sample processing
process() {
    local sample="$1"
    local raw_r1="$2"
    local raw_r2="$3"

    local fastp_r1="${OUT_FASTP}/${sample}.fastp_R1.fastq.gz"
    local fastp_r2="${OUT_FASTP}/${sample}.fastp_R2.fastq.gz"
    local bbduk_r1="${OUT_BBDUK}/${sample}.bbduk_R1.fastq.gz"
    local bbduk_r2="${OUT_BBDUK}/${sample}.bbduk_R2.fastq.gz"
    local kraken_pattern="${OUT_KRAKEN2}/${sample}.unclassified#.fastq"
    local kraken_r1="${OUT_KRAKEN2}/${sample}.unclassified_1.fastq"
    local kraken_r2="${OUT_KRAKEN2}/${sample}.unclassified_2.fastq"

    # skip if already cleaned
    if [[ -f "${CLEANED_DIR}/${sample}_R1.clean.fastq.gz" && -f "${CLEANED_DIR}/${sample}_R2.clean.fastq.gz" ]]; then
        echo "${sample} already done, skipping"
        return 0
    fi

    fastqc \
        --threads "${THREADS}" \
        --outdir "${OUT_FASTQC_INIT}" \
        "${raw_r1}" "${raw_r2}" \
        >>"${OUT_LOGS}/${sample}.initfastqc.log" 2>&1

    bbduk.sh \
        in1="${raw_r1}" \
        in2="${raw_r2}" \
        out1="${bbduk_r1}" \
        out2="${bbduk_r2}" \
        ref="${ADAPTER_REF},${RRNA_REF}" \
        ktrim=r \
        k=23 \
        mink=11 \
        hdist=1 \
        tpe \
        tbo \
        qtrim=rl \
        trimq=8 \
        minlength="${MINLEN}" \
        threads="${THREADS}" \
        stats="${OUT_BBDUK}/${sample}.bbduk.stats.txt" \
        >>"${OUT_LOGS}/${sample}.bbduk.log" 2>&1

    kraken2 \
        --db "${KRAKEN2_DB}" \
        --threads "${THREADS}" \
        --gzip-compressed \
        --paired \
        --use-names \
        --report "${OUT_KRAKEN2}/${sample}.kraken2.report" \
        --output /dev/null \
        --unclassified-out "${kraken_pattern}" \
        "${bbduk_r1}" "${bbduk_r2}" \
        >>"${OUT_LOGS}/${sample}.kraken2.log" 2>&1

    fastp \
        --in1 "${kraken_r1}" \
        --in2 "${kraken_r2}" \
        --out1 "${fastp_r1}" \
        --out2 "${fastp_r2}" \
        --detect_adapter_for_pe \
        --thread "${THREADS}" \
        --compression 6 \
        --length_required "${MINLEN}" \
        --average_qual "${FASTP_AVG_QUAL}" \
        --n_base_limit "${FASTP_N_BASE_LIMIT}" \
        --qualified_quality_phred "${FASTP_QUALIFIED_QUAL}" \
        --html "${OUT_FASTP}/${sample}.fastp.html" \
        --json "${OUT_FASTP}/${sample}.fastp.json" \
        --verbose \
        >>"${OUT_LOGS}/${sample}.fastp.log" 2>&1

    fastqc \
        --threads "${THREADS}" \
        --outdir "${OUT_FASTQC_POST}" \
        "${fastp_r1}" "${fastp_r2}" \
        >>"${OUT_LOGS}/${sample}.postfastqc.log" 2>&1

    ln -sf "$(realpath "${fastp_r1}")" "${CLEANED_DIR}/${sample}_R1.clean.fastq.gz"
    ln -sf "$(realpath "${fastp_r2}")" "${CLEANED_DIR}/${sample}_R2.clean.fastq.gz"

    # remove large intermediate files to conserve space
    rm -f "${bbduk_r1}" "${bbduk_r2}" "${kraken_r1}" "${kraken_r2}"

    echo "${sample} completed"
}

main() {
    preflight
    populate_samples_array
    if [[ ${#SAMPLES[@]} -eq 0 ]]; then
        echo "ERROR: no paired samples found in ${RNA_DIR}" >&2
        exit 1
    fi

    for sample in "${SAMPLES[@]}"; do
        process "${sample}" "${SAMPLE_R1[${sample}]}" "${SAMPLE_R2[${sample}]}"
    done

    multiqc \
        --force \
        --title "RNA-seq QC (prior to EviAnn annotation)" \
        --filename "rna-seq-qc-multiqc.html" \
        -o "${OUT_MULTIQC}" \
        "${OUT_FASTQC_INIT}" "${OUT_FASTQC_POST}" "${OUT_BBDUK}" "${OUT_KRAKEN2}" "${OUT_FASTP}" "${OUT_RRNA}" "${OUT_FINAL}"
}

main "$@"
