#!/usr/bin/env bash
#SBATCH --job-name=rna-pre-anno
#SBATCH -p icelake
#SBATCH -A honours
#SBATCH -N 1
#SBATCH --cpus-per-task=24
#SBATCH --time=08:00:00
#SBATCH --mem=120GB
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com

# RNA-seq pre-processing for EviAnn genome annotation input

set -euo pipefail
module purge

source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
conda activate /hpcfs/users/a1864358/miniconda/envs/rnaqc


RNA_DIR="/scratchdata1/users/a1864358/sanders_lab/asm/files/annotation/pre-anno/rna_input_test"
OUT="/scratchdata1/users/a1864358/sanders_lab/asm/files/annotation/pre-anno/rna-qc"
CLEANED_DIR="${OUT}/cleaned"
THREADS=${SLURM_CPUS_PER_TASK}
OUT_LOGS="${OUT}/logs"

# databases
KRAKEN2_DB="/scratchdata1/users/a1864358/dbs/kraken2/k2_test_db"
ADAPTER_REF="/hpcfs/users/a1864358/miniconda/envs/general/share/bbmap/resources/adapters.fa"
RRNA_REF="/scratchdata1/users/a1864358/dbs/rrna_refs/rrna_combined.fa"

# parameters
MINLEN=50
FASTP_AVG_QUAL=15
FASTP_QUALIFIED_QUAL=15
FASTP_N_BASE_LIMIT=5

# outputs
OUT_FASTQC_INIT="${OUT}/fastqc-init"
OUT_FASTQC_POST="${OUT}/fastqc-post"
OUT_BBDUK="${OUT}/bbduk"
OUT_KRAKEN2="${OUT}/kraken2"
OUT_FASTP="${OUT}/fastp"
OUT_RRNA="${OUT}/rrna"
OUT_FINAL="${OUT}/final"
OUT_MULTIQC="${OUT}/multiqc"

mkdir -p "${OUT}" "${CLEANED_DIR}" "${OUT_LOGS}"

# find pairs
find_pair() {
    local r1="$1"
    local r2="${r1/_R1./_R2.}"
    [[ "${r2}" == "${r1}" ]] && r2="${r1/_1./_2.}"
    printf '%s' "${r2}"
}

declare -a SAMPLES=()
declare -A SAMPLE_R1=()
declare -A SAMPLE_R2=()

sample_id() {
    basename "${1}" | sed -E 's/[._](R)?1\.(fastq|fq)(\.gz)?$//'
}

populate_samples_array() {
    for file in "${RNA_DIR}"/*; do
        [[ -f "${file}" ]] || continue
        [[ "$file" =~ (_R2\.|_2\.) ]] && continue

        if [[ "$file" =~ (_R1\.|_1\.) ]]; then
            r1="${file}"
            r2="$(find_pair "${r1}")"
            if [[ -f "${r2}" ]]; then
                sample=$(sample_id "${r1}")
                SAMPLES+=("${sample}")
                SAMPLE_R1["${sample}"]="${r1}"
                SAMPLE_R2["${sample}"]="${r2}"
                echo "paired file found: ${r1} and ${r2}"
            else
                echo "unpaired r1: ${r1}"
            fi
        else 
            # single-end ?
            echo "single-end file: ${file}"
        fi 
    done
}

main() {
    mkdir -p \
        "${OUT_FASTQC_INIT}" \
        "${OUT_FASTQC_POST}" \
        "${OUT_BBDUK}" \
        "${OUT_KRAKEN2}" \
        "${OUT_FASTP}" \
        "${OUT_RRNA}" \
        "${OUT_FINAL}" \
        "${OUT_MULTIQC}" 
    
    for sample in "${SAMPLES[@]}"; do
        raw_r1="${SAMPLE_R1[${sample}]}"
        raw_r2="${SAMPLE_R2[${sample}]}"
        process "${sample}" "${raw_r1}" "${raw_r2}"
    done

    multiqc \
        --force \
        --title "RNA-seq QC (prior to EviAnn annotation)" \
        --filename "rna-seq-qc-multiqc.html" \
        -o "${OUT_MULTIQC}" \
        "${OUT_FASTQC_INIT}" "${OUT_FASTQC_POST}" "${OUT_BBDUK}" "${OUT_KRAKEN2}" "${OUT_FASTP}" "${OUT_RRNA}" "${OUT_FINAL}"
}


process() {
    local sample="$1"
    local raw_r1="$2"
    local raw_r2="$3"

    local fastp_r1="${OUT_FASTP}/${sample}.fastp_R1.fastq.gz"
    local fastp_r2="${OUT_FASTP}/${sample}.fastp_R2.fastq.gz"
    local bbduk_r1="${OUT_BBDUK}/${sample}.bbduk_R1.fastq.gz"
    local bbduk_r2="${OUT_BBDUK}/${sample}.bbduk_R2.fastq.gz"
    local kraken_pattern="${OUT_KRAKEN2}/${sample}.unclassified_#.fastq"
    local kraken_r1="${OUT_KRAKEN2}/${sample}.unclassified_1.fastq"
    local kraken_r2="${OUT_KRAKEN2}/${sample}.unclassified_2.fastq"

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
        --output "${OUT_KRAKEN2}/${sample}.kraken2.output" \
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

        echo "${sample} completed"
}

populate_samples_array
main "$@"