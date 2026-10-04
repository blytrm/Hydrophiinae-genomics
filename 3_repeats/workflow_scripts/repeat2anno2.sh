#!/usr/bin/env bash
#SBATCH --job-name=teforest_merge
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --time=2-00:00:00
#SBATCH --mem=200G
#SBATCH --chdir=/hpcfs/users/a1864358/sanders_lab/asm/files/repeat-anno/ra2
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=13billy.trim13@gmail.com


set -euo pipefail

BASE=/hpcfs/users/a1864358/sanders_lab/asm/files/repeat-anno/ra2
cd "$BASE"

THREADS=${SLURM_CPUS_PER_TASK:-72}
SAMPLE=SRR16961054


for d in */; do printf "%10s  %s\n" "$(find "$d" -xdev | wc -l)" "$d"; done | sort -rn


GENOME="$BASE/data/10-horn-final.renamed.fa"
R1="$BASE/data/${SAMPLE}_1.fastq.gz"
R2="$BASE/data/${SAMPLE}_2.fastq.gz"

CONSENSUS_TES="$BASE/results/hybridTEs.fasta"
REF_TE_BED="$BASE/results/reference_TEs_teforest.bed"            # 7-col reference TEs (family in col7)
EUCHROMATIN="$BASE/results/euchromatin_wholegenome.bed"
RM_OUT=/hpcfs/users/a1864358/sanders_lab/asm/files/repeat-anno/edta_10-horn/repeatmasker_out/10-horn-final.renamed.fa.out
TEFOREST_MODEL="$BASE/models/TEforest_Nonreference50X.pkl"
TEFOREST_MODEL_REF="$BASE/models/TEforest_reference50X.pkl"

TEFOREST_PY=/hpcfs/users/a1864358/sanders_lab/prog/TEforest/TEforest.py
TEFOREST_WORKFLOW=/hpcfs/users/a1864358/sanders_lab/prog/TEforest/workflow
RELOCATE2_PY=/hpcfs/users/a1864358/sanders_lab/prog/RelocaTE2/scripts/relocaTE2.py
MERGE_PY="$BASE/merge_te_annotations.py"

RELOCATE2_OUTDIR="$BASE/results/02_relocate2"
TEFOREST_OUTDIR="$BASE/results/03_teforest"
MERGE_OUTDIR="$BASE/results/04_merged"
mkdir -p logs "$RELOCATE2_OUTDIR" "$TEFOREST_OUTDIR" "$MERGE_OUTDIR"

TEFOREST_BED="$TEFOREST_OUTDIR/output/${SAMPLE}_TEforest_bps_nonredundant.bed"
RELOCATE2_GFF="$RELOCATE2_OUTDIR/repeat/results/ALL.all_nonref_insert.gff"

source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh

echo "=== [$(date)] Pipeline start | threads=$THREADS | sample=$SAMPLE ==="

# Stage 1: TEforest (short-read non-reference + reference TE genotyping)
if [[ -s "$TEFOREST_BED" ]]; then
    echo "=== [$(date)] TEforest: output exists, skipping ==="
else
    echo "=== [$(date)] TEforest: running ==="
    conda run --no-capture-output -n TEforest python "$TEFOREST_PY" \
        --workflow_dir "$TEFOREST_WORKFLOW" \
        --workdir "$TEFOREST_OUTDIR" \
        --threads "$THREADS" \
        --samples "$SAMPLE" \
        --fq_base_path "$BASE/data" \
        --consensusTEs "$CONSENSUS_TES" \
        --ref_genome "$GENOME" \
        --ref_te_locations "$REF_TE_BED" \
        --euchromatin "$EUCHROMATIN" \
        --model "$TEFOREST_MODEL" \
        --ref_model "$TEFOREST_MODEL_REF" \
        --cleanup_intermediates
    if [[ ! -s "$TEFOREST_BED" ]]; then
        echo "ERROR: TEforest finished but missing $TEFOREST_BED"
        exit 1
    fi
    echo "=== [$(date)] TEforest: done ==="

    # Scrub huge regenerable intermediates; keep output/ finals + config.yaml
    echo "=== [$(date)] TEforest: scrubbing intermediates ==="
    rm -rf \
        "$TEFOREST_OUTDIR/fastp" \
        "$TEFOREST_OUTDIR/aligned" \
        "$TEFOREST_OUTDIR/.snakemake" \
        "$TEFOREST_OUTDIR/benchmarks" \
        "$TEFOREST_OUTDIR/logs" \
        "$TEFOREST_OUTDIR/ref_genome"
    echo "=== [$(date)] TEforest size after scrub: $(du -sh "$TEFOREST_OUTDIR" 2>/dev/null | cut -f1) ==="
fi

# Stage 2: merge into a unified call set (EDTA reference loci + TEforest non-ref).
conda run --no-capture-output -n TEforest python "$MERGE_PY" \
    --teforest "$TEFOREST_BED" \
    --ref_bed "$REF_TE_BED" \
    --window 100 \
    --out_gff3 "$MERGE_OUTDIR/TE_annotation_unified.gff3" \
    --out_bed "$MERGE_OUTDIR/TE_annotation_unified.bed"

