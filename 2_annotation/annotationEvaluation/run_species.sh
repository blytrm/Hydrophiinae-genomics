#!/usr/bin/env bash
# Launch postEviann-eval on one species' EviAnn outputs.
set -euo pipefail

LABEL="${1:?usage: run_species.sh <label> (hcure|hcurw|hcy|hmaj|horn)}"

NF_DIR="/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/nf"
EVI="/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/eviann_results/${LABEL}"
GENOME="/hpcfs/users/a1864358/sanders_lab/asm/files/final-asms/10-${LABEL}-final.renamed.fa"
OMARK_DB="/scratchdata1/users/a1864358/dbs/omamer/LUCA.h5"
ETE_DB="/scratchdata1/users/a1864358/dbs/etetoolkit/taxa.sqlite"
COMPLEASM_LIB="/scratchdata1/users/a1864358/dbs/compleasm/mb_downloads"

BASE="/scratchdata1/users/a1864358/sanders_lab/annotation"
WORK="${BASE}/nf_work_${LABEL}"                                        # transient, on scratch
# Publish final results to persistent hpcfs storage (one dir per species),
# matching the hmaj precedent and avoiding scratch's inode limits / purge.
OUTDIR="/hpcfs/users/a1864358/sanders_lab/annotation/nf_results/${LABEL}"
LOGDIR="${BASE}/nf_logs"
CONDA_CACHE="${BASE}/conda"          # shared env cache (already built)
NXF_HOME_DIR="${BASE}/nxf_home"

PROT="${EVI}/${LABEL}.fa.functional_note.proteins.fasta"
CDS="${EVI}/${LABEL}.fa.functional_note.transcripts.fasta"
GFF="${EVI}/${LABEL}.fa.functional_note.pseudo_label.gff"

mkdir -p "$WORK" "$OUTDIR" "$LOGDIR" "$CONDA_CACHE" "$NXF_HOME_DIR"

for f in "$PROT" "$CDS" "$GFF" "$GENOME" "$OMARK_DB" "$ETE_DB"; do
  [[ -s "$f" ]] || { echo "ERROR[$LABEL]: missing/empty input: $f" >&2; exit 1; }
done
[[ -d "$COMPLEASM_LIB" ]] || { echo "ERROR[$LABEL]: missing compleasm library: $COMPLEASM_LIB" >&2; exit 1; }

export NXF_HOME="$NXF_HOME_DIR"
export NXF_CONDA_CACHEDIR="$CONDA_CACHE"
export NXF_OPTS="-Xms1g -Xmx4g"
export NXF_ANSI_LOG=false

cd "$NF_DIR"

nextflow run main.nf \
  -profile slurm,conda \
  -work-dir "$WORK" \
  -resume \
  --proteins      "$PROT" \
  --cds           "$CDS" \
  --gff3          "$GFF" \
  --genome_fasta  "$GENOME" \
  --omark_db      "$OMARK_DB" \
  --ete_db        "$ETE_DB" \
  --compleasm_lib "$COMPLEASM_LIB" \
  --busco_lineage squamata \
  --outdir        "$OUTDIR" \
  -with-report    "${LOGDIR}/report_${LABEL}.html" \
  -with-timeline  "${LOGDIR}/timeline_${LABEL}.html" \
  -with-trace     "${LOGDIR}/trace_${LABEL}.txt"
