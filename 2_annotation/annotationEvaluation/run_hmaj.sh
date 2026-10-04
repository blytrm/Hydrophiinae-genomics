#!/usr/bin/env bash
# Launch postEviann-eval on hmaj post-annotation outputs
set -euo pipefail

NF_DIR="/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/nf"
POST="/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/post-anno-res"
GENOME="/hpcfs/users/a1864358/sanders_lab/asm/files/final-asms/10-hmaj-final.renamed.fa"
OMARK_DB="/scratchdata1/users/a1864358/dbs/omamer/LUCA.h5"
ETE_DB="/scratchdata1/users/a1864358/dbs/etetoolkit/taxa.sqlite"
COMPLEASM_LIB="/scratchdata1/users/a1864358/dbs/compleasm/mb_downloads"
WORK="/scratchdata1/users/a1864358/sanders_lab/annotation/nf_work"
OUTDIR="/scratchdata1/users/a1864358/sanders_lab/annotation/nf_results"
LOGDIR="/scratchdata1/users/a1864358/sanders_lab/annotation/nf_logs"
CONDA_CACHE="/scratchdata1/users/a1864358/sanders_lab/annotation/conda"
NXF_HOME_DIR="/scratchdata1/users/a1864358/sanders_lab/annotation/nxf_home"

mkdir -p "$WORK" "$OUTDIR" "$LOGDIR" "$CONDA_CACHE" "$NXF_HOME_DIR"

if [[ ! -s "$OMARK_DB" ]]; then
  echo "ERROR: missing OMAmer DB: $OMARK_DB" >&2
  exit 1
fi
if [[ ! -s "$ETE_DB" ]]; then
  echo "ERROR: missing ete3 taxonomy DB: $ETE_DB" >&2
  exit 1
fi
if [[ ! -d "$COMPLEASM_LIB" ]]; then
  echo "ERROR: missing compleasm library: $COMPLEASM_LIB" >&2
  exit 1
fi

export NXF_HOME="$NXF_HOME_DIR"
export NXF_CONDA_CACHEDIR="$CONDA_CACHE"
export NXF_OPTS="-Xms1g -Xmx4g"
export NXF_ANSI_LOG=false

cd "$NF_DIR"

nextflow run main.nf \
  -profile slurm,conda \
  -work-dir "$WORK" \
  -resume \
  --proteins      "${POST}/hmaj.fa.functional_note.proteins.fasta" \
  --cds           "${POST}/hmaj.fa.functional_note.transcripts.fasta" \
  --gff3          "${POST}/hmaj.fa.functional_note.pseudo_label.gff" \
  --genome_fasta  "$GENOME" \
  --omark_db      "$OMARK_DB" \
  --ete_db        "$ETE_DB" \
  --compleasm_lib "$COMPLEASM_LIB" \
  --busco_lineage squamata \
  --outdir        "$OUTDIR" \
  -with-report    "${LOGDIR}/report.html" \
  -with-timeline  "${LOGDIR}/timeline.html" \
  -with-trace     "${LOGDIR}/trace.txt" \
  "$@"
