#!/usr/bin/env bash
# Functional annotation of the curated V2R proteins with eggNOG-mapper.
#
#   bash eggnog/run_emapper.sh [results_dir] [threads]
#
# Writes results/eggnog/v2r.emapper.annotations, which make_tables.py folds into
# v2r_master.tsv as the GO / KEGG / COG / PFAM columns on the next run.
set -euo pipefail

OUT="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/results}"
THREADS="${2:-16}"
SP="${3:-hmaj}"
DATA=/scratchdata1/users/a1864358/sanders_lab/prog/eggnog_db

source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
conda activate eggnog
export EGGNOG_DATA_DIR="$DATA"

for f in eggnog.db eggnog_proteins.dmnd; do
    [[ -s "$DATA/$f" ]] || { echo "ERROR: $DATA/$f missing - run fetch_db.sh" >&2; exit 1; }
done

EGG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/out/${SP}"
mkdir -p "${EGG}"

# Only the retained V2R proteins: the point is to characterise the repertoire,
# not to re-annotate every candidate the panel rejected.
awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) h[$i]=i; next}
    $h["is_v2r"]=="yes" {print $h["transcript_id"]}' "${OUT}/v2r_master.tsv" \
  | sort -u > "${EGG}/v2r_ids.txt"
seqkit grep -f "${EGG}/v2r_ids.txt" "${OUT}/candidates.faa" \
  > "${EGG}/v2r_proteins.faa"
echo "annotating $(grep -c '^>' "${EGG}/v2r_proteins.faa") V2R proteins"

# diamond mode: fast and sufficient for a few hundred proteins. --tax_scope auto
# lets eggNOG pick the narrowest orthologous group that still resolves.
emapper.py \
    -i "${EGG}/v2r_proteins.faa" --itype proteins \
    -m diamond --sensmode more-sensitive \
    --cpu "${THREADS}" --override \
    --output v2r --output_dir "${EGG}" \
    > "${EGG}/emapper.log" 2>&1

echo "annotated: $(grep -vc '^#' "${EGG}/v2r.emapper.annotations") rows"
echo "-> ${EGG}/v2r.emapper.annotations"
