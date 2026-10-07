#!/usr/bin/env bash
# One MCscan pairwise comparison in an isolated subdir (safe for parallel runs).
# diamond builds pairwise protein hits (.last = BLAST m8); jcvi reuses them for
# synteny chaining -> anchors -> .simple blocks.  No `last` binary required.
# Usage: 00_align_pair.sh <A> <B> [threads]
set -uo pipefail
A=$1; B=$2; T=${3:-6}
PROJ=/hpcfs/users/a1864358/sanders_lab/asm/files/synteny/mcscan-synteny
DIAMOND=/hpcfs/users/a1864358/miniconda/envs/eggnog/bin/diamond
WD="$PROJ/work/pairs/${A}__${B}"
mkdir -p "$WD" "$PROJ/anchors"
source /home/a1864358/jcvi-venv/bin/activate

cd "$WD" || exit 1
cp -f "$PROJ/bed/$A.bed" "$A.bed"; cp -f "$PROJ/bed/$B.bed" "$B.bed"
ln -sf "$PROJ/pep/$A.pep" "$A.pep"; ln -sf "$PROJ/pep/$B.pep" "$B.pep"

# 1. diamond pairwise protein alignment -> A.B.last (BLAST tabular / m8)
if [[ ! -s "$A.$B.last" ]]; then
  "$DIAMOND" makedb --in "$B.pep" -d "$B" --quiet 2>/dev/null
  "$DIAMOND" blastp -q "$A.pep" -d "$B" -o "$A.$B.last" \
    --outfmt 6 --evalue 1e-5 --max-target-seqs 5 -p "$T" --quiet 2>/dev/null
fi
touch "$A.$B.last"   # newer than inputs so jcvi skips the last aligner

# 2. jcvi: filter + synteny scan + anchors (reuses the existing .last)
python -m jcvi.compara.catalog ortholog --dbtype prot --no_strip_names --notex "$A" "$B" 2>&1 | tail -3

# 3. collapse to collinear blocks for the ribbon plot
python -m jcvi.compara.synteny screen --minspan=30 --simple "$A.$B.anchors" "$A.$B.anchors.new" 2>&1 | tail -2

cp -f "$A.$B.anchors" "$A.$B.anchors.simple" "$A.$B.last" "$PROJ/anchors/" 2>/dev/null
echo "[done] $A $B hits=$(wc -l < "$A.$B.last" 2>/dev/null) anchors=$(wc -l < "$A.$B.anchors" 2>/dev/null) simple=$(wc -l < "$A.$B.anchors.simple" 2>/dev/null)"
