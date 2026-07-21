#!/bin/bash
# Run Tandem Repeats Finder (TRF v4.09.1) on the horn genome.
# https://github.com/Benson-Genomics-Lab/TRF
#
# Recommended genome-wide parameters (from TRF README):
#   trf <file> Match Mismatch Delta PM PI Minscore MaxPeriod
#   2 7 7 80 10 50 2000   (PM=80 PI=10 are the recommended best-performance probs)
#   -l 10  : longest expected tandem array = 10 Mbp (large centromeric satellites)
#   -h     : suppress bulky per-sequence HTML
#   -ngs   : compact, parseable single-stream output (implies -d/-h)
#
# TRF is single-threaded -> split the genome per-sequence 
# in parallel across chromosomes, then concatenate the .dat output.
set -euo pipefail

WORKDIR=/scratchdata1/users/a1864358/sanders_lab/asm/files/final-asms/tandem-repeat-anal
GENOME=/scratchdata1/users/a1864358/sanders_lab/asm/files/final-asms/10-horn-final.renamed.fa
JOBS=12                 
PARAMS="2 7 7 80 10 50 2000"
OPTS="-l 10 -h -ngs"

CONDA=/hpcfs/users/a1864358/miniconda/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
conda activate general          
export PATH=/gpfs/apps/icl/software/TRF/4.09.1-GCCcore-11.2.0/bin:$PATH

cd "$WORKDIR"
OUT=trf_horn
mkdir -p "$OUT/split" "$OUT/dat"

echo "[$(date)] Splitting genome per sequence ..."
seqkit split -i -O "$OUT/split" "$GENOME" >/dev/null 2>&1 || \
  seqkit split2 -i -O "$OUT/split" "$GENOME"

echo "[$(date)] Running TRF ($JOBS parallel) with params: $PARAMS $OPTS"
ls "$OUT"/split/*.fa* | xargs -I{} -P "$JOBS" bash -c '
  f="$1"; out="$2"
  base=$(basename "$f"); base="${base%.*}"
  trf "$f" '"$PARAMS"' '"$OPTS"' > "$out/dat/${base}.dat" 2>/dev/null || true
  echo "  done: ${base}"
' _ {} "$OUT"

echo "[$(date)] Concatenating per-sequence .dat -> $OUT/horn.trf.ngs.dat"
cat "$OUT"/dat/*.dat > "$OUT/horn.trf.ngs.dat"

echo "[$(date)] TRF complete."
echo "Repeats found: $(grep -vc '^@' "$OUT/horn.trf.ngs.dat" 2>/dev/null || echo NA) data lines"
echo "Output: $OUT/horn.trf.ngs.dat"
