#!/bin/bash
# Run GenMap (v1.3.0) mappability on two horn genomes.

# (30,2)-mappability is the documented example; outputs text/wig/bedGraph.

set -uo pipefail

WORKDIR=/scratchdata1/users/a1864358/sanders_lab/asm/files/final-asms/tandem-repeat-anal
GENMAPDIR="$WORKDIR/genmap"
CONDA=/hpcfs/users/a1864358/miniconda/miniconda3
K=30
E=2
THREADS=8

source "$CONDA/etc/profile.d/conda.sh"
conda activate genmap

mkdir -p "$GENMAPDIR"
export TMPDIR="$GENMAPDIR/tmp"; mkdir -p "$TMPDIR"
genmap --version 2>/dev/null | head -1 || true

# name  fasta-path
run_one() {
  local name="$1" fasta="$2"
  local idx="$GENMAPDIR/${name}_index"
  local out="$GENMAPDIR/${name}_map_K${K}_E${E}"
  echo "=================================================================="
  echo "[$(date)] GENOME: $name  ($fasta)"

  if [ -d "$idx" ] && [ -f "$idx/index.txt.concat" ]; then
    echo "[$(date)] index exists, skipping build: $idx"
  else
    rm -rf "$idx"   # genmap requires the index dir to not pre-exist
    echo "[$(date)] Building index (-A divsufsort, ~20GB RAM) -> $idx"
    genmap index -F "$fasta" -I "$idx" -A divsufsort
  fi

  mkdir -p "$out"
  echo "[$(date)] Computing (${K},${E})-mappability -> $out"
  genmap map -K "$K" -E "$E" -I "$idx" -O "$out/${name}" -t -w -bg -T "$THREADS"
  echo "[$(date)] DONE $name. Outputs:"
  ls -la "$out"
}

run_one final  /scratchdata1/users/a1864358/sanders_lab/asm/files/final-asms/10-horn-final.renamed.fa
run_one og     /hpcfs/users/a1864358/sanders_lab/asm/files/horn/horn-og.fa

echo "[$(date)] All GenMap runs complete."
