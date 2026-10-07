#!/usr/bin/env bash
# Run all 11 adjacent pairwise MCscan comparisons in parallel.
set -uo pipefail
PROJ=/hpcfs/users/a1864358/sanders_lab/asm/files/synteny/mcscan-synteny
cd "$PROJ"
mkdir -p logs/pairs

PAIRS=(
  "horn hmaj" "hmaj hcure" "hcure hcurw" "hcurw hcy" "hcy nscut"
  "nscut ptext" "ptext nhel" "nhel tele" "tele vberus" "vberus caspera" "caspera vkomodo"
)
for pr in "${PAIRS[@]}"; do
  set -- $pr; A=$1; B=$2
  bash scripts/00_align_pair.sh "$A" "$B" 6 > "logs/pairs/${A}__${B}.log" 2>&1 &
done
wait
echo "=== ALL PAIRS DONE ==="
for pr in "${PAIRS[@]}"; do
  set -- $pr; A=$1; B=$2
  tail -1 "logs/pairs/${A}__${B}.log"
done