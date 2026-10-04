#!/usr/bin/env bash
# Run postEviann-eval on every EviAnn species sequentially.
# Default set is all five snakes (hmaj first). Each species gets its own work
# dir (scratch) and results dir (hpcfs) so runs never collide.
set -uo pipefail

# allow either "run_all.sh" (defaults) or "run_all.sh hcy horn"
if [[ $# -gt 0 ]]; then SPECIES=("$@"); else SPECIES=(hmaj hcure hcurw hcy horn); fi

NF_DIR="/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/nf"
LOGDIR="/scratchdata1/users/a1864358/sanders_lab/annotation/nf_logs"
mkdir -p "$LOGDIR"

cd "$NF_DIR"
echo "=== run_all started $(date) ; species: ${SPECIES[*]} ==="

for sp in "${SPECIES[@]}"; do
  LOG="${LOGDIR}/run_${sp}_$(date +%Y%m%d_%H%M%S).log"
  echo "$LOG" > "${LOGDIR}/LATEST_${sp}"
  echo "--- [$sp] launching $(date) -> $LOG ---"
  if bash run_species.sh "$sp" > "$LOG" 2>&1; then
    echo "--- [$sp] SUCCESS $(date) ---"
  else
    echo "--- [$sp] FAILED (exit $?) $(date); see $LOG ---"
  fi
done

echo "=== run_all finished $(date) ==="
