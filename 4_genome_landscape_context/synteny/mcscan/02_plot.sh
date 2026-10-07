#!/usr/bin/env bash
# Build karyotype inputs (seqids, layout, correspondence) then render the ribbon PNG.
set -uo pipefail
PROJ=/hpcfs/users/a1864358/sanders_lab/asm/files/synteny/mcscan-synteny
source /home/a1864358/jcvi-venv/bin/activate
cd "$PROJ"

# 1. ordering + layout + correspondence table (uses .anchors.simple blocks)
python scripts/build_plot_inputs.py

# 2. render karyotype ribbon plot
cd "$PROJ/plot"
export MPLBACKEND=Agg
python -m jcvi.graphics.karyotype \
  --chrstyle=roundrect --basepair \
  --figsize=16x11 --dpi=300 --format=png \
  --outfile="$PROJ/plot/hydrophis_outgroup_synteny.png" \
  seqids layout 2>&1 | tail -6

echo "=== output ==="
ls -la "$PROJ/plot/hydrophis_outgroup_synteny.png" 2>/dev/null
