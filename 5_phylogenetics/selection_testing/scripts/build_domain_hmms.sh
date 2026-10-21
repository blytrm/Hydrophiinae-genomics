#!/usr/bin/env bash
# download Pfam domain HMMs and run hmmsearch on backbone.faa
# Usage: bash build_domain_hmms.sh [V2R_proteins.faa] [out_dir]
set -euo pipefail

FAA=${1:-/scratchdata1/users/a1864358/sanders_lab/phylogenomics/v2r_tree/workspace/data/familyC/backbone.faa}
H=${2:-/scratchdata1/users/a1864358/sanders_lab/phylogenomics/v2r_tree/workspace/data/familyC/domain_hmm}
S=/scratchdata1/users/a1864358/sanders_lab/phylogenomics/v2r_tree/workspace/scripts
mkdir -p "$H"


for acc in PF01094 PF07562 PF00003; do
    echo "downloading $acc..."
    curl -sL "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/${acc}?annotation=hmm" | gunzip > "$H/$acc.hmm"
done

# vft=PF01094 (ANF_receptor, VFT ectodomain)
# ncd3g=PF07562 (NCD3G, 9-cysteine-rich)
# tm7=PF00003 (7tm_3, 7-transmembrane core)
for f in vft:PF01094 ncd3g:PF07562 tm7:PF00003; do
    name=${f%%:*}; acc=${f##*:}
    echo "hmmsearch $name ($acc)..."
    hmmsearch --domtblout "$H/$name.domtbl" --noali --cpu 8 "$H/$acc.hmm" "$FAA" > "$H/$name.out"
done

# split into per-sequence domain boundaries
python3 -I "$S/split_domains.py" \
    "$H/vft.domtbl" "$H/ncd3g.domtbl" "$H/tm7.domtbl" \
    "$FAA" "$H/split"

echo "done — split.tsv at $H/split.tsv"
