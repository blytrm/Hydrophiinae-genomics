#!/usr/bin/env bash
# generate the two primary evidence tracks 

set -euo pipefail

SP="$1"
THREADS="${2:-16}"
ROOT=/scratchdata1/users/a1864358/sanders_lab/asm/files/annotation/eviann_results
ANNO=/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/eviann_results
GENOME=/hpcfs/users/a1864358/sanders_lab/asm/files/final-asms/10-${SP}-final.renamed.fa
QUERY=${ROOT}/hmaj/query_nr.faa
V2R=/hpcfs/users/a1864358/sanders_lab/v2r_hmm/hmm/v2r_full.hmm
OUT=${ROOT}/${SP}

mkdir -p "${OUT}"

# primary track 1 — spliced alignment of 1,516 vomeronasal receptor proteins.
# -G 20k because V2R introns run to ~20 kb; the default would split real genes.
if [[ ! -s "${OUT}/miniprot.v2r.gff" ]]; then
    miniprot -t "${THREADS}" --gff -G 20k "${GENOME}" "${QUERY}" \
        > "${OUT}/miniprot.v2r.gff" 2> "${OUT}/miniprot.err"
fi

# primary track 2 — the V2R profile against this species' EviAnn proteome
if [[ ! -s "${OUT}/${SP}.prot.domtbl" ]]; then
    hmmsearch --cpu "${THREADS}" -E 1e-3 --domtblout "${OUT}/${SP}.prot.domtbl" \
        "${V2R}" "${ANNO}/${SP}/${SP}.fa.functional_note.proteins.fasta" > /dev/null
fi

printf '%s ready: %s miniprot alignments, %s proteome hits\n' "${SP}" \
    "$(awk -F'\t' '$3=="mRNA"' "${OUT}/miniprot.v2r.gff" | wc -l)" \
    "$(grep -vc '^#' "${OUT}/${SP}.prot.domtbl")"
