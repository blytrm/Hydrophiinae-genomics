#!/usr/bin/env bash
# positive control on Hydrophis cyanocinctus.
#
#   bash processing/01_positive_control.sh [species]      (default hcy)
#
# Generates the two primary evidence tracks for the species, runs the curation
# pipeline on it, and scores the result against the independent -> Li et al. 2021 H. cyanocinctus benchmark.
#
# Why this species: nothing from H. cyanocinctus is in the V2R profile's training
# set, and it already has an answer derived from a *different* annotation
# (funannotate, 30,750 proteins) — 458 V2R-family / 426 confident / ~408 V2R.

set -euo pipefail
export LC_ALL=C

SP="${1:-hcy}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ROOT="/scratchdata1/users/a1864358/sanders_lab/asm/files/annotation/eviann_results"
ANNO="/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/eviann_results"
V2RHMM="/hpcfs/users/a1864358/sanders_lab/v2r_hmm"
GENOME="/hpcfs/users/a1864358/sanders_lab/asm/files/final-asms/10-${SP}-final.renamed.fa"
QUERY="${ROOT}/hmaj/query_nr.faa"
BASE="${ROOT}/${SP}"
REPORT="${HERE}/processing/01_control_report.md"

mkdir -p "${BASE}"

# miniprot -> primary track : spliced alignment of V2R query proteins
if [[ ! -s "${BASE}/miniprot.v2r.gff" ]]; then
    echo "[01] miniprot on ${SP} (~20 min)"
    miniprot -t "${THREADS}" --gff -G 20k "${GENOME}" "${QUERY}" \
        > "${BASE}/miniprot.v2r.gff" 2> "${BASE}/miniprot.err"
fi
echo "[01] miniprot alignments: $(grep -cP '\tmRNA\t' "${BASE}/miniprot.v2r.gff")"

# V2R profile vs the EviAnn proteome 
if [[ ! -s "${BASE}/${SP}.prot.domtbl" ]]; then
    echo "[01] hmmsearch on the ${SP} proteome"
    hmmsearch --cpu "${THREADS}" -E 1e-3 --domtblout "${BASE}/${SP}.prot.domtbl" \
        "${V2RHMM}/hmm/v2r_full.hmm" \
        "${ANNO}/${SP}/${SP}.fa.functional_note.proteins.fasta" > /dev/null
fi
echo "[01] proteome pHMM hits: $(grep -vc '^#' "${BASE}/${SP}.prot.domtbl")"

# run the curation pipeline 

bash "${HERE}/run_curation.sh" "${SP}"

# score against the benchmark 
OUT="${HERE}/results/${SP}"
read -r FAM CONF ARCH GENES INTACT < <(awk -F'\t' '
    NR==1 {for (i=1; i<=NF; i++) h[$i]=i; next}
    $h["classC_family"]=="V2R_OlfC" {fam++}
    $h["is_v2r"]=="yes" {
        conf++; g[$h["gene_id"]]
        if ($h["eggnog_arch_consistent"]=="yes") arch++
        if ($h["coverage_class"]=="intact") intact++
    }
    END {print fam, conf, arch, length(g), intact}' "${OUT}/v2r_master.tsv")

{
    printf '# 01 — positive control: %s\n\n' "${SP}"
    printf 'Benchmark: `v2r_hmm/V2R_hcyanocinctus_benchmark_report.md`, derived from the\n'
    printf 'funannotate proteome (30,750 proteins) with the same V2R profile.\n\n'
    printf '| Metric | Benchmark | This pipeline | Difference |\n|---|---|---|---|\n'
    printf '| V2R-family proteins | 458 | %s | %+.1f%% |\n' "$FAM" \
        "$(echo "scale=4; 100*($FAM-458)/458" | bc)"
    printf '| confident V2R transcripts | 426 | %s | %+.1f%% |\n' "$CONF" \
        "$(echo "scale=4; 100*($CONF-426)/426" | bc)"
    printf '| after removing non-V2R Class C | ~408 | %s | %+.1f%% |\n' "$ARCH" \
        "$(echo "scale=4; 100*($ARCH-408)/408" | bc)"
    printf '| **V2R genes** | **~408** | **%s** | **%+.1f%%** |\n' "$GENES" \
        "$(echo "scale=4; 100*($GENES-408)/408" | bc)"
    printf '| intact models only | — | %s | — |\n\n' "$INTACT"
    printf 'The benchmark counted **proteins from a one-protein-per-gene annotation**, so\n'
    printf 'the gene row is the like-for-like comparison; the transcript rows are inflated\n'
    printf 'by EviAnn isoforms.\n'
} > "${REPORT}"

echo "[01] wrote ${REPORT}"
cat "${REPORT}"
