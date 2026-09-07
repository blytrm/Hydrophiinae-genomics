#!/usr/bin/env bash
# Quantify filtered V2R TBLASTN hits per sequence (chromosome / scaffold / contig),
# using the same filters as filter-tblastn-7tm-hits.r:
#   full query coverage (qstart==1, qend==length), length>=200, evalue<=1e-10, pident>=60
# Counts are merged non-overlapping intervals per seqid (bedtools merge), matching the old script logic.
#
# Inputs: old.7tm.tblastn.outfmt6 (00-hcy-li-original / GenBank+JAAZTL), new.7tm.tblastn.outfmt6 (10-hcy-final ch*/contig*)
# Reference naming: 00-hcy-li-original.fa.fai, 10-hcy-final.fa.fai
#
# Usage: from this directory, or set V2R_DIR:
#   V2R_DIR=/path/to/v2r ./tblastn-v2r-quant.sh
#
set -euo pipefail

V2R_DIR="${V2R_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
BEDTOOLS="${BEDTOOLS:-bedtools}"
TOP_N="${TOP_N:-25}"

OLD_OUTFMT6="${V2R_DIR}/old.7tm.tblastn.outfmt6"
NEW_OUTFMT6="${V2R_DIR}/new.7tm.tblastn.outfmt6"

for f in "$OLD_OUTFMT6" "$NEW_OUTFMT6"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing $f" >&2
    exit 1
  fi
done
command -v "$BEDTOOLS" >/dev/null || { echo "ERROR: bedtools not found (set BEDTOOLS=...)" >&2; exit 1; }

# outfmt6 columns: qseqid sseqid pident length ... qstart qend sstart send evalue bitscore (1-12)
# Emit BED (0-based start): sseqid, start, end
outfmt6_to_sorted_bed() {
  awk -F'\t' '
    $7 == 1 && $8 == $4 && $4 >= 200 && $11 <= 1e-10 && $3 >= 60 {
      s = ($9 < $10) ? $9 : $10
      e = ($9 < $10) ? $10 : $9
      print $2 "\t" (s - 1) "\t" e
    }
  ' "$1" | sort -k1,1 -k2,2n -k3,3n
}

count_merged_per_seq() {
  # stdin: sorted BED -> stdout: seqid, seq_class, merged_interval_count (all seqs)
  "$BEDTOOLS" merge -i - | awk -F'\t' '
    function classify(s) {
      if (s ~ /^ch[0-9]+$/) return "chromosome"
      if (s ~ /^contig[0-9]+$/) return "contig"
      if (s ~ /^CM[0-9]+\.[0-9]+$/) return "chromosome"
      if (s ~ /^JAAZTL/) return "scaffold"
      return "other"
    }
    { c[$1]++ }
    END {
      for (s in c) print s "\t" classify(s) "\t" c[s]
    }
  ' | sort -t $'\t' -k3,3nr -k1,1
}

whole_genome_merged_count() {
  outfmt6_to_sorted_bed "$1" | "$BEDTOOLS" merge -i - | wc -l | tr -d ' '
}

write_per_assembly() {
  local label="$1"
  local outfmt6="$2"
  local bed_out="${V2R_DIR}/${label}_v2r_filtered.bed"
  local per_seq="${V2R_DIR}/quant_${label}_per_seq.tsv"
  local top25="${V2R_DIR}/quant_${label}_top${TOP_N}.tsv"

  outfmt6_to_sorted_bed "$outfmt6" >"$bed_out"
  count_merged_per_seq <"$bed_out" >"$per_seq"
  head -n "$TOP_N" "$per_seq" >"$top25"

  echo "${label}_whole_genome_merged_loci	$(whole_genome_merged_count "$outfmt6")"
  echo "${label}_seqs_with_hits	$(awk -F'\t' '$3>0{c++} END{print c+0}' "$per_seq")"
}

echo "# V2R quantification (merged intervals per seq); TOP_N=${TOP_N}" >"${V2R_DIR}/quant_compare_summary.tsv"
write_per_assembly "old" "$OLD_OUTFMT6" >>"${V2R_DIR}/quant_compare_summary.tsv"
write_per_assembly "new" "$NEW_OUTFMT6" >>"${V2R_DIR}/quant_compare_summary.tsv"

# Side-by-side top lists (same row index = rank by count within each assembly; names are not homologous)
TOP_OLD="${V2R_DIR}/quant_old_top${TOP_N}.tsv"
TOP_NEW="${V2R_DIR}/quant_new_top${TOP_N}.tsv"
SIDE="${V2R_DIR}/quant_top${TOP_N}_old_vs_new.tsv"

{
  echo -e "rank\told_seqid\told_class\told_merged_loci\tnew_seqid\tnew_class\tnew_merged_loci"
  paste -d '\t' \
    <(awk -F'\t' -v OFS='\t' '{print NR,$1,$2,$3}' "$TOP_OLD") \
    <(awk -F'\t' -v OFS='\t' '{print $1,$2,$3}' "$TOP_NEW") \
  | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}'
} >"$SIDE"

echo "Wrote:"
echo "  ${V2R_DIR}/old_v2r_filtered.bed"
echo "  ${V2R_DIR}/new_v2r_filtered.bed"
echo "  ${V2R_DIR}/quant_old_per_seq.tsv  (all sequences, sorted by count)"
echo "  ${V2R_DIR}/quant_new_per_seq.tsv"
echo "  ${V2R_DIR}/quant_old_top${TOP_N}.tsv"
echo "  ${V2R_DIR}/quant_new_top${TOP_N}.tsv"
echo "  ${V2R_DIR}/quant_top${TOP_N}_old_vs_new.tsv  (ranks side-by-side; not 1:1 chromosome mapping)"
echo "  ${V2R_DIR}/quant_compare_summary.tsv"
cat "${V2R_DIR}/quant_compare_summary.tsv"
