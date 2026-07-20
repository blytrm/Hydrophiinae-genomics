#!/usr/bin/env bash
#SBATCH --job-name=v2r-tblastn
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
# Run tblastn of V2R 7tm protein query against all final assemblies, then
# filter + merge HSPs into loci and count per sequence + genome-wide.
# Filters mirror filter-tblastn-7tm-hits.r / tblastn-v2r-quant.sh:
#   qstart==1, qend==length (full query cov), length>=200, evalue<=1e-10, pident>=60
# Loci = bedtools merge of passing HSPs per seqid (0 gap).
set -euo pipefail

ENVBIN=/hpcfs/users/a1864358/miniconda/envs/general/bin
export PATH="$ENVBIN:$PATH"

ASM_DIR=/scratchdata1/users/a1864358/sanders_lab/asm/files/final-asms
QUERY=/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/v2r/hma_7tm.fa
OUT=/scratchdata1/users/a1864358/sanders_lab/asm/files/annotation/v2r/final
DB=$OUT/db
THREADS=${SLURM_CPUS_PER_TASK:-24}
mkdir -p "$DB"

declare -A ASMS=(
  [hcurw]=$ASM_DIR/10-hcurw-final.fa
  [hcy]=$ASM_DIR/10-hcy-final.fa
  [hmaj]=$ASM_DIR/10-hmaj-final.fa
  [horn]=$ASM_DIR/10-horn-final.fa
)

filter_count() {
  # $1 outfmt6 -> writes loci bed, per-seq counts; echoes genome total
  local o6="$1" label="$2"
  local bed="$OUT/${label}_v2r_loci.bed"
  local perseq="$OUT/${label}_v2r_per_seq.tsv"
  awk -F'\t' '$7==1 && $8==$4 && $4>=200 && $11<=1e-10 && $3>=60 {
      s=($9<$10)?$9:$10; e=($9<$10)?$10:$9; print $2"\t"(s-1)"\t"e
    }' "$o6" | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > "$bed"
  awk -F'\t' '{c[$1]++} END{for(s in c) print s"\t"c[s]}' "$bed" \
    | sort -t$'\t' -k2,2nr -k1,1 > "$perseq"
  wc -l < "$bed" | tr -d ' '
}

SUMMARY=$OUT/v2r_loci_summary.tsv
echo -e "assembly\ttotal_v2r_loci\tseqs_with_loci" > "$SUMMARY"

for label in hcurw hcy hmaj horn; do
  ref="${ASMS[$label]}"
  echo "[$label] makeblastdb"
  if [[ ! -f "$DB/$label.nin" && ! -f "$DB/$label.nsq" ]]; then
    makeblastdb -in "$ref" -dbtype nucl -out "$DB/$label" >/dev/null
  fi
  echo "[$label] tblastn"
  o6="$OUT/${label}.7tm.tblastn.outfmt6"
  tblastn -query "$QUERY" -db "$DB/$label" -out "$o6" -outfmt 6 -num_threads "$THREADS"
  total=$(filter_count "$o6" "$label")
  nseq=$(wc -l < "$OUT/${label}_v2r_per_seq.tsv" | tr -d ' ')
  echo -e "${label}\t${total}\t${nseq}" >> "$SUMMARY"
  echo "[$label] done: $total loci across $nseq seqs"
done

echo "=== SUMMARY ==="
cat "$SUMMARY"
echo "DONE_TBLASTN_ALL"
