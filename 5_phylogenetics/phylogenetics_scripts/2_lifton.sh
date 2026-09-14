#!/bin/bash
#SBATCH -J lifton
#SBATCH -a 0-1    # 0=nscut 1=ptext
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o /hpcfs/users/a1864358/sanders_lab/resources/terstrl-elapids_liftover-anno/logs/lifton_%A_%a.out

# Lift 2018 RefSeq GFF -> 2025 assembly with LiftOn (Liftoff DNA + miniprot protein)

set -euo pipefail

RES=/hpcfs/users/a1864358/sanders_lab/resources
OUT=$RES/terstrl-elapids_liftover-anno
T=${SLURM_CPUS_PER_TASK:-16}

SP=(nscut ptext)
s=${SP[$SLURM_ARRAY_TASK_ID]}

case $s in
  nscut) REF=$RES/nscut/nscut.fa;  GFF=$RES/nscut/genomic.gff;  TGT=$RES/nscut/2025/nscut_2025.fa ;;
  ptext) REF=$RES/ptext/ptext.fa;  GFF=$RES/ptext/genomic.gff;  TGT=$RES/ptext/2025/ptext_2025.fa ;;
  *) echo "bad species $s"; exit 1 ;;
esac

o=$OUT/$s
mkdir -p "$o"
LIFTED=$o/${s}.lifton.gff3

source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh

# pre: source annotation stats 
set +u; conda activate general; set -u
agat_sp_statistics.pl --gff "$GFF" -o "$o/${s}_source.agat_stats.txt"

set +u; conda activate lifton; set -u
cd "$o"
lifton \
  -g "$GFF" \
  -o "$LIFTED" \
  -u "$o/${s}.unmapped.txt" \
  -a 0.50 -s 0.50 -flank 0.1 \
  -copies -sc 0.95 \
  -cds -polish \
  --validate-output \
  -t "$T" \
  "$TGT" "$REF"

set +u; conda activate general; set -u
agat_sp_statistics.pl --gff "$LIFTED" -o "$o/${s}_lifted.agat_stats.txt"

# flag genes with overlapping coords after lift (paralog collapse check)
conda activate general
awk '$3=="gene"{print $1"\t"$4-1"\t"$5"\t"$9}' "$LIFTED" | sort -k1,1 -k2,2n > "$o/${s}.genes.bed"
bedtools merge -i "$o/${s}.genes.bed" -c 1 -o count > "$o/${s}.genes.merged_count.bed"
awk -F'\t' '$4>1' "$o/${s}.genes.merged_count.bed" > "$o/${s}.overlapping_gene_loci.bed" || true

echo "$s done. lifted=$LIFTED  overlaps=$(wc -l < "$o/${s}.overlapping_gene_loci.bed") loci"
