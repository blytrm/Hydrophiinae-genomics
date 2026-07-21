#!/bin/bash
# pHMM(miniprot gene-level) ∩ tBLASTn locus overlap, name-fixed



set -uo pipefail
source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
conda activate general
WD=/hpcfs/users/a1864358/sanders_lab/v2r_hmm
A=/hpcfs/users/a1864358/sanders_lab/asm/files/final-asms
T=/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/tblastn/ll
TMP=$WD/validation/multi_scan/overlap_beds
mkdir -p $TMP
printf "%-6s %-6s %-7s %-6s %-7s | %-7s %-7s %-7s\n" sp all ov_tb tb ov_ph intact ov_tb ov_ph
for s in hcurw hcy hmaj horn; do   # hcure: no tBLASTn bed yet
  awk -F'\t' 'NR==FNR{if(FNR>1)m[$1]=$2; next}{if($1 in m)$1=m[$1]; print $1"\t"$2"\t"$3}' OFS='\t' \
    $A/10-$s-final.rename_map.tsv $T/${s}_v2r_loci.bed | sort -k1,1 -k2,2n > $TMP/${s}_tb_ren.bed
  awk -F, 'NR>1{print $1"\t"$2"\t"$3}'                  $WD/validation/intact/${s}_v2r_genes.csv | sort -k1,1 -k2,2n > $TMP/${s}_all.bed
  awk -F, 'NR>1 && $9=="intact"{print $1"\t"$2"\t"$3}'  $WD/validation/intact/${s}_v2r_genes.csv | sort -k1,1 -k2,2n > $TMP/${s}_intact.bed
  nall=$(wc -l < $TMP/${s}_all.bed); nint=$(wc -l < $TMP/${s}_intact.bed); ntb=$(wc -l < $TMP/${s}_tb_ren.bed)
  a_ot=$(bedtools intersect -u -a $TMP/${s}_all.bed    -b $TMP/${s}_tb_ren.bed | wc -l)
  a_op=$(bedtools intersect -u -a $TMP/${s}_tb_ren.bed -b $TMP/${s}_all.bed    | wc -l)
  i_ot=$(bedtools intersect -u -a $TMP/${s}_intact.bed -b $TMP/${s}_tb_ren.bed | wc -l)
  i_op=$(bedtools intersect -u -a $TMP/${s}_tb_ren.bed -b $TMP/${s}_intact.bed | wc -l)
  printf "%-6s %-6s %-7s %-6s %-7s | %-7s %-7s %-7s\n" $s $nall $a_ot $ntb $a_op $nint $i_ot $i_op
done
