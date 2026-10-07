#!/usr/bin/env bash
# Apply the same V2R filter as annotation/v2r/filter-tblastn-7tm-hits.r:
#   qstart==1 AND qend==length AND length>=200 AND evalue<=1e-10 AND pident>=60
# Then bedtools merge (gap=0) for loci, count per chromosome.

set -eu
PATH=/hpcfs/users/a1864358/miniconda/envs/general/bin:$PATH
SYNT=/hpcfs/users/a1864358/sanders_lab/asm/files/synteny

filter_to_bed() {
  # outfmt6 cols: 1=qseqid 2=sseqid 3=pident 4=length 7=qstart 8=qend 9=sstart 10=send 11=evalue
  awk 'BEGIN{OFS="\t"}
       $7==1 && $8==$4 && $4>=200 && $11<=1e-10 && $3>=60 {
         s=$9; e=$10; if (s>e){t=s;s=e;e=t}
         print $2, s-1, e
       }' "$1" | sort -k1,1 -k2,2n
}

echo -e "species\tasm\tchr\tloci"
for sp in hcurw hmaj horn hcy hcure; do
  for w in new old; do
    OUT=$SYNT/$sp/v2r.${w}.tblastn.outfmt6
    [ -f $OUT ] || continue
    F=$SYNT/$sp/v2r.${w}.filtered.bed
    filter_to_bed $OUT > $F
    M=$SYNT/$sp/v2r.${w}.filtered.merged.bed
    bedtools merge -i $F > $M
    awk -v sp="$sp" -v w="$w" '{c[$1]++} END{for (k in c) printf "%s\t%s\t%s\t%d\n", sp, w, k, c[k]}' $M
  done
done | sort -k1,1 -k2,2 -k3,3
