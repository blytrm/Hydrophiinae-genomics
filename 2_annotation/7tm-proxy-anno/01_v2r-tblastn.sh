#!/usr/bin/env bash
#SBATCH --job-name=v2r-tblastn
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --mem=32GB
#SBATCH -o %x_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com

set -eu
PATH=/hpcfs/users/a1864358/miniconda/envs/general/bin:$PATH

SYNT=/hpcfs/users/a1864358/sanders_lab/asm/files/synteny
QUERY=/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/v2r/hma_7tm.fa
SPECIES=(hcurw hmaj horn hcy)

THREADS=4

for sp in "${SPECIES[@]}"; do
  for which in new old; do
    DB=$SYNT/$sp/${which}.fa
    OUT=$SYNT/$sp/v2r.${which}.tblastn.outfmt6
    [ -f ${DB}.nhr ] || makeblastdb -input_type fasta -in $DB -dbtype nucl -parse_seqids
    tblastn -query $QUERY -db $DB -out $OUT -outfmt 6 -num_threads $THREADS -evalue 1e-10 &
  done
done
wait

for sp in "${SPECIES[@]}"; do
  for which in new old; do
    OUT=$SYNT/$sp/v2r.${which}.tblastn.outfmt6
    BED=$SYNT/$sp/v2r.${which}.bed
    awk 'BEGIN{OFS="\t"} {
      s=$9; e=$10; if (s>e){t=s;s=e;e=t}
      print $2, s-1, e
    }' $OUT | sort -k1,1 -k2,2n > $BED
    echo "[$sp/$which] hits: $(wc -l < $OUT)  bed rows: $(wc -l < $BED)"
  done
done

echo "DONE"
