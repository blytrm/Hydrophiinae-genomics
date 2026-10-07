#!/usr/bin/env bash
#SBATCH --job-name=v2r-tblastn-full
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH -o %x_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com

# tblastn v2r query vs FULL new + original assemblies (all seqs, incl scaffolds).
set -eu
PATH=/hpcfs/users/a1864358/miniconda/envs/general/bin:$PATH

FILES=/hpcfs/users/a1864358/sanders_lab/asm/files
SYNT=$FILES/synteny
QUERY=$FILES/annotation/v2r/hma_7tm.fa
OUTDIR=$SYNT/v2r-full
mkdir -p $OUTDIR

declare -A NEW OLD
NEW[hcurw]=$FILES/final-asms/10-hcurw-final.fa
NEW[hmaj]=$FILES/final-asms/10-hmaj-final.fa
NEW[horn]=$FILES/final-asms/10-horn-final.fa
NEW[hcy]=$FILES/final-asms/10-hcy-final.fa
OLD[hcurw]=$FILES/hcur-w/hydrophis_curtus-west.fa
OLD[hmaj]=$FILES/hmaj/hydrophis_major.fa
OLD[horn]=$FILES/horn/horn-og.fa
OLD[hcy]=$FILES/hcy/files/hcy.backup0.backup0/00-hcy-li-original.fa

# 8 jobs (4 sp x 2 asms) x 4 threads = 32 cores
THREADS=4

run_one() {
  local sp=$1 which=$2 fa=$3
  local db=$OUTDIR/${sp}.${which}
  local out=$OUTDIR/${sp}.${which}.outfmt6
  [ -f $fa.fai ] || samtools faidx $fa
  cp -f $fa.fai $OUTDIR/${sp}.${which}.fai
  if [ ! -f ${db}.nin ] && [ ! -f ${db}.nal ]; then
    makeblastdb -in $fa -dbtype nucl -out $db >$OUTDIR/${sp}.${which}.mkdb.log 2>&1
  fi
  tblastn -query $QUERY -db $db -out $out -outfmt 6 -num_threads $THREADS -evalue 1e-10
  echo "[$sp/$which] hits: $(wc -l < $out)"
}

for sp in hcurw hmaj horn hcy; do
  run_one $sp new "${NEW[$sp]}" &
  run_one $sp old "${OLD[$sp]}" &
done
wait
echo "ALL TBLASTN DONE"
