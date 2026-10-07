#!/usr/bin/env bash
#SBATCH --job-name=nucmer-all
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=48
#SBATCH --time=24:00:00
#SBATCH --mem=128GB
#SBATCH -o %x_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com

set -eu
PATH=/hpcfs/users/a1864358/miniconda/envs/general/bin:$PATH

SYNT=/hpcfs/users/a1864358/sanders_lab/asm/files/synteny
SPECIES=(hcurw hmaj horn hcy)

HEADER="reference_start,reference_end,query_start,query_end,reference_alignment_length,query_alignment_length,percent_identity,reference_length,query_length,reference_chromosome,query_chromosome"

# Run all species sequentially (each nucmer uses all 48 threads — fast on icelake).
for sp in "${SPECIES[@]}"; do
  cd $SYNT/$sp
  echo "=== nucmer $sp ==="
  date
  # ref = old, query = new (matches plot-synteny-run.R orientation)
  nucmer --prefix algn old.fa new.fa -t ${SLURM_CPUS_PER_TASK}
  show-coords -lTH algn.delta > algn.coords
  (echo "${HEADER}" && awk '{$1=$1}1' OFS="," algn.coords) > algn.csv
  awk -F',' 'NR==1 || ($5>=500 && $6>=500)' algn.csv > filtered-500-algn.csv
  echo "[$sp] csv rows: $(wc -l < algn.csv)  filtered-500: $(wc -l < filtered-500-algn.csv)"
done

echo "DONE"
date
