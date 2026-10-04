#!/usr/bin/env bash
#SBATCH --job-name=tetrimmer
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=3-00:00:00
#SBATCH --mem=180G
#SBATCH --chdir=/hpcfs/users/a1864358/sanders_lab/asm/files/repeat-anno/ra2
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail

BASE=/hpcfs/users/a1864358/sanders_lab/asm/files/repeat-anno/ra2
OUTDIR="$BASE/results/01_tetrimmer"
ENV=/hpcfs/users/a1864358/miniconda/envs/tetrimmer

mkdir -p logs "$OUTDIR"

export PERL5LIB="$ENV/lib/perl5/site_perl/5.22.0/x86_64-linux-thread-multi:$ENV/lib/perl5/site_perl/5.22.0:$ENV/share/pfam_scan-1.6-5${PERL5LIB:+:$PERL5LIB}"

source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
conda run --no-capture-output -n tetrimmer TEtrimmer \
    --input_file /hpcfs/users/a1864358/sanders_lab/asm/files/repeat-anno/edta_10-horn/10-horn-final.renamed.fa.mod.EDTA.TElib.fa \
    --genome_file "$BASE/data/10-horn-final.renamed.fa" \
    --output_dir "$OUTDIR" \
    --pfam_dir /hpcfs/users/a1864358/sanders_lab/prog/pfam_db \
    --num_threads "${SLURM_CPUS_PER_TASK:-48}" \
    --classify_all \
    --dedup
