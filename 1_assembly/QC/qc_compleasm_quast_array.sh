#!/usr/bin/env bash
#SBATCH --job-name=quast_compleasm
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=40
#SBATCH --time=00:59:00
#SBATCH --mem=50GB
#SBATCH -o %x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com

set -euo pipefail
module purge

export PATH="/scratchdata1/users/a1864358/miniconda/bin:$PATH"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /hpcfs/users/a1864358/miniconda/envs/compleasm

WORKDIR="/hpcfs/users/a1864358/sanders_lab/asm/files/hcur-E"
cd "$WORKDIR"

asms=(
    'hcurE-denovo.fa'
    'hcure-ins.fa'
    'hcure-polish1.fa'
    'hcure-haphic.fa'
    'hcure-autohic.fa'
    'hcure-curated.fa'
    'hcure-finalpolish.fa'
)

for i in "${asms[@]}"; do 
    base="${i%.fa}"
    q_out="quast_${base}"
    b_out="busc_${base}"
    
    quast "$i" --threads 40 --eukaryote --plots-format png --no-icarus -o "$q_out"

    compleasm run -a "$i" -o "$b_out" -l squamata -t 40

    rm -rf mb_downloads 2>/dev/null || true
done
