#!/bin/bash
#SBATCH --job-name=lai_10horn
#SBATCH --partition=dm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=80G
#SBATCH --time=2-00:00:00
#SBATCH --output=lai_%j.log
#SBATCH --error=lai_%j.err

set -euo pipefail

source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
conda activate EDTA

PROJ=/scratchdata1/users/a1864358/sanders_lab/repeatdensity-atv2rs
RUN=$PROJ/edta_run
cd "$RUN"

G=10-horn-final.renamed.fa.mod
GENOME=$PROJ/10-horn-final.renamed.fa                 
PASS=$RUN/$G.EDTA.raw/LTR/$G.pass.list                
RMOUT=$RUN/$G.EDTA.anno/$G.out                        

LAI -genome "$GENOME" -intact "$PASS" -all "$RMOUT" -t 36

echo "DONE -> result: $RUN/$G.out.LAI"
