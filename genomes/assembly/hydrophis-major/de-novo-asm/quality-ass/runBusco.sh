#!/bin/bash
#SBATCH --job-name=compleasm
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH -p icelake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com

module purge
module load Python
module load Anaconda3

source $(conda info --base)/etc/profile.d/conda.sh
conda activate /hpcfs/users/a1864358/miniconda/envs/compleasm

cd /hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm

ASM='hmaj_denovo.fa'
LINEAGE_DIR='/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/mb_downloads/'
OUTPUT='/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/compleasm'
OUTPUT_ASM='/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/compleasm/asm'

mkdir -p ${OUTPUT}

compleasm run -a ${ASM} -o ${OUTPUT_ASM} -l squamata -L ${LINEAGE_DIR} -t ${SLURM_CPUS_PER_TASK} 

echo "Compleasm finished."
