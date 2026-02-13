#!/bin/bash
#SBATCH --job-name=hifiasm_hydmaj
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH -p icelake
#SBATCH --output=hifiasm_%j.out
#SBATCH --error=hifiasm_%j.err
#SBATCH --mail-user=13billy.trim13@gmail.com
#SBATCH --mail-type=ALL

module purge
module load Anaconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /hpcfs/users/a1864358/miniconda/envs/hifiasm

cd $SLURM_SUBMIT_DIR

HIFI="/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/hydmaj.fastq.gz"
HIC_R1="/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/350845_R1.fastq.gz"
HIC_R2="/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/350845_R2.fastq.gz"
OUTDIR="/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm"
mkdir -p $OUTDIR
cd $OUTDIR

ls -la hydmaj.fastq.gz 350845_R1.fastq.gz 350845_R2.fastq.gz
file hydmaj.fastq.gz

hifiasm \
    -o hmaj_26 \
    -t ${SLURM_CPUS_PER_TASK} \
    --primary \
    --dual-scaf \
    --hom-cov 30 \
    -u 1 \
    --h1 350845_R1.fastq.gz \
    --h2 350845_R2.fastq.gz \
    hydmaj.fastq.gz
