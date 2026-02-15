#!/usr/bin/env bash
#SBATCH --job-name=hifi_filt
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32     
#SBATCH --mem=64G              
#SBATCH --time=04:00:00        
#SBATCH --output=hifi_filt_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=13billy.trim13@gmail.com


module purge
module load Anaconda3
module load BLAST+
module load Python

source $(conda info --base)/etc/profile.d/conda.sh
conda activate /hpcfs/users/a1864358/miniconda/envs/general

export PATH=$PATH:/hpcfs/users/a1864358/miniconda/envs/general/bin/hifiadapterfilt.sh
export PATH=$PATH:/hpcfs/users/a1864358/miniconda/envs/general/bin/DB

cd /hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/hifilt

bash hifiadapterfilt.sh -t "${SLURM_CPUS_PER_TASK}"
