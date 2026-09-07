#!/usr/bin/env bash
#SBATCH --job-name=flye-hor
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH -c 40
#SBATCH --time=72:00:00
#SBATCH --mem=600GB
#SBATCH -o %x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

# flye run on Hydrophis ornatus 
DIR="/home/${USER}/al-biohub/billy_trim/honours"
FQ="${DIR}/data/sra-data/hydrophis_ornatus-SAMN23222008/nanopore_wgs/hor.fastq.gz"
OUT="${DIR}/01_genome_assembly/hydrophis_ornatus/contig/assembly"

flye \
  --nano-hq "${FQ}" \
  --out-dir "${OUT}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --genome-size 2g \
  --iterations 2
