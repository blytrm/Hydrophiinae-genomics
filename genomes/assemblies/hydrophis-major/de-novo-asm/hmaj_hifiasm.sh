#!/bin/bash
#SBATCH --job-name=hifiasm_hydmaj
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=750G
#SBATCH --time=3-00:00:00
#SBATCH -p highmem
#SBATCH --output=hifiasm_%j.out
#SBATCH --error=hifiasm_%j.err
#SBATCH --mail-user=13billy.trim13@gmail.com
#SBATCH --mail-type=ALL

module purge
cd $SLURM_SUBMIT_DIR

HIFI="/uofaresstor/sanders_lab/projects/billy/genome-assembly/seq/public/hydrophis_major/hifi/hydmaj.fastq.gz"
HIC_R1="/uofaresstor/sanders_lab/projects/billy/genome-assembly/seq/public/hydrophis_major/hic-trim/350845_R1.fastq.gz"
HIC_R2="/uofaresstor/sanders_lab/projects/billy/genome-assembly/seq/public/hydrophis_major/hic-trim/350845_R2.fastq.gz"
OUTDIR="/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm"
mkdir -p $OUTDIR
cd $OUTDIR

rm -f hydmaj.fastq.gz 350845_R1.fastq.gz 350845_R2.fastq.gz
ln -s $HIFI hydmaj.fastq.gz
ln -s $HIC_R1 350845_R1.fastq.gz
ln -s $HIC_R2 350845_R2.fastq.gz

hifiasm \\
    -o hmaj_26 \\
    -t ${SLURM_CPUS_PER_TASK} \\
    --primary \\
    --dual_scaf \\
    --hom-cov 30 \\
    -u 1 \\
    --h1 350845_R1.fastq.gz --h2 350845_R2.fastq.gz \\
    hydmaj.fastq.gz
