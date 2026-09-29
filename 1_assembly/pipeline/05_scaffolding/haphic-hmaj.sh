#!/usr/bin/env bash
#SBATCH --job-name=haphic-hmaj
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=48
#SBATCH --time=24:00:00
#SBATCH --mem=240GB
#SBATCH -o %x_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com

# haphic for hydrophis major

module purge

export PATH="/scratchdata1/users/a1864358/miniconda/bin:$PATH"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /hpcfs/users/a1864358/miniconda/envs/haphic

cd /hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hmaj
mkdir -p haphic

# bwa index hmaj_filt_denovo.fa

# bwa mem -5SP -t ${SLURM_CPUS_PER_TASK} hmaj_filt_denovo.fa ./350845_R1.fastq.gz ./350845_R2.fastq.gz | samblaster | samtools view - -@ 14 -S -h -b -F 3340 -o HiC-hmaj.bam

# /hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/HapHiC/utils/filter_bam HiC-hmaj.bam 1 --nm 3 --threads 20 | samtools view - -b -@ 14 -o HiC-hmaj.filtered.bam

FA="/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/hmaj_filt_denovo.fa"
BAM="/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/HiC-hmaj.filtered.bam"
OUT="/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/haphic"
SW="/hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hifiasm/HapHiC"


# purge dups
${SW}/utils/purge_dups -2 -T ${SLURM_CPUS_PER_TASK} -c ${FA}.fai -b ${BAM} -o ${BAM}.pbed

${SW}/haphic pipeline \
    ${FA} ${BAM} 18 \
    --RE "GATC,GANTC,CTNAG,TTAA" \
    --outdir ${OUT}/n18_micro \
    --Nx 100 \
    --min_inflation 1.1 \
    --max_inflation 7.0 \
    --inflation_step 0.1 \
    --correct_nrounds 2 \
    --threads ${SLURM_CPUS_PER_TASK} \
    --processes ${SLURM_CPUS_PER_TASK} \
    --verbose


rg "recommend_inflation|grouped together|maximum number" \
  ${OUT}/n18_micro/01.cluster/HapHiC_cluster.log > ${OUT}/n18_micro/01.cluster/HapHiC_cluster_summary.txt