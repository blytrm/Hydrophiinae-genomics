#!/usr/bin/env bash
#SBATCH --job-name=ins
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=40
#SBATCH --time=2-00:00:00
#SBATCH --mem=240GB
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com

WORKDIR="/hpcfs/users/a1864358/sanders_lab/asm/files/horn"
CONTIG="${WORKDIR}/horn-denovo.fa"
READS_RAW="${WORKDIR}/hor.fastq.gz"
READS_SUB="${WORKDIR}/hor_subsamp.fastq.gz"
INS_OUT="${WORKDIR}/insp_horn_out"
CORRECTED_OUT="${WORKDIR}/ins-corrected"

INSPECTOR_BIN='/hpcfs/users/a1864358/sanders_lab/Inspector'

set -euxo pipefail
module purge

export PATH="/hpcfs/users/a1864358/miniconda/envs/ins/bin:$PATH"
source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
conda activate /hpcfs/users/a1864358/miniconda/envs/ins

cd "${WORKDIR}"

mkdir -p "${INS_OUT}"

python "${INSPECTOR_BIN}/inspector.py" \
    -c "${CONTIG}" \
    -r "${READS_SUB}" \
    --datatype nanopore \
    -o "${INS_OUT}" \
    --noplot \
    -t "${SLURM_CPUS_PER_TASK}"

python "${INSPECTOR_BIN}/inspector-correct.py" \
    -i "${INS_OUT}/" \
    --datatype nano-raw \
    -o "${CORRECTED_OUT}/" \
    -t "${SLURM_CPUS_PER_TASK}"

ASM="${CORRECTED_OUT}/contig_corrected.fa"
CLR="${READS_RAW}"

if [[ ! -s "${ASM}" ]]; then
    echo "Expected corrected assembly missing or empty: ${ASM}" >&2
    exit 1
fi
