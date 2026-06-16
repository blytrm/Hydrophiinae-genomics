#!/usr/bin/env bash
#SBATCH --job-name=morehic
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=42
#SBATCH --time=24:00:00
#SBATCH --mem=240GB
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com
 
# AutoHiC morehic pipeline: detect and correct assembly errors using Hi-C.

set -euo pipefail

ASM="/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/files/hcy.backup0.backup0/hcy-scaff1.fasta"
OUT="/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/files/hcy.backup0.backup0/resources-autohic/try2/1"
BAM="/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/files/hcy.backup0.backup0/hcy-1.bam"

TOOLS_DIR="/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/files/hcy.backup0.backup0/autohic-work"
JUICEBOX_SCRIPTS="${TOOLS_DIR}/juicebox_scripts/juicebox_scripts"
THREE_D_DNA="${TOOLS_DIR}/3d-dna"
AUTOHIC_DIR="${TOOLS_DIR}/AutoHiC"
ERROR_MODEL="${TOOLS_DIR}/models/error_model.pth"

MOREHIC_ENV="/hpcfs/users/a1864358/miniconda/envs/morehic"
AUTOHIC_ENV="/hpcfs/users/a1864358/miniconda/envs/autohic"

module purge
export PATH="/scratchdata1/users/a1864358/miniconda/bin:$PATH"
source "$(conda info --base)/etc/profile.d/conda.sh"

PREFIX=$(basename "${ASM}" .fa)
PREFIX=$(basename "${PREFIX}" .fasta)

mkdir -p "${OUT}"
cd "${OUT}"

# fasta -> agp -> .assembly
conda activate "${MOREHIC_ENV}"

if [[ -f "${OUT}/${PREFIX}.assembly" ]]; then
else
    python3 "${JUICEBOX_SCRIPTS}/makeAgpFromFasta.py" \
        "${ASM}" "${OUT}/${PREFIX}.agp"

    python3 "${JUICEBOX_SCRIPTS}/agp2assembly.py" \
        "${OUT}/${PREFIX}.agp" "${OUT}/${PREFIX}.assembly"
fi
echo ">>>   ${OUT}/${PREFIX}.assembly"

# BAM -> .hic
if [[ -f "${OUT}/${PREFIX}.hic" ]]; then
else
    if [[ -f "${OUT}/${PREFIX}.sorted.links.txt" ]]; then
    else
        matlock bam2 juicer "${BAM}" "${OUT}/${PREFIX}.links.txt"

        sort -k2,2 -k6,6 --parallel="${SLURM_CPUS_PER_TASK}" -S 50% -T "${OUT}" \
            "${OUT}/${PREFIX}.links.txt" > "${OUT}/${PREFIX}.sorted.links.txt"

        rm -f "${OUT}/${PREFIX}.links.txt"
    fi

    bash "${THREE_D_DNA}/visualize/run-assembly-visualizer.sh" \
        -p false \
        "${OUT}/${PREFIX}.assembly" \
        "${OUT}/${PREFIX}.sorted.links.txt"
fi
echo ">>>   ${OUT}/${PREFIX}.hic"

# AutoHiC : OneHiC correction
conda activate "${AUTOHIC_ENV}"

python3 "${AUTOHIC_DIR}/onehic.py" \
    -hic "${OUT}/${PREFIX}.hic" \
    -asy "${OUT}/${PREFIX}.assembly" \
    -autohic "${AUTOHIC_DIR}" \
    -p "${ERROR_MODEL}" \
    -out "${OUT}/" \
    -t "${SLURM_CPUS_PER_TASK}"

echo ">>>   ${OUT}/adjusted.assembly"

# .assembly -> fasta
conda activate "${MOREHIC_ENV}"

python3 "${JUICEBOX_SCRIPTS}/juicebox_assembly_converter.py" \
    -a "${OUT}/adjusted.assembly" \
    -f "${ASM}"


