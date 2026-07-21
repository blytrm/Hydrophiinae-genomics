#!/bin/bash
# Run RepeatOBserverV1 on the horn genome.

set -euo pipefail

WORKDIR=/scratchdata1/users/a1864358/sanders_lab/asm/files/final-asms/tandem-repeat-anal
GENOME=/scratchdata1/users/a1864358/sanders_lab/asm/files/final-asms/10-horn-final.renamed.fa
SPECIES=Horn           
HAP=H0
CPU=4                   
MEM=40000               
CONDA=/hpcfs/users/a1864358/miniconda/miniconda3

source "$CONDA/etc/profile.d/conda.sh"
conda activate repeatobserver

cd "$WORKDIR"
ln -sf "$GENOME" ./horn_genome.fasta

chmod +x Setup_Run_Repeats.sh
echo "[$(date)] Launching RepeatOBserver: -i $SPECIES -h $HAP -c $CPU -m $MEM -g FALSE"
bash Setup_Run_Repeats.sh -i "$SPECIES" -f horn_genome.fasta -h "$HAP" -c "$CPU" -m "$MEM" -g FALSE

echo "[$(date)] RepeatOBserver complete."
echo "Summary output: $WORKDIR/output_chromosomes/${SPECIES}_${HAP}-AT/Summary_output/"
