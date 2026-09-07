#!/bin/bash
#SBATCH --job-name=rm_10horn
#SBATCH --partition=dm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=38
#SBATCH --mem=80G
#SBATCH --time=2-00:00:00
#SBATCH --output=rm_%j.log
#SBATCH --error=rm_%j.err

set -euo pipefail

source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
conda activate EDTA

PROJ=/scratchdata1/users/a1864358/sanders_lab/repeatdensity-atv2rs
GENOME=$PROJ/10-horn-final.renamed.fa                                          
LIB=$PROJ/edta_run/10-horn-final.renamed.fa.mod.EDTA.TElib.fa                  
OUT=$PROJ/repeatmasker_out

mkdir -p "$OUT"

RepeatMasker \
  -engine rmblast \
  -pa 9 \
  -lib "$LIB" \
  -xsmall \
  -gff \
  -a \
  -no_is \
  -gccalc \
  -dir "$OUT" \
  "$GENOME"

cd "$OUT"
ALIGN=$(basename "$GENOME").align

calcDivergenceFromAlign.pl -s "$(basename "$GENOME").divsum" "$ALIGN"

GENOME_SIZE=$(grep -v ">" "$GENOME" | tr -d '\n' | wc -c)

createRepeatLandscape.pl -div "$(basename "$GENOME").divsum" -g "$GENOME_SIZE" > "$(basename "$GENOME").landscape.html"

# extract BED / TE class from RepeatMasker .out
OUT_BASE=$(basename "$GENOME")
RM_OUT="${OUT_BASE}.out"
TMP_CLASS=.tmp_repeat_classes.bed

# BED6 (col5=score=0) + col7=class/family; .out cols: 5=query 6=start 7=end 9=strand(C/+) 10=repeat 11=class
awk 'NR>3 && $5!="" {
  start = $6 - 1; if (start < 0) start = 0;
  strand = ($9 == "C") ? "-" : "+";
  print $5"\t"start"\t"$7"\t"$10"\t0\t"strand          > "all_repeats.bed";
  print $5"\t"start"\t"$7"\t"$10"\t0\t"strand"\t"$11   > TMP "";
}' TMP="$TMP_CLASS" "$RM_OUT"

# split by class/family (class is col7 in the tmp file)
cut -f7 "$TMP_CLASS" | sort -u | while read -r CLASS; do
    SAFE=$(echo "$CLASS" | tr '/ ' '__')
    awk -v c="$CLASS" -F'\t' '$7 == c {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' \
        "$TMP_CLASS" > "repeats_${SAFE}.bed"
    echo "  ${CLASS}: $(wc -l < "repeats_${SAFE}.bed") instances"
done

rm -f "$TMP_CLASS"