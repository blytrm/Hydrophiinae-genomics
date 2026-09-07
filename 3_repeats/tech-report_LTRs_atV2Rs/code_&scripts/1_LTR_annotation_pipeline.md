# 01 - Repeat annotation + filtering

- **Why EDTA first:** build a species-informed TE library (`EDTA.TElib.fa`) including intact LTRs, then feed that library to RepeatMasker for genome-wide coords.
- **Why LAI:** assembly-quality index from intact vs all LTR content; used later as a confound / mask check, not as the primary biology metric.
- **Why restrict to ch2 + chZ:** V2R arrays concentrate on the autosome (ch2) and sex chromosome (chZ); chromosomes are never pooled for H1

---

## 1. EDTA (TE library + intact LTRs)

```bash
#!/bin/bash
#SBATCH --job-name=edta_10horn
#SBATCH --partition=dm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180G
#SBATCH --time=8-00:00:00
#SBATCH --output=edta_%j.log
#SBATCH --error=edta_%j.err

set -euo pipefail

source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
conda activate EDTA

PROJ=/scratchdata1/users/a1864358/sanders_lab/repeatdensity-atv2rs
GENOME=$PROJ/10-horn-final.renamed.fa
CURATED=$PROJ/repeats/hydrophis_ornatus-garvin.fa.mod.EDTA.TElib.fa
THREADS=36
RUNDIR=$PROJ/edta_run

mkdir -p "$RUNDIR"
cd "$RUNDIR"
ln -sf "$GENOME" .
GBASE=$(basename "$GENOME")

EDTA.pl \
  --genome "$GBASE" \
  --species others \
  --step all \
  --curatedlib "$CURATED" \
  --sensitive 1 \
  --anno 1 \
  --threads "$THREADS" \
  --overwrite 0
```

- Curated Garvin *Hydrophis* lib seeds annotation (better than de novo-only for this clade)
- Canonical HPC tree (not shipped here): `/hpcfs/.../asm/files/repeat-anno/edta_10-horn/`.

---

## 2. LAI (LTR Assembly Index)

```bash
#!/bin/bash
#SBATCH --job-name=lai_10horn
#SBATCH --partition=dm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=80G
#SBATCH --time=2-00:00:00
#SBATCH --output=lai_%j.log
#SBATCH --error=lai_%j.err

set -euo pipefail

source /hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh
conda activate EDTA

PROJ=/scratchdata1/users/a1864358/sanders_lab/repeatdensity-atv2rs
RUN=$PROJ/edta_run
cd "$RUN"

G=10-horn-final.renamed.fa.mod
GENOME=$PROJ/10-horn-final.renamed.fa                 
PASS=$RUN/$G.EDTA.raw/LTR/$G.pass.list                
RMOUT=$RUN/$G.EDTA.anno/$G.out                        

LAI -genome "$GENOME" -intact "$PASS" -all "$RMOUT" -t 36

echo "DONE -> result: $RUN/$G.out.LAI"
```

- Scores how complete LTR annotation looks relative to intact elements.

---

## 3. RepeatMasker + class BED split

```bash
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
```

- Class column feeds `te_classification.tsv`

---

## 4. Workspace filtering (C2/CZ + LTR subset)

```r
# From scripts/v2r_priority_analysis.R — essential steps

# Restrict repeats to ch2 + chZ; join TE class; keep LTR/*
repeats_dt <- fread("repeats_c2cZ.bed", ...)  # 0-based BED
repeats_dt[, start := as.integer(start) + 1L]
repeats_dt <- merge(repeats_dt, te_class, by = "te_id", all.x = TRUE)
repeats_dt[chr == "ch6", chr := "chZ"]
repeats_dt <- repeats_dt[chr %in% c("ch2", "chZ")]

repeats_ltr_dt <- repeats_dt[grepl("^LTR/", te_class)]

# Gaps as permutation mask (Null B / C base)
mask_gaps <- make_gr(fread("gaps.bed", ...), bed0 = TRUE)
mask_c2cZ <- mask_gaps[seqnames(mask_gaps) %in% c("ch2", "chZ")]
```

- `**LTR/unknown`:** dominant near V2Rs but short fragments (~solo LTRs / fossils), not intact Gypsy/Copia — validated in report (see `figures/ltr_unknown_length_distribution.png`).

---

## Outputs used downstream

| Asset                   | Role                                |
| ----------------------- | ----------------------------------- |
| `te_classification.tsv` | TE id → class/family                |
| `gaps.bed`              | assembly-gap mask                   |
| `repeats_c2cZ.bed`      | all repeat intervals on C2/CZ       |
| `gc_1kb.bedgraph`       | GC for Null C match-then-test       |
| EDTA intact / LAI       | confound + intact-LTR side analyses |


