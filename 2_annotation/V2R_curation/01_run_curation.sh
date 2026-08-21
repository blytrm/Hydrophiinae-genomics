#!/usr/bin/env bash
#SBATCH --job-name=v2r_curation
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --time=8:00:00
#SBATCH --mem=64GB
#SBATCH -o ./joblogs/%x_%j.log

# V2R curation for the hmaj EviAnn annotation.
# 01_run_curation.sh = evidence tracks -> candidate loci -> pHMM scan

set -euo pipefail
export LC_ALL=C   

SPECIES="${1:-hmaj}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="/scratchdata1/users/a1864358/sanders_lab/asm/files/annotation/eviann_results"
ANNO="/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/eviann_results"
V2RHMM="/hpcfs/users/a1864358/sanders_lab/v2r_hmm"
BASE="${ROOT}/${SPECIES}"
OUT="${HERE}/results/${SPECIES}"

require_file() {
    [[ -s "$1" ]] || { echo "ERROR: missing required input $1" >&2; exit 1; }
}

sort_bed() {
    sort -k1,1 -k2,2n "$@"
}

# REQUIRED INPUTS 
EVIANN_GFF="${BASE}/${SPECIES}_tidy.gff"
[[ -s "${EVIANN_GFF}" ]] || EVIANN_GFF="${ANNO}/${SPECIES}/${SPECIES}.fa.functional_note.pseudo_label.gff"
EVIANN_PROT="${ANNO}/${SPECIES}/${SPECIES}.fa.functional_note.proteins.fasta"
MINIPROT="${BASE}/miniprot.v2r.gff" # primary = spliced protein aln
PHMM_PROT="${BASE}/${SPECIES}.prot.domtbl" # primary = v2r pHMM vs proteome


# OPTIONAL INPUTS
G2OR_DIR="${BASE}/g2or" # control = non-V2R OR finder (genome2or)

# SPECIES-INDEPENDENT REFERENCES
DECOYS="${V2RHMM}/filtered/02_decoy_sequences_no_v2r.fasta" # decoys = non-V2R Class C seeds
V2R_HMM="${V2RHMM}/hmm/v2r_full.hmm" # pre-validated V2R profile
PFAM="/hpcfs/users/a1864358/sanders_lab/prog/pfam_db/Pfam-A.hmm" # source of the 7TM profile
SEVENTM_QUERY="${V2RHMM}/hma_7tm.fa" # fixed 98-query, 250aa 7TM domain slices; same set for every species
ASSEMBLY="/hpcfs/users/a1864358/sanders_lab/asm/files/final-asms/10-${SPECIES}-final.renamed.fa"

# MERGE
    # miniprot + pHMM intervals merged
MERGE_GAP="${MERGE_GAP:-1000}" #bp -> merge loci closer than
MIN_LEN="${MIN_LEN:-400}" # bp; drop merged regions shorter than


# PANEL
PANEL_DIR="${HERE}/panel" # shared class C HMM panel 
# /results/<species>/tracks/ = per species evidence BEDs
# /joblogs = slurm job logs
mkdir -p "${OUT}/tracks" "${PANEL_DIR}" "${HERE}/joblogs"

# files = gff annotation, protein fasta from annotation, miniprot primary: v2r proteins aligned to genome, pHMM primary: v2r profile vs the proteome, decoys (HMMS), validated v2r profile HMM
for f in "${EVIANN_GFF}" "${EVIANN_PROT}" "${MINIPROT}" "${PHMM_PROT}" \
         "${DECOYS}" "${V2R_HMM}" "${SEVENTM_QUERY}" "${ASSEMBLY}"; do
    require_file "$f"
done

TRACKS=(miniprot phmm_prot)
echo "species=${SPECIES} threads=${THREADS} merge_gap=${MERGE_GAP} min_len=${MIN_LEN}"


# read annotation gff -> from mRNA extract
#   intervals + gene ID -> eviann_mrna.bed
#   exon counts -> exon_counts.tsv
#   mrna attributes -> mrna_attrs.tsv
#   v2r labelled transcripts from EviAnn -> eviann_v2r_labelled.txt

awk -F'\t' -v out="${OUT}" 'BEGIN{OFS="\t"} $3=="mRNA" {
    id=""; gene=""; ev="NA"; sc="NA"; pc="NA"; ne="NA"; bt="NA"; note="NA"
    if (match($9, /Note=[^;]*/)) note = substr($9, RSTART+5, RLENGTH-5)
    n = split($9, kv, ";")
    for (i = 1; i <= n; i++) {
        split(kv[i], a, "="); key = tolower(a[1])
        if      (a[1] == "ID")           id   = a[2]
        else if (a[1] == "Parent")       gene = a[2]
        else if (key == "evidence")      ev   = a[2]
        else if (key == "startcodon")    sc   = a[2]
        else if (key == "stopcodon")     pc   = a[2]
        else if (key == "num_exons")     ne   = a[2]
        else if (a[1] == "gene_biotype") bt   = a[2]
    }
    if (id == "") next
    print $1, $4-1, $5, id, ".", $7, gene  > (out "/mrna.bed.tmp")
    print id, ev, sc, pc, ne, bt, note     > (out "/attrs.tmp")
    if (ne != "NA") print "Parent=" id, ne > (out "/exon_counts.tsv")
    if (tolower(note) ~ /vomeronasal|vmn2r|v2r/) print id > (out "/v2r.tmp")
}' "${EVIANN_GFF}"

sort_bed "${OUT}/mrna.bed.tmp" > "${OUT}/eviann_mrna.bed"
sort -u "${OUT}/v2r.tmp" > "${OUT}/eviann_v2r_labelled.txt"
{ printf 'transcript_id\tevidence\tstart_codon\tstop_codon\tnum_exons\tgene_biotype\tnote\n'
  cat "${OUT}/attrs.tmp"; } > "${OUT}/mrna_attrs.tsv"
rm -f "${OUT}"/mrna.bed.tmp "${OUT}"/attrs.tmp "${OUT}"/v2r.tmp
EXON_COUNTS="${OUT}/exon_counts.tsv"

echo "EviAnn mRNAs: $(wc -l < "${OUT}/eviann_mrna.bed")"
echo "EviAnn-labelled V2R transcripts: $(wc -l < "${OUT}/eviann_v2r_labelled.txt")"



# EVIDENCE TRACKS -> BED
#  PRIMARY (miniprot, phmm_prot) -> create locus
#  OTHER (tblastn / 7TM_PFAM) -> additional info
    # only annotate loci a primary track already found
# CONTROL (genome2or) -> flag loci that look like other non-V2R ORs instead

# miniprot: spliced alignments of validated V2R proteins
awk -F'\t' 'BEGIN{OFS="\t"} $3=="mRNA" {print $1, $4-1, $5, "miniprot", ".", $7}' \
    "${MINIPROT}" | sort_bed > "${OUT}/tracks/miniprot.bed"

# phmm_prot: domtbl targets are mRNA IDs, so coordinates come from the annotation
grep -v '^#' "${PHMM_PROT}" | awk '{print $1}' | sort -u > "${OUT}/tracks/phmm_prot.ids"
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {keep[$1]; next} ($4 in keep) {print $1,$2,$3,$4,".",$6}' \
    "${OUT}/tracks/phmm_prot.ids" "${OUT}/eviann_mrna.bed" > "${OUT}/tracks/phmm_prot.bed"



# SEVEN TRANSMEMBRANE DOMAIN
    # 2 independent detectors -> 2 cols each in v2r_master.tsv
    #   tBLASTn 7tm hit [yes/no ; genomic coords] -- own genome, catches loci EviAnn never modelled
    #   pfam    7tm hit [yes/no ; protein coords] -- candidates.faa only

# SKIP_TBLASTN7TM=1 -> skip the genome-wide search, leave hit=no / coords=NA
if [[ "${SKIP_TBLASTN7TM:-0}" == "1" ]]; then
    echo "tBLASTn 7TM: SKIPPED (SKIP_TBLASTN7TM=1)"
    : > "${OUT}/tblastn_7tm.outfmt6"
else
    BLASTDB="${OUT}/blastdb_genome"
    [[ -s "${BLASTDB}.ndb" ]] || makeblastdb -in "${ASSEMBLY}" -dbtype nucl \
        -out "${BLASTDB}" > "${OUT}/makeblastdb.log" 2>&1
    [[ -s "${OUT}/tblastn_7tm.outfmt6" ]] || tblastn -query "${SEVENTM_QUERY}" \
        -db "${BLASTDB}" -evalue 1e-3 -seg no -num_threads "${THREADS}" \
        -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' \
        -out "${OUT}/tblastn_7tm.outfmt6"
fi

# HSPs -> BED, then union span of overlapping HSPs per transcript
awk -F'\t' 'BEGIN{OFS="\t"} NF==12 {
    s=($9<$10)?$9:$10; e=($9<$10)?$10:$9
    print $2, s-1, e, $1, $11, ($9<$10?"+":"-")
}' "${OUT}/tblastn_7tm.outfmt6" | sort_bed > "${OUT}/tracks/tblastn_7tm.bed"
{
    printf 'transcript_id\ttblastn_7tm_hit\ttblastn_7tm_coords\n'
    bedtools intersect -wa -wb -a "${OUT}/eviann_mrna.bed" -b "${OUT}/tracks/tblastn_7tm.bed" \
      | awk -F'\t' 'BEGIN{OFS="\t"} {
            k=$4
            if (!(k in seen) || $9 < lo[k]) lo[k]=$9
            if (!(k in seen) || $10 > hi[k]) hi[k]=$10
            chrom[k]=$8; seen[k]=1
        } END {for (k in seen) print k, "yes", chrom[k]":"lo[k]"-"hi[k]}'
} > "${OUT}/tblastn_7tm_hits.tsv"
echo "tBLASTn 7TM: $(wc -l < "${OUT}/tblastn_7tm.outfmt6") HSPs, " \
     "$(( $(wc -l < "${OUT}/tblastn_7tm_hits.tsv") - 1 )) candidate transcripts hit"

# pfam 7tm_3 (PF00003): protein-space, candidates.faa only -- the blind spot tBLASTn covers
[[ -s "${OUT}/pfam_7tm.hmm" ]] || hmmfetch "${PFAM}" 7tm_3 > "${OUT}/pfam_7tm.hmm"
hmmsearch --cpu "${THREADS}" -E 1e-3 --domtblout "${OUT}/pfam_7tm.domtbl" \
    "${OUT}/pfam_7tm.hmm" "${OUT}/candidates.faa" > /dev/null 2> "${OUT}/pfam_7tm.err"

# miniprot -> in-frame stop codons and frameshifts in its alignments.
    # EviAnn gene models cannot carry an internal stop .. only emits complete ORFs (ATG..stop)
awk -F'\t' 'BEGIN{OFS="\t"} $3=="mRNA" {
    fs = 0; st = 0; idy = 0
    if (match($9, /Frameshift=[0-9]+/)) fs  = substr($9, RSTART+11, RLENGTH-11)
    if (match($9, /StopCodon=[0-9]+/))  st  = substr($9, RSTART+10, RLENGTH-10)
    if (match($9, /Identity=[0-9.]+/))  idy = substr($9, RSTART+9,  RLENGTH-9)
    print $1, $4-1, $5, fs, st, idy
}' "${MINIPROT}" | sort_bed > "${OUT}/tracks/miniprot_detail.bed"
{
    printf 'locus_id\tminiprot_alns\tminiprot_frameshifts\tminiprot_stops\tminiprot_max_identity\n'
    bedtools intersect -wa -wb -a "${OUT}/loci.bed" -b "${OUT}/tracks/miniprot_detail.bed" \
      | awk -F'\t' 'BEGIN{OFS="\t"} {n[$4]++; f[$4]+=$10; s[$4]+=$11; if ($12+0 > m[$4]) m[$4]=$12+0}
            END {for (l in n) print l, n[l], f[l], s[l], m[l]}'
} > "${OUT}/locus_miniprot.tsv"


# genome2or
    # outputs FASTA, not BED --> coordinated in headers 
        # "hit1_ch1_141116198_141115253_-_iso1" or "ch1_141494637_141493709_-_INTER"
HAS_CONTROL=0
if compgen -G "${G2OR_DIR}/genome2or_itera*_final_func_dna_ORs.fasta" > /dev/null; then
    HAS_CONTROL=1
    cat "${G2OR_DIR}"/genome2or_itera*_final_func_dna_ORs.fasta \
        "${G2OR_DIR}"/genome2or_itera*_final_pseu_ORs.fasta \
      | grep '^>' | sed 's/^>//' \
      | awk -F'_' 'BEGIN{OFS="\t"} {
            i = ($1 ~ /^hit[0-9]+$/) ? 2 : 1
            s=$(i+1); e=$(i+2); st=$(i+3)
            if (s !~ /^[0-9]+$/ || e !~ /^[0-9]+$/) next
            lo=(s<e)?s:e; hi=(s<e)?e:s
            print $i, lo-1, hi, "OR", ".", (st=="+"||st=="-") ? st : "."
        }' | sort_bed \
      | bedtools merge -i - -c 4,5,6 -o distinct,distinct,distinct > "${OUT}/tracks/genome2or.bed"
fi



# candidate loci 
    # miniprot + phmm_prot intervals merged within MERGE_GAP, then length filtered
PREFIX="$(echo "${SPECIES}" | tr '[:lower:]' '[:upper:]')V2R"
    # make id prefix -> hmaj -> HMAJV2R
cat "${OUT}/tracks/miniprot.bed" "${OUT}/tracks/phmm_prot.bed" \
  | sort_bed \
  | bedtools merge -d "${MERGE_GAP}" -i - \
  | awk -v p="${PREFIX}" -v m="${MIN_LEN}" 'BEGIN{OFS="\t"} $3-$2 >= m {
        n++; print $1, $2, $3, sprintf("%s_%s_%06d", p, $1, n), 0, "."
    }' > "${OUT}/loci.bed"

# - `cat` the 2 beds
	# - `sort_bed` -> sort by chrom, then start & `bedtools merge -d 1000 -i -` -> merge nearby/overlapping ... `-d 1000` => merge gap ... `-i -` = bread BED from stdin
		# - `| awk ...` -> drop short merges + assignn locus ID

echo "candidate loci: $(wc -l < "${OUT}/loci.bed")"



# which tracks support each locus
# which loci flagged as OR (control)

# keep loci that overlap the track, print locus ID + trackname
for t in "${TRACKS[@]}"; do
    bedtools intersect -u -a "${OUT}/loci.bed" -b "${OUT}/tracks/${t}.bed" \
        | awk -v t="$t" 'BEGIN{OFS="\t"} {print $4, t}'
done | sort -k1,1 -k2,2 > "${OUT}/locus_track_pairs.tsv" #  all pairs sorted

: > "${OUT}/locus_control.txt" # create/empty control list
if [[ "${HAS_CONTROL}" == 1 ]]; then
    bedtools intersect -u -a "${OUT}/loci.bed" \
        -b "${OUT}/tracks/genome2or.bed" | cut -f4 | sort -u > "${OUT}/locus_control.txt"
        # if genome2or -> intersect loci -> keep overlapping loci ids (cut -f4) -> unique sort 
fi

# collapse to one row per locus: locus_id, n_primary_tracks, tracks, control
awk -F'\t' 'BEGIN{OFS="\t"; print "locus_id","n_primary_tracks","tracks","control_overlap"}
    ARGIND==1 {ctrl[$1]=1; next}
    {tracks[$1] = ($1 in tracks) ? tracks[$1] "," $2 : $2
     if ($2=="miniprot" || $2=="phmm_prot") n[$1]++}
    END {for (l in tracks) print l, n[l]+0, tracks[l], (l in ctrl ? "genome2or" : "")}' \
    "${OUT}/locus_control.txt" "${OUT}/locus_track_pairs.tsv" > "${OUT}/locus_tracks.tsv"
echo "loci flagged by OR control: $(wc -l < "${OUT}/locus_control.txt")"
# locus_control.txt -> loci flagged as OR (control)
# locus_track_pairs.tsv -> sep=, tracks list; count miniprot/phmm_prot (n_primary_tracks)



# Candidate transcripts = every EviAnn mRNA overlapping a locus (not only
# pHMM hits). Column 4 = phmm_hit ->  win if a transcript spans two loci
{
    printf 'locus_id\ttranscript_id\tgene_id\tphmm_hit\n'
    bedtools intersect -wa -wb -a "${OUT}/loci.bed" -b "${OUT}/eviann_mrna.bed" \
      | awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {hit[$1]; next}
            {print $4, $10, $13, ($10 in hit) ? 1 : 0}' \
            "${OUT}/tracks/phmm_prot.ids" - \
      | sort -u
} > "${OUT}/locus_transcripts.tsv"
    # => every overlap 
    # Intersect candidate loci with every EviAnn mRNA (-wa -wb keeps both sides).
    # Join against phmm_prot.ids: if that transcript was a proteome-pHMM hit → phmm_hit=1, else 0.
    # Columns: locus_id, transcript_id, gene_id, phmm_hit
    # One transcript can appear on more than one locus if it spans a boundary.

    # one locus per transcript: pHMM-supported first, then lowest locus ID

{
    printf 'transcript_id\tlocus_id\n'
    awk -F'\t' 'BEGIN{OFS="\t"} NR>1 {print $2, $1, $4}' "${OUT}/locus_transcripts.tsv" \
      | sort -k1,1 -k3,3r -k2,2 \
      | awk -F'\t' 'BEGIN{OFS="\t"} !seen[$1]++ {print $1, $2}'
} > "${OUT}/mrna_locus.tsv"
    # one locus per transcript ... prefer phmm suppported loci ==> each transcript credited to 1 locus

{ awk -F'\t' 'NR>1 {print $2}' "${OUT}/locus_transcripts.tsv"
  cat "${OUT}/eviann_v2r_labelled.txt"; } | sort -u > "${OUT}/candidate_ids.txt"
seqkit grep -f "${OUT}/candidate_ids.txt" "${EVIANN_PROT}" > "${OUT}/candidates.faa"
echo "candidate transcripts: $(wc -l < "${OUT}/candidate_ids.txt")  proteins retrieved: $(grep -c '^>' "${OUT}/candidates.faa")"
    # candidate_ids.txt + candidates.faa -> what is scored
        # All overlapping transcript IDs, plus anything EviAnn already labelled V2R (even if it sits outside every locus).
        # pulls proteins from proteome -> candidates.faa




# Per-model quality evidence. None of this creates or supports a locus --
    # it only annotates one --> cannot change any count.

# v2r_domains.tsv -> v2r profile coverage from the pHMM domtbl
    # $1 = transcript id / $6 = v2r profile length / $16/17 = which part of profile was hit -> hmm_from / hmm_to
{
    printf 'transcript_id\tmodel_len\thmm_from\thmm_to\n'
    grep -v '^#' "${PHMM_PROT}" | awk 'BEGIN{OFS="\t"} {print $1, $6, $16, $17}'
} > "${OUT}/v2r_domains.tsv"


# non-V2R Class C panel
    # A V2R profile alone cannot say "V2R rather than GRM", so each Class C family
    # gets its own profile and assignment is argmax bitscore with a margin.
SEED_DIR="${PANEL_DIR}/seeds"
PANEL="${PANEL_DIR}/classC_panel.hmm"
if [[ -s "${PANEL}.h3i" ]]; then
    echo "Class C panel: reusing ${PANEL}"
else
    rm -rf "${SEED_DIR}"; mkdir -p "${SEED_DIR}"
    seqkit fx2tab "${DECOYS}" | awk -F'\t' -v dir="${SEED_DIR}" '
    BEGIN {
        IGNORECASE = 1
        split("GRM CASR GPRC6A TAS1R GABBR GPRC5", fams, " ")   # order = tie-break order
        pat["GRM"]    = "\\yGRM[1-8]\\y\t\\ymGluR[1-8]\\y\tmetabotropic glutamate"
        pat["CASR"]   = "\\yCASR\\y\tcalcium[- ]sensing"
        pat["GPRC6A"] = "\\yGPRC6A?\\y\treceptor family C group 6"
        pat["TAS1R"]  = "\\yTAS1R[1-3]\\y\ttaste receptor type 1"
        pat["GABBR"]  = "\\yGABBR[12]\\y\tGABA[-_ ]?B\\y\tgamma-aminobutyric acid type B\tgamma-aminobutyric acid \\(GABA\\) B"
        pat["GPRC5"]  = "\\yGPRC5[A-D]?\\y\treceptor family C group 5\tretinoic acid[- ]induced (protein )?3"
    }
    length($2) >= 200 {
        best = ""; bestscore = 0
        for (i = 1; i <= 6; i++) {
            f = fams[i]; n = split(pat[f], p, "\t"); score = 0
            for (j = 1; j <= n; j++) if ($1 ~ p[j]) score++
            if (score > bestscore) { bestscore = score; best = f }
        }
        if (best != "") print ">" $1 "\n" $2 > (dir "/" best ".faa")
    }'

    # One profile per family: cd-hit removes near-duplicates, mafft aligns, trimal
    # drops gappy columns, hmmbuild builds. -n names the profile after the family.
    : > "${PANEL}"
    printf 'family\tn_seeds\tn_after_cdhit\taln_len\n' > "${PANEL_DIR}/panel.key.tsv"
    for FAM in GRM CASR GPRC6A TAS1R GABBR GPRC5; do
        W="${PANEL_DIR}/_${FAM}"; mkdir -p "${W}"
        N_SEED=$(grep -c '^>' "${SEED_DIR}/${FAM}.faa")
        # -T 1: cd-hit's cluster order depends on thread scheduling, which shifts
        # profile bitscores by ~1 bit between runs. 
        cd-hit -i "${SEED_DIR}/${FAM}.faa" -o "${W}/cdhit.faa" -c 0.95 -n 5 \
               -T 1 -M 0 > "${W}/cdhit.log" 2>&1
        # --thread 1 for the same reason as cd-hit above: multi-threaded mafft breaks
        # alignment ties differently between runs.
        mafft --localpair --maxiterate 1000 --anysymbol --thread 1 \
              "${W}/cdhit.faa" > "${W}/aln.faa" 2> "${W}/mafft.log"
        trimal -gt 0.5 -cons 60 -in "${W}/aln.faa" -out "${W}/aln.trim.faa" 2> "${W}/trimal.log"
        hmmbuild --amino -n "${FAM}" --cpu "${THREADS}" "${W}/${FAM}.hmm" \
                 "${W}/aln.trim.faa" > "${W}/hmmbuild.log" 2>&1
        cat "${W}/${FAM}.hmm" >> "${PANEL}"
        printf '%s\t%s\t%s\t%s\n' "${FAM}" "${N_SEED}" \
            "$(grep -c '^>' "${W}/cdhit.faa")" \
            "$(awk '/^LENG /{print $2; exit}' "${W}/${FAM}.hmm")" >> "${PANEL_DIR}/panel.key.tsv"
    done

    # The V2R profile is the pre-validated one
    sed '0,/^NAME /s/^NAME .*/NAME  V2R_OlfC/' "${V2R_HMM}" >> "${PANEL}"
    printf 'V2R_OlfC\tprebuilt\tprebuilt\t%s\n' \
        "$(awk '/^LENG /{print $2; exit}' "${V2R_HMM}")" >> "${PANEL_DIR}/panel.key.tsv"
    hmmpress -f "${PANEL}" > /dev/null
fi
cat "${PANEL_DIR}/panel.key.tsv"

# Score every candidate protein against the whole panel at once.
hmmscan --cpu "${THREADS}" -E 1e-3 --tblout "${OUT}/classC.tblout" \
        "${PANEL}" "${OUT}/candidates.faa" > /dev/null 2> "${OUT}/classC.err"


# Tables
python3 "${HERE}/make_tables.py" "${OUT}" "${EXON_COUNTS}"

# Functional annotation,
if [[ -s /scratchdata1/users/a1864358/sanders_lab/prog/eggnog_db/eggnog.db ]]; then
    if [[ ! -s "${HERE}/eggnog/out/${SPECIES}/v2r.emapper.annotations" ]]; then
        bash "${HERE}/eggnog/run_emapper.sh" "${OUT}" "${THREADS}" "${SPECIES}"
    fi
    python3 "${HERE}/make_tables.py" "${OUT}" "${EXON_COUNTS}" > /dev/null
fi
python3 "${HERE}/plot_genome.py" "${OUT}" "${ASSEMBLY}.fai"

echo
echo "Done. Results in ${OUT}/"
    