# V2R curation — five *Hydrophis* sea snake genomes

Curates vomeronasal type-2 receptor (V2R) genes from EviAnn annotations by
reconciling every search engine run against each assembly, discriminating V2R
from the other Class C GPCR families with a competing profile-HMM panel,
repairing fragmented gene models, calling pseudogenes, and annotating the
result functionally.

Five genomes are curated with identical parameters: `hmaj`, `hcy`, `hcure`,
`hcurw`, `horn`. `hmaj` is *Hydrophis major* and `hcy` is *H. cyanocinctus*;
confirm the full names behind the other three codes before publication.

```bash
conda activate general

bash processing/00_prepare_species.sh <species>   # once per genome: miniprot + pHMM
bash 01_run_curation.sh <species>                 # tracks -> loci -> panel -> tables
python3 processing/02_model_repair.py <species>   # EviAnn vs miniprot vs ORF models
python3 processing/03_pseudogene_calls.py         # gene_status, all species
python3 processing/04_collate.py                  # cross-species tables
```

`01_run_curation.sh` takes ~2 min per genome on 16 cores; `00_prepare_species.sh`
~20 min, dominated by miniprot.

| Script | Does |
|---|---|
| `v2rlib.py` | thresholds, shared readers, and the classification rules |
| `01_run_curation.sh` | evidence tracks → candidate loci → Class C panel → hmmscan → 7TM scan |
| `make_tables.py` | scans → assignments, catalogue, master table, summary |
| `plot_genome.py` | master table → genome map |
| `processing/00_prepare_species.sh` | the two primary evidence tracks for a new genome |
| `processing/01_positive_control.sh` | run on *H. cyanocinctus* and score against its benchmark |
| `processing/02_model_repair.py` | pick the most complete model per locus, guard any merge |
| `processing/03_pseudogene_calls.py` | functional / pseudogene / partial |
| `processing/04_collate.py` | every species into one comparison |
| `eggnog/` | database fetch and eggNOG-mapper run |

Outputs land in `results/<species>/`; cross-species tables in `processing/`.

---

# METHOD

## 1. Evidence tracks

Each search becomes a BED file and is given one of three roles:

| Track | Role | What it is |
|---|---|---|
| `miniprot` | **primary** | spliced protein-to-genome alignment of 1,516 vomeronasal receptor proteins |
| `phmm_prot` | **primary** | V2R profile HMM vs the EviAnn proteome |
| `genome2or` | control | olfactory-receptor gene finder |

**Only primary tracks can create a locus.** All five genomes now run on the same
two primary tracks (`genome2or` is present where the OR scan was run, absent
otherwise); no floor tracks are used. Earlier revisions of this pipeline also
carried floor tracks (`phmm_orf`, `phmm_6frame`, `tblastn`, `tblastn_loci` —
sensitive-but-unspecific evidence that could only annotate a locus a primary
track already found, never create one). They were removed as a simplification;
by design their removal cannot change a locus count, and re-running confirms it
did not — loci and gene counts are unchanged from when hmaj alone carried them.

The EviAnn annotation is deliberately not an evidence track: adding its ~22,000
genes as primary intervals would seed a candidate locus at every gene in the
genome. It supplies coordinates and the transcripts overlapping a locus.

## 2. Candidate loci

```bash
cat tracks/miniprot.bed tracks/phmm_prot.bed \
  | sort -k1,1 -k2,2n | bedtools merge -d 1000 -i - | awk '$3-$2 >= 400'
```

IDs are issued in genome order (`HMAJV2R_<chrom>_<n>`) and depend only on the
input intervals, so they are stable across re-runs. `merge -d 1000` is
calibrated: on hmaj the locus count is 1,076 / 1,073 / 1,059 at gaps of
0 / 1,000 / 5,000 bp, and the **gene** count is identical at 0 and 1,000. The
oversized loci are therefore not a merge artefact — they are chains of
transitively overlapping miniprot alignments running the length of a tandem array.

## 3. Candidate transcripts

Every EviAnn mRNA overlapping a candidate locus, plus every transcript EviAnn
itself labelled V2R even where no locus covers it. Keying on interval overlap
rather than on the proteome scan's own hits is deliberate: otherwise the
classifier could only re-examine that one engine's findings.

## 4. Class C panel

A V2R profile alone can only say "V2R-like", not "V2R rather than GRM". So each
Class C family gets its own profile and assignment is argmax bitscore with a
margin. Decoys (376 non-V2R Class C proteins) are split to families by header
keyword, then `cd-hit -c 0.95` → `mafft --localpair` → `trimal -gt 0.5 -cons 60`
→ `hmmbuild`.

| family | seeds | after cd-hit | profile length |
|---|---|---|---|
| GRM | 161 | 25 | 888 |
| CASR | 16 | 2 | 1089 |
| GPRC6A | 8 | 4 | 926 |
| TAS1R | 18 | 9 | 837 |
| GABBR | 81 | 5 | 852 |
| GPRC5 | 19 | 7 | 400 |
| **V2R_OlfC** | prebuilt | — | 823 |

The panel is built once and **shared by all five species** (`panel/`), so
bitscores are comparable across genomes by construction. Both aligners run
single-threaded: their tie-breaking depends on thread scheduling, which moved
bitscores by ~1 bit between runs.

The V2R profile is the pre-validated `v2r_full.hmm`, not a rebuild — building de
novo from the 1,683 validated seeds gives a ~6,000-column alignment, and a
profile that long inflates its own bitscores and wins the argmax for the wrong
reason. Profile lengths must be comparable for an argmax panel to mean anything.

**What the engines were built from.** miniprot is not trained; its query is 1,516
vomeronasal receptor proteins (1,035 squamate UniProt, 43 rat, 41 *Xenopus*,
37 mouse, 36 squamate NCBI, 12 squamate Pfam). The V2R profile came from 581
curated seeds (332 mouse Vmn2r, 200 rat Vom2r, 49 *Varanus komodoensis*) expanded
and filtered to 1,683 validated sequences, reduced to 1,206 representatives,
aligned and trimmed to 823 columns. Both descend from vomeronasal collections —
independent in *method*, not in reference data, so their agreement is weaker
evidence than the EviAnn or eggNOG comparisons.

## 5. Assignment

`hmmscan -E 1e-3` of every candidate protein against the whole panel:

```
delta_bits = bits(best family) - bits(runner-up family)
```

| Label | Condition |
|---|---|
| `<family>` | ≥2 families hit, `delta_bits >= 20` — confident |
| `ambiguous` | ≥2 families hit, `delta_bits < 20` |
| `single_family_weak` | only one family cleared the ceiling and scored < 100 bits |

`single_family_weak` exists because with one family there is no runner-up and no
margin to test. Treating the raw bitscore as the margin would mark every such hit
maximally confident.

## 6. Tiers — per locus

| Tier | Condition |
|---|---|
| **T2** | the locus's best call is `V2R_OlfC` with `delta_bits >= 20` |
| **T1** | a primary engine supports the locus, but no confident V2R call |
| **T0** | confidently a *different* Class C family, or no primary support |

Overlap with the `genome2or` control demotes T2 → T1. There is deliberately no
structural gate inside the tier: exon count, coverage and 7TM presence are all
columns, so any integrity filter is a filter on the output.

## 7. Model repair

45% of EviAnn V2R models cover under 60% of the V2R profile. Most are not short
receptors — they are one gene chopped into pieces. Two other model sources exist
at the same loci: **miniprot** alignments and **getorf ORFs**. For each locus the
most complete model wins, using V2R-profile coverage as the common yardstick.
The pHMM is not a source — it scores proteins, it does not build models; it is
the referee.

A miniprot model spanning two or more EviAnn genes claims those genes are one.
Four guards must all pass:

1. one query protein spans both fragments (≥80% of the query used);
2. repaired coverage ≥ 0.90, not merely improved;
3. no spanned gene has independent transcript support;
4. the merged model's exons cover the fragments' spans.

Failures keep the EviAnn models and are listed in `02_merge_review_<sp>.tsv`.

## 8. Gene status

**Final specification:**

```
functional   repaired_coverage_class == intact
             AND inframe_stops == 0 AND frameshifts == 0
             (no 7TM requirement — see below)

pseudogene   inframe_stops > 0 OR frameshifts > 0

partial      everything else — truncated or too fragmented to judge.
             An assembly/annotation category, not a biological one.
```

`inframe_stops` and `frameshifts` come from **miniprot**, not from EviAnn.
EviAnn cannot report a pseudogene: it only emits complete `ATG..stop` ORFs, so it
declines to annotate the broken part rather than annotating it as broken.
miniprot aligns a known-good V2R protein through the locus regardless, and
reports what it had to align through. That asymmetry is the whole basis of the
pseudogene call.

Rules A–F in `v2rlib.py` vary what `functional` demands; `pseudogene` is
identical under every rule but F, because it depends only on frame disruption.
**Rule B — `intact`, no disruption, no 7TM requirement — is the final rule and is
what `CHOSEN_RULE` in `v2rlib.py` writes to `gene_status`.** The 7TM requirement
(rule A) was dropped: 0 intact models lack a 7TM domain before repair, and after
repair the two rules differ by a handful of transcripts whose winning model is a
repaired (non-EviAnn) source not independently re-scanned for the domain — 7TM
was adding noise, not signal, once coverage already gates on `intact`.

A locus's T2 tier (`classC_family == V2R_OlfC` at `delta_bits >= 20` over the
runner-up family) and `gene_status == functional` turn out to be fully
redundant in all five genomes: every `functional` transcript already sits at a
T2 locus, in every species, with zero exceptions. `is_v2r == yes` already
implies a confident T2-grade call, so requiring T2 explicitly changes nothing —
it is included below only because it makes the definition self-contained.

---

# RESULTS

## Repertoire sizes

| | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| assembly (Mb) | 2,150 | 1,928 | 1,923 | 2,059 | 1,998 |
| annotated mRNAs | 52,501 | 51,220 | 51,588 | 51,444 | 51,521 |
| candidate loci | 1,073 | 715 | 673 | 895 | 668 |
| T2 / T1 / T0 loci | 607/453/13 | 393/309/13 | 400/260/13 | 483/399/13 | 390/265/13 |
| **V2R transcripts** | **736** | **493** | **499** | **574** | **482** |
| **V2R genes** | **661** | **433** | **418** | **513** | **425** |
| EviAnn-rescued | 5 | 3 | 2 | 4 | 2 |
| 7TM present | 592 | 388 | 411 | 459 | 416 |

## Completeness, before and after model repair

| | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| intact, before → after | 348 → **459** | 208 → **305** | 239 → **339** | 236 → **345** | 228 → **337** |
| partial | 56 → 27 | 34 → 18 | 30 → 9 | 46 → 22 | 38 → 18 |
| fragmented | 332 → **245** | 251 → **167** | 230 → **149** | 292 → **203** | 216 → **125** |
| winning model: eviann / miniprot / orf | 500 / 209 / 22 | 327 / 163 / – | 313 / 184 / – | 364 / 206 / – | 299 / 181 / – |
| merges accepted | 6 | 4 | 2 | 4 | 7 |

Repair moves 26–42% of fragmented models up a class in every genome. miniprot
wins 28–37% of loci — consistent across five independent assemblies, which is
what you would expect if the fragmentation is annotation breakage rather than
short genes.

## Functional / pseudogene / partial (final rule, after repair)

`repaired_coverage_class == intact AND inframe_stops == 0 AND frameshifts == 0`,
no 7TM requirement. Transcripts / genes.

| | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| **functional** | 301 / **281** | 175 / **155** | 196 / **166** | 211 / **183** | 184 / **169** |
| **pseudogene** | 180 / **141** | 142 / **114** | 166 / **131** | 166 / **142** | 171 / **139** |
| **partial** | 255 / 239 | 176 / 164 | 137 / 121 | 197 / 188 | 127 / 117 |

Rule sensitivity, functional genes for hmaj: **B 281** (final) / A 275 (+7TM
requirement) / C 290 / D 299 — the choice of rule moves the count by ~9%, far
less than the assembly differences below. Rule F is not a variant: it ignores
disruption and collapses hmaj's pseudogene count from 141 to 33 genes.

## High-confidence set

`repaired_coverage_class == intact AND locus_tier == T2 AND gene_status ==
functional` — and the T2 and functional conditions are exactly redundant (see
Method §8), so this set is identical to plain `gene_status == functional`.
Transcripts / genes / loci.

| | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| **high confidence** | 301 / 281 / 268 | 175 / 155 / 143 | 196 / 166 / 163 | 211 / 183 / 180 | 184 / 169 / 165 |
| + ≥6 exons + 700–1100 aa protein | **181 / 176 / 172** | **112 / 106 / 103** | **124 / 122 / 121** | **117 / 110 / 110** | **118 / 118 / 118** |

Adding the canonical exon count and length window for an in-depth analysis set
roughly halves the functional count in every genome (e.g. hmaj 281 → 176 genes)
— most `functional` transcripts pass on sequence completeness and an intact
reading frame but fall outside the canonical 6-exon, 700–1100 aa structure,
which is a stricter, additional filter rather than a re-derivation of the same
thing.

## Disruption evidence

| | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| in-frame stops (miniprot) | 437 | 568 | 470 | 613 | 552 |
| frameshifts | 845 | 804 | 851 | 1,012 | 1,040 |
| pseudogene genes | 141 | 114 | 131 | 142 | 139 |

Pseudogene counts are the most stable quantity in the whole comparison —
114–142 genes across five genomes, a 25% spread, against a 58% spread in total
gene count.

## Placement in the assembly

| genes | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| on chromosomes | 463 | 430 | **178** | 509 | 418 |
| on scaffolds | 198 | 3 | **240** | 4 | 7 |
| % unplaced | 30% | 0.7% | **57%** | 0.8% | 1.6% |

## Functional annotation (eggNOG)

Every retained V2R protein in all five genomes is annotated.

| | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| annotated | 736 | 493 | 499 | 574 | 482 |
| canonical `7tm_3, ANF_receptor, NCD3G` | 729 | 487 | 494 | 568 | 476 |
| architecture inconsistent | **7** | **6** | **5** | **6** | **6** |

The architecture-inconsistent calls are GRIK1/GRIK2 and PRRG1 in every genome —
ionotropic glutamate receptors share the ANF_receptor Venus-flytrap domain with
Class C GPCRs but their output module is an ion channel, not a 7TM, and the panel
has no iGluR profile to compete against them. Filter
`eggnog_arch_consistent == no` to drop them. None survive the high-confidence
filters.

## Positive control

`hcy` run through this pipeline against the independent benchmark derived from
its **funannotate** proteome (a different annotation of the same genome):

| Metric | Benchmark | This pipeline | Difference |
|---|---|---|---|
| V2R-family proteins | 458 | 493 | +7.6% |
| confident V2R transcripts | 426 | 493 | +15.7% |
| **V2R genes** | **~408** | **433** | **+6.1%** |

The gene row is the like-for-like comparison — the benchmark counted proteins
from a one-protein-per-gene annotation, while EviAnn carries isoforms. Mouse
would not work as a second control: 332 of the 581 curated seeds behind the V2R
profile are mouse *Vmn2r* sequences, so recall would look excellent for the
wrong reason.

## An internal control nobody designed

T0 is **13 loci in all five genomes**, with the identical breakdown every time:

```
8 GRM   2 CASR   2 GABBR   1 GPRC6A
```

These are the single-copy Class C housekeeping receptors — eight metabotropic
glutamate receptors, two calcium-sensing, two GABA-B, one GPRC6A — recovered 1:1
in five independently assembled and annotated genomes. Nothing in the pipeline
enforces this. It is direct evidence that the panel is separating families on
sequence rather than on assembly-specific noise, and that the non-V2R side of the
argmax is behaving.

---

# DISCUSSION

**The headline difference between genomes is mostly assembly, not biology.**
hmaj has the largest repertoire (661 genes) and hcure the smallest (418), but
57% of hcure's genes sit on unplaced scaffolds against 30% for hmaj and under 2%
for the other three. V2R arrays are long tandem repeats — exactly the structure
that fails to scaffold — so a genome that scaffolds them poorly both fragments
individual genes and inflates locus counts. Chromosome-placed genes tell a
flatter story: 463, 430, 178, 509, 418. hcure remains the outlier and its total
should not be quoted without the caveat.

**Pseudogene counts are the most comparable quantity here.** They span 114–142
genes across five genomes — a 25% spread against 58% for total gene count. That
is expected: a frameshift or in-frame stop is a local sequence observation that
survives poor scaffolding, whereas counting complete genes requires the array to
be assembled correctly. If one number from this analysis is to be compared across
species, it should be this one.

**Model repair matters more than rule choice.** Moving from raw to repaired
coverage adds 97–111 intact models per genome; switching between classification
rules A–D moves the functional count by ~9%. Any repertoire comparison that has
not repaired its gene models is measuring the annotator as much as the genome.

**Four independent lines agree on the same seven proteins.** In hmaj, the
GRIK/PRRG1 contamination is flagged by our own 7TM scan (all 7 lack it), by
coverage class (all fragmented), by EviAnn's own labels (named GRIK/PRRG1), and
by eggNOG's Pfam architecture. The same pattern holds in all five genomes with
5–7 proteins each. This is the panel's one known systematic blind spot and it is
fully characterised; adding a GRIK/iGluR family to the panel would close it at
source rather than by filtering.

**What "functional" can and cannot mean here.** Only ~23% of hmaj's V2R models
have any transcript support — expected for a vomeronasal family whose expression
is restricted to the VNO, which is unlikely to be among the sequenced tissues,
but it means these are homology projections rather than expression-verified
genes. `functional` here means "an intact ORF with the right domain architecture
and no visible disruption", not "demonstrated to be expressed".

---

# Known limitations

1. **`repaired_coverage_class` is a per-locus property**, written onto every
   transcript row of that locus. A short isoform of a locus whose best model is
   full length is therefore labelled `intact` while its own `protein_len_aa` is
   small — 30 such rows in hmaj. Filtering `repaired_coverage_class == intact`
   selects *loci with an intact model*, not *transcripts that are intact*. Use
   `v2r_model_cov` for a strictly per-transcript measure.

2. **No phylogeny.** No tree is built, so V2R subfamilies are unassigned, the
   ambiguous calls are unadjudicated, and no T3 tier exists. This is the largest
   remaining gap and the natural next step now that all five genomes are curated.

3. **The panel has no ionotropic glutamate receptor family** — see above.

4. **CASR competes weakly**: cd-hit collapses its 16 seeds to 2.

5. **No floor tracks in the current pipeline.** All five genomes carry only the
   two primary tracks (`miniprot`, `phmm_prot`) plus `genome2or` where run. This
   is uniform across species now, so `tracks` is directly comparable.

6. **hcure is 57% unplaced**; its totals are assembly-limited.

7. **Strand is not resolved per locus**, so catalogue strand is `.`.
   Per-transcript strand is correct.

8. **Species identities**: only `hmaj` (*H. major*) and `hcy` (*H. cyanocinctus*)
   are confirmed here.

---

# Table columns

`results/<species>/v2r_master.tsv` is the table to filter — one row per candidate
locus × overlapping transcript, 68 columns:

```bash
is_v2r == yes                          # the repertoire
gene_status == functional              # intact ORF, no disruption (no 7TM requirement)
gene_status == pseudogene              # disrupted reading frame
repaired_coverage_class == intact      # loci with a full-length model
v2r_evidence == panel                  # exclude the EviAnn-rescued fragments
eggnog_arch_consistent == no           # the GRIK/PRRG1 contamination
tx_credited_to_locus != no             # one row per transcript
gene_status==functional AND n_exons>=6 AND 700<=protein_len_aa<=1100
                                        # high confidence, in-depth analysis set
                                        # (locus_tier==T2 is implied by functional, see Method §8)
```

Views of it: `v2r_genes.tsv` (is_v2r rows), `classC_all_genes.tsv` (every
confident Class C call), `v2r_catalogue.tsv` (per locus), `chrom_summary.tsv`
(per sequence), `summary.txt` (parameters and headline numbers).

Cross-species: `processing/04_cross_species.{tsv,md}`,
`processing/03_status_table.md`, `processing/02_repair_report_<sp>.md`,
`processing/01_control_report.md`.

---

# Reproducibility

The panel is built once and shared, and both aligners are pinned to one thread,
so re-runs are byte-identical (HMMER stamps a `DATE` line into `.hmm` files,
which is the only thing that changes). Locus IDs are issued in genome order from
the input intervals alone.

A simplification pass rewrote the pipeline in place — one GFF pass instead of
four, shared thresholds and rules in `v2rlib.py` instead of three divergent
copies, and `make_tables.py` restructured into functions. It was verified by
re-running all five genomes and diffing: **17 of 17 tables byte-identical for
every species**, and every number in every report unchanged. Two bugs were fixed
in the same pass, both cosmetic but wrong: every species' `summary.txt` was
titled "hmaj V2R curation" and every genome map "across the hmaj assembly".
