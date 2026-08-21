# V2R repertoires across five *Hydrophis* assemblies

All five sea snakes curated with one pipeline, one Class C profile panel, one set
of thresholds. This document is the method, the numbers, and — where the numbers
disagree with each other — what is causing it.

**Species:** `hmaj` (*H. major*), `hcy` (*H. cyanocinctus*), `hcure`, `hcurw`
(*H. curtus*, two assemblies), `horn` (*H. ornatus*).

---

# HEADLINE

| | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| **V2R transcripts** | **736** | 493 | 499 | 574 | 482 |
| **V2R genes** | **661** | 433 | 418 | 513 | 425 |
| non-redundant at 99% identity | **485** | 379 | 375 | 443 | 383 |
| functional (rule A) | **275** | 153 | 165 | 180 | 168 |
| pseudogene | 141 | 114 | 131 | 142 | 139 |
| partial | 247 | 170 | 128 | 194 | 123 |
| high confidence | **237** | 148 | 182 | 161 | 174 |

**Between 418 and 661 V2R genes per snake — but that spread is not all biology.**
The rest of this document is mostly about which parts of it are.

The single most defensible cross-species number is the **99%-non-redundant count
(375–485)**, because it is the one that does not move when an assembly retains
uncollapsed heterozygosity. The raw gene count spread is 1.58× (418→661); the
non-redundant spread is 1.29× (375→485).

![cross-species comparison](05_cross_species.png)

---

# METHOD

Identical for every species. Only two inputs are species-specific.

## What goes in

| Input | Per species? | Source |
|---|---|---|
| EviAnn annotation (GFF) + proteome | yes | existing EviAnn runs |
| miniprot alignment of 1,516 V2R proteins | yes — generated here | `processing/00_prepare_species.sh` |
| V2R profile vs proteome (domtbl) | yes — generated here | same script |
| V2R profile HMM (823 columns) | **no** | `v2r_hmm/hmm/v2r_full.hmm` |
| Class C decoy panel | **no** | built once into `panel/`, shared |
| Pfam 7tm_3, eggNOG 5.0.2 | **no** | shared references |

Floor tracks (getorf ORFs, 6-frame, tBLASTn) and the genome2or control exist only
for hmaj. They are **optional by design**: a floor track can annotate a locus a
primary engine already found but can never create one, so their absence cannot
change a locus count. Confirmed empirically — hcy, hcure, hcurw and horn ran on
the two primary engines alone.

**The panel is now built once and shared.** Previously each species rebuilt it,
which produced identical profiles but only by coincidence; sharing makes
cross-species bitscores comparable by construction and removed 55 MB of duplicated
HMM files.

## The chain

```
miniprot alignments  ┐
                     ├─ merge within 1 kb, keep ≥400 bp ──► candidate loci
V2R pHMM vs proteome ┘                                          │
                                                                ▼
                          every EviAnn mRNA overlapping a locus (+ EviAnn's own
                          V2R-labelled transcripts, wherever they sit)
                                                                │
                                                                ▼
                   hmmscan vs the 7-family Class C panel ──► family argmax
                                                                │
        ┌───────────────────────────────────────────────────────┤
        ▼                        ▼                              ▼
  tier (per locus)     coverage class (per model)      gene_status (per transcript)
  T2 / T1 / T0         intact / partial / fragmented   functional / pseudogene / partial
```

Three orthogonal axes, deliberately kept separate:

- **Tier** answers *is this locus a V2R?* — a sequence question.
- **Coverage class** answers *is the model complete?* — an annotation-quality question.
- **gene_status** answers *is the gene intact?* — a biological question.

Conflating them is the standard way repertoire sizes stop being comparable
between papers.

## Definitions

**Tiers** — per locus, first match wins:

| Tier | Condition |
|---|---|
| T2 | best call is `V2R_OlfC` with `delta_bits >= 20` over the runner-up family |
| T1 | a primary engine supports the locus, but no confident V2R call |
| T0 | confidently a *different* Class C family, or no primary support |

Override: overlap with the genome2or olfactory-receptor control demotes T2 → T1
(hmaj only; no other species has that track).

**Coverage classes** — per model, on repaired coverage:

```
intact       repaired_model_cov >= 0.80
partial      0.60 <= repaired_model_cov < 0.80
fragmented   repaired_model_cov < 0.60
```

`repaired_model_cov` = fraction of the 823-column V2R profile covered by the best
of the EviAnn, miniprot and ORF models at that locus.

**gene_status** — per transcript:

```
functional   repaired_coverage_class == intact
             AND has_7tm == yes
             AND inframe_stops == 0 AND frameshifts == 0
pseudogene   inframe_stops > 0 OR frameshifts > 0
partial      everything else
```

`inframe_stops` and `frameshifts` come from miniprot, not EviAnn. **EviAnn cannot
report a pseudogene**: it only emits complete `ATG..stop` ORFs, so it declines to
annotate the broken part rather than annotating it as broken. miniprot aligns a
known-good V2R protein through the locus regardless, and reports what it had to
align through. That asymmetry is the whole basis of the pseudogene call.

---

# RESULTS

## Assemblies — the context for everything below

| | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| assembly size (Mb) | **2150** | 1928 | 1923 | 2059 | 1998 |
| sequences | **1111** | 48 | 361 | 31 | 196 |
| % of bases in chromosomes | **92.4** | 99.9 | 94.4 | 99.9 | 99.5 |
| annotated mRNAs | 52,501 | 51,220 | 51,588 | 51,444 | 51,521 |

Annotation depth is essentially constant (51.2–52.5 k mRNAs). Assembly quality is
not: **hmaj is both the largest assembly and by far the most fragmented**, and
hcure is the second most fragmented. Hold that thought.

## Evidence recovered

| | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| miniprot alignments | 2929 | 2674 | 2863 | 2993 | 2876 |
| proteome pHMM hits | **939** | 633 | 630 | 729 | 607 |
| candidate loci | **1073** | 715 | 673 | 895 | 668 |

**Cause and effect.** miniprot alignment counts are nearly flat (2674–2993) — the
genome-level V2R signal is similar in all five. The proteome pHMM hits are not
(607–939), and candidate loci track the pHMM hits almost exactly. So the locus
spread is driven by **how many V2R-like proteins EviAnn annotated**, not by how
much V2R-like sequence exists. hmaj has 40–55% more annotated V2R-like proteins
than the others from a comparable amount of genomic signal.

## Tiers

| | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| T2 | 607 | 393 | 400 | 483 | 390 |
| T1 | 453 | 309 | 260 | 399 | 265 |
| T0 | 13 | 13 | 13 | 13 | 13 |

**T0 is 13 in every single species.** That is not a coincidence and it is a useful
sanity check: T0 is dominated by loci confidently assigned to a *different* Class
C family (GRM, GABBR, CASR, GPRC6A), and those are single-copy housekeeping
receptors present once per genome. Recovering exactly 13 five times independently
says the panel's non-V2R arm is behaving identically across species.

T1/T2 ratio is 0.75 (hmaj), 0.79 (hcy), 0.65 (hcure), 0.83 (hcurw), 0.68 (horn) —
no species is unusually full of unresolvable loci.

## The V2R set

| | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| transcripts | 736 | 493 | 499 | 574 | 482 |
| genes | 661 | 433 | 418 | 513 | 425 |
| loci | 608 | 393 | 400 | 483 | 390 |
| isoforms per gene | 1.11 | 1.14 | 1.19 | 1.12 | 1.13 |
| EviAnn-rescued | 5 | 3 | 2 | 4 | 2 |
| 7TM present | 592 | 388 | 411 | 459 | 416 |
| eggNOG architecture inconsistent | 7 | 6 | 5 | 6 | 6 |

Isoform rate is constant, so the transcript spread is not an isoform artefact —
it is a gene-count spread.

**eggNOG flags 5–7 architecture-inconsistent calls in every species.** These are
the same contamination in all five: GRIK1/GRIK2 kainate receptors and PRRG1. They
share the ANF_receptor Venus-flytrap domain with Class C GPCRs but their output
module is an ion channel, not a 7TM, and the panel has no ionotropic glutamate
receptor profile to compete against them. Consistent across species, small, and
removable with `eggnog_arch_consistent == no`.

## Gene-model repair

| | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| intact before → after | 348 → **459** | 208 → **305** | 239 → **339** | 236 → **345** | 228 → **337** |
| partial | 56 → 27 | 34 → 18 | 30 → 9 | 46 → 22 | 38 → 18 |
| fragmented | 332 → **245** | 251 → **167** | 230 → **149** | 292 → **203** | 216 → **125** |
| winning model: eviann | 500 | 327 | 313 | 364 | 299 |
| winning model: miniprot | 209 | 163 | 184 | 206 | 181 |
| winning model: orf | 22 | – | – | – | – |
| merges accepted | 6 | 4 | 2 | 4 | 7 |

**Cause and effect.** Fragmentation falls 26–42% in every species, and miniprot
wins 28–37% of loci everywhere. That reproducibility is the point: if EviAnn
fragmentation were random noise, miniprot would not systematically recover a
better model at a third of loci in five independent annotations. It is a
consistent property of how the annotator handles tandem arrays.

Only 23 merges were accepted across all five species. **Almost all of the
improvement is replacing one fragmented model with a better model at the same
locus, not merging two genes into one** — so repair improves completeness without
materially deflating gene counts. Four guards must pass before a merge; the ones
that fired were g3 (a spanned gene has independent transcript support) and g4 (the
merged model's exons do not cover a fragment).

**Asymmetry to be aware of:** ORF models exist only for hmaj, where they won 22
loci. Running getorf on the other four would let their coverage improve slightly
further. The effect is bounded at ~3% of loci and would *raise* the others toward
hmaj, so it cannot explain hmaj's excess.

## Disruption evidence

| | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| in-frame stops (miniprot) | 437 | 568 | 470 | 613 | 552 |
| frameshifts | 845 | 804 | 851 | 1012 | 1040 |
| pseudogene genes | 141 | 114 | 131 | 142 | 139 |

**Pseudogene counts are the most stable thing in this analysis: 114–142 genes,
a 1.25× spread, against a 1.58× spread in total genes.** Every one of these
assemblies carries roughly 130 decayed V2R copies. That is what a tandem array
looks like at equilibrium, and it is consistent with these being closely related
species sharing an ancestral array structure.

Note hmaj has the *fewest* in-frame stops (437) despite the most genes. If hmaj's
extra genes were assembly noise you might expect them to be full of disruptions;
they are not, which argues the extra copies are sequence-intact — see below for
what they actually are.

## gene_status — all six rules, all five species

Genes (not transcripts). Rule A is what is written to `gene_status`.

### functional

| rule | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| **A** intact + 7TM + no disruption | **275** | **153** | **165** | **180** | **168** |
| B intact + no disruption | 281 | 155 | 166 | 183 | 169 |
| C intact\|partial + 7TM + no disruption | 290 | 165 | 170 | 193 | 181 |
| D intact\|partial + no disruption | 299 | 167 | 171 | 196 | 184 |
| E 7TM + no disruption | 400 | 238 | 227 | 283 | 240 |
| F intact + 7TM, **ignores disruption** | 388 | 246 | 271 | 284 | 285 |

### pseudogene

| rule | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| A–E (identical) | 141 | 114 | 131 | 142 | 139 |
| F | 33 | 24 | 28 | 42 | 25 |

### partial

| rule | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| A, F | 247 | 170 | 128 | 194 | 123 |
| B | 239 | 164 | 121 | 188 | 117 |
| C | 232 | 158 | 123 | 181 | 110 |
| D | 221 | 152 | 116 | 175 | 102 |
| E | 123 | 85 | 67 | 91 | 51 |

**Reading this table.**

- **`pseudogene` is identical under A, B, C, D and E** in all five species. It
  depends only on whether the reading frame is disrupted, never on coverage or
  7TM. It is the one term that is a biological claim rather than a quality
  threshold, and it is the most robust number here.
- **A→D moves 9–24 genes** (letting 0.60–0.80 coverage count as functional).
  Small, defensible either way. A is the conservative choice.
- **A→E moves 59–125 genes.** E drops the coverage requirement entirely, so any
  7TM-bearing fragment becomes "functional". That is why E's `partial` collapses
  to 51–123: it is not finding more genes, it is relabelling fragments.
- **F is a different measurement, not a variant.** It ignores disruption, so
  pseudogene collapses from ~130 to ~30 in every species. Its functional count
  (271–388) is inflated by exactly the genes A calls pseudogenes. **Do not compare
  an F count to an A count.**

The ordering of species is stable across A, B, C and D — hmaj > hcurw > horn ≈
hcure > hcy — so the choice of rule does not change which snake has the largest
repertoire. Under E and F the ordering shuffles, which is another reason not to
use them.

## High-confidence set

7TM present + ≥6 exons + 700–1100 aa protein.

| | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| transcripts | 246 | 157 | 188 | 170 | 177 |
| genes | 237 | 148 | 182 | 161 | 174 |

Note this ordering differs from rule A: **hcure ranks second here but fourth by
functional count.** hcure has proportionally more models that pass all three
structural filters — consistent with its higher `complete` evidence rate (129
transcripts with both transcript and protein support, the highest of the five).

---

# INVESTIGATION — why does hmaj have 1.6× the genes of hcure?

Three candidate explanations. Two are testable with what is already computed.

## Hypothesis 1: hmaj really has more V2Rs

Possible — *H. major* is a distinct species and V2R arrays expand and contract
fast. But it should be believed only after the artefacts are ruled out.

## Hypothesis 2: uncollapsed heterozygosity — **supported**

hmaj has the **largest assembly (2150 Mb, 90–230 Mb more than the others) and the
most sequences (1111)**. Both are classic signatures of a diploid assembly that
retained haplotypes instead of collapsing them. If a V2R array is represented on
both a chromosome and an unplaced scaffold, every gene in it is counted twice.

**Test:** cluster each species' V2R proteins at 99% identity and count scaffold
genes that fall in a cluster with a chromosome gene.

| | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| genes on scaffolds | 198 | 3 | **240** | 4 | 7 |
| of those, ≥99% identical to a chromosome copy | **63** | 1 | **6** | 1 | 3 |
| as % of scaffold genes | **32%** | – | **2.5%** | – | – |

**63 of hmaj's 198 scaffold genes are near-identical duplicates of a chromosome
copy.** No other species shows this — including hcure, which has *more* scaffold
genes (240) but only 6 duplicates. hcure's unplaced V2Rs are genuinely distinct
sequences that simply failed to scaffold; hmaj's are substantially redundant.

**Corrected hmaj estimate: 661 − 63 ≈ 598 genes.** Still the largest, but the gap
to hcurw narrows from 148 to 85 genes.

The 99%-non-redundant counts make the same point more directly: **485 / 443 / 383
/ 379 / 375**. The spread drops from 1.58× to 1.29×, and hcure — last by raw gene
count — moves to essentially tied with hcy and horn.

## Hypothesis 3: annotation depth — **not supported**

If hmaj simply had a more thorough annotation, its total mRNA count would be
higher. It is 52,501 against 51,220–51,588 — a 2% difference against a 58%
difference in V2R genes. Annotation depth does not explain it.

But note the pHMM hit counts: 939 for hmaj versus 607–729. hmaj's proteome
contains 30–55% more V2R-like proteins from a similar amount of genomic V2R
signal (miniprot alignments are flat at 2674–2993). **That is exactly what
retained haplotypes would produce** — the same locus annotated twice, giving two
proteins, two pHMM hits, and two genes. Hypotheses 2 and 3 are the same
observation seen from two directions.

## Is fragmentation alone the cause? — **no**

Panel C of the figure plots assembly sequence count against V2R genes. There is no
clean relationship: hcurw has the *fewest* sequences (31) and the second *most*
genes (513), while hcure has 361 sequences and the fewest (418). Fragmentation is
not a general predictor. **hmaj is a specific case — a large, fragmented,
haplotype-retaining assembly — not an instance of a trend.**

## What this means for the headline

| Estimate | hmaj | hcy | hcure | hcurw | horn | spread |
|---|---|---|---|---|---|---|
| raw genes | 661 | 433 | 418 | 513 | 425 | 1.58× |
| minus scaffold duplicates | 598 | 432 | 412 | 512 | 422 | 1.45× |
| 99% non-redundant | 485 | 379 | 375 | 443 | 383 | **1.29×** |
| functional (rule A) | 275 | 153 | 165 | 180 | 168 | 1.80× |
| pseudogene | 141 | 114 | 131 | 142 | 139 | **1.25×** |

**Quote the non-redundant count for cross-species comparison, and the raw count
only with the assembly caveat attached.** The functional count has the widest
spread of all and is the most sensitive to assembly quality, because a duplicated
intact gene is counted twice while a duplicated pseudogene often is not (the two
copies decay differently and one may fail the intact-coverage test).

## Loose ends worth closing

1. **Run getorf on the four non-hmaj species.** ORF models won 22 hmaj loci; the
   others have no ORF track, so their coverage repair is slightly handicapped.
   Bounded effect, but it removes an asymmetry.
2. **Add a GRIK/iGluR profile to the panel.** It would remove the 5–7
   architecture-inconsistent calls per species at source instead of by filter.
3. **The hmaj duplicate scaffolds should be purged** (purge_dups or equivalent)
   and the annotation re-run, if hmaj's repertoire size is going to be a headline
   claim.
4. **Phylogeny** — now the natural next step, with all five curated. A single tree
   of the intact sets with decoy outgroups would assign subfamilies, resolve the
   ambiguous calls, and show whether hmaj's extra copies form species-specific
   clades (real expansion) or pair up with existing genes (residual duplication).

---

# PIPELINE AND FILES

## Directory

```
curation/
  run_curation.sh          evidence tracks -> loci -> panel scan -> 7TM -> tables
  make_tables.py           scans -> assignments, catalogue, master table
  plot_genome.py           master table -> per-species genome map
  panel/                   the shared Class C pHMM panel (built once)
  eggnog/                  DB fetch, emapper runner, per-species output
  processing/              00-05 scripts, reports, work dirs
  results/<species>/       one directory per species, identical layout
```

## Scripts, in order

| Script | Does |
|---|---|
| `processing/00_prepare_species.sh <sp>` | miniprot + proteome pHMM for a new species |
| `run_curation.sh <sp>` | the pipeline |
| `processing/01_positive_control.sh <sp>` | 00 + pipeline + score against the hcy benchmark |
| `processing/02_model_repair.py <sp>` | EviAnn vs miniprot vs ORF, merge guards |
| `processing/03_pseudogene_calls.py <sp...>` | `gene_status`, rule comparison |
| `processing/04_collate.py <sp...>` | cross-species tables |
| `processing/05_compare_plot.py` | the comparison figure |

Full five-species run from scratch:

```bash
for sp in hmaj hcy hcure hcurw horn; do
    bash processing/00_prepare_species.sh $sp
    bash run_curation.sh $sp
    python3 processing/02_model_repair.py $sp
done
python3 processing/03_pseudogene_calls.py hmaj hcy hcure hcurw horn
python3 processing/04_collate.py hmaj hcy hcure hcurw horn
python3 processing/05_compare_plot.py
```

miniprot is ~70 s per genome; the whole thing is under 30 minutes for five species.

## Outputs

`results/<sp>/v2r_master.tsv` — 68 columns, every candidate locus × transcript.
Everything else is a view of it:

```bash
is_v2r == yes                                   # the repertoire
gene_status == functional                       # rule A
repaired_coverage_class == intact               # complete models
eggnog_arch_consistent == no                    # the GRIK/PRRG1 contamination
v2r_evidence == panel                           # drop EviAnn-rescued fragments
has_7tm==yes AND n_exons>=6 AND 700<=protein_len_aa<=1100   # high confidence
```

Cross-species: `processing/04_cross_species.tsv` / `.md`,
`processing/05_cross_species.png`.

## Caveats carried forward

- `gene_status` is per transcript. A gene whose isoforms disagree appears in more
  than one bucket, so the per-gene columns sum to slightly more than the gene
  total (hmaj: 663 vs 661).
- Only hmaj has floor and control tracks. Floor tracks cannot change locus counts;
  the genome2or control can demote a T2, and only hmaj has it (1 locus affected).
- 23–35% of models are still `partial` — an assembly and annotation limit, not
  biology. It is the honest bucket for models that cannot be judged.
- Transcript (RNA-seq) support is low everywhere: `complete` evidence covers
  64–129 transcripts per species. These are mostly homology projections, expected
  for a vomeronasal-organ-restricted family.
- The positive control (hcy vs the published funannotate proteome) is **+6.1% at
  gene level**. Since this assembly was rebuilt from the same reads as the
  published one, that is best framed as *recovering 6.1% more V2R genes from the
  same data*, not as agreement between two independent genomes.
