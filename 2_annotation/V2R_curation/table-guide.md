# Table guide

Every file in `results/<species>/`, and every column in it. The same 20 files
exist for all five genomes — `hmaj`, `hcy`, `hcure`, `hcurw`, `horn` — with
identical columns; only the row counts differ.

Row counts quoted below are hmaj. The cross-species table at the end gives the
rest.

**If you only read one thing:** `v2r_master.tsv` is the table. Everything else is
either a view of it, an input to it, or a summary of it.

## Reading the source column

Every table below has a **source** column naming where each field's value comes
from. Six kinds of source appear:

| source | means |
|---|---|
| **EviAnn GFF** | copied straight out of the annotation's `Note=`/`ID=`/attribute fields — this pipeline never edits it |
| **EviAnn protein FASTA** | from the `>...Name:"..."` header of the protein file EviAnn also produced |
| **miniprot** | from `miniprot.v2r.gff`, a spliced protein-to-genome alignment run *before* this pipeline (`processing/00_prepare_species.sh`) |
| **pHMM scan** | from `hmmsearch`/`hmmscan` output — either the V2R profile vs the proteome, or the Class C panel, or the Pfam 7tm_3 scan |
| **eggNOG** | from `eggnog/run_emapper.sh`, run after the panel call, on the retained V2R proteins only |
| **derived** | computed by this pipeline's own scripts from the above — never present in any upstream file |

Script names are `01_run_curation.sh` (bash, tracks → loci → scans),
`make_tables.py` (assignment + all tables), `02_model_repair.py`,
`03_pseudogene_calls.py`. Section numbers refer to the `# N.` comment banners
inside each script.

---

# 1. `v2r_master.tsv` — the table

One row per **candidate locus × overlapping transcript**. Every candidate locus
appears whether or not it turned out to be V2R, so this is the full search space,
not the answer. 68 columns, 1,294 rows.

A locus with three overlapping transcripts contributes three rows sharing the
same locus columns. A locus with no annotated transcript contributes one row with
`NA` in every transcript column (373 such rows in hmaj).

## Locus block — columns 1–19

These describe the genomic region and repeat identically across every row of that
locus.

| # | column | source | meaning |
|---|---|---|---|
| 1 | `locus_id` | **derived** — `01_run_curation.sh` §3 | stable ID, issued in genome order: `HMAJV2R_ch2_000252`. `UNPLACED` for EviAnn-labelled V2Rs that no locus covers. Not an EviAnn ID — EviAnn has no concept of a locus, only genes/transcripts (`<SPECIES>V2R_<chrom>_<counter>`) |
| 2 | `chrom` | **derived**, from `bedtools merge` of the tracks below | sequence name |
| 3 | `seq_type` | **derived** — `make_tables.py::seq_type()` | `chromosome` if the name starts with `ch`, else `scaffold` |
| 4–5 | `locus_start`, `locus_end` | **derived** | 0-based half-open interval, the merged extent of the primary tracks |
| 6 | `locus_len_bp` | **derived** | end − start |
| 7 | `locus_tier` | **derived** — `make_tables.py::catalogue()` | `T2` confidently V2R / `T1` supported but no confident call / `T0` confidently another Class C family or unsupported |
| 8 | `locus_family` | **derived**, from the pHMM Class C panel scan | the family of the locus's strongest confident call, `NA` if none |
| 9 | `locus_delta_bits` | **pHMM scan** (Class C panel `hmmscan`) | that call's margin over the runner-up family |
| 10 | `n_primary_tracks` | **derived** — `01_run_curation.sh` §4 | how many of the two primary engines (miniprot, pHMM-vs-proteome) support it (1 or 2). Floor tracks deliberately excluded |
| 11 | `tracks` | **derived** | every track supporting the locus, comma-separated, floor tracks included |
| 12 | `control_overlap` | **derived**, from the genome2or control track (external tool, not part of this pipeline) | `genome2or` if the olfactory-receptor control overlaps — the locus was demoted T2→T1 |
| 13–14 | `n_transcripts_in_locus`, `n_genes_in_locus` | **derived**, from `eviann_mrna.bed` (EviAnn GFF) | EviAnn features overlapping. >1 gene means the locus spans a tandem array |
| 15 | `miniprot_alns` | **miniprot** | miniprot alignments landing here |
| 16 | `frameshifts` | **miniprot** | **frameshifts summed over those alignments** — miniprot's own `Frameshift=` GFF attribute |
| 17 | `inframe_stops` | **miniprot** | **in-frame stop codons summed over those alignments** — miniprot's own `StopCodon=` GFF attribute |
| 18 | `miniprot_max_identity` | **miniprot** | best alignment identity, 0–1, from miniprot's `Identity=` attribute |
| 19 | `locus_flags` | **derived** — `make_tables.py::catalogue()` | `unclassified`, `assigned_GRM`, `control_overlap:genome2or`, `no_candidate_locus` |

> **Columns 15–18 are per *locus*, not per transcript.** A transcript inherits the
> disruption counts of its whole locus. This is why a clean-looking transcript can
> come out `pseudogene` — see the worked example below.

## Transcript block — columns 20–41

| # | column | source | meaning |
|---|---|---|---|
| 20 | `transcript_id` | **EviAnn GFF** | mRNA `ID=` attribute, verbatim |
| 21 | `gene_id` | **EviAnn GFF** | mRNA `Parent=` attribute, verbatim. **Count genes on this, never on `locus_id`** |
| 22–25 | `tx_start`, `tx_end`, `strand`, `tx_span_bp` | **EviAnn GFF** (22–24 verbatim; 25 derived = end − start) | the transcript's own coordinates |
| 26 | `n_exons` | **EviAnn GFF** | the `Num_exons=` attribute, verbatim |
| 27 | `protein_len_aa` | **EviAnn protein FASTA** | length of the sequence in `candidates.faa`, trailing `*` stop excluded |
| 28 | `v2r_model_cov` | **derived**, from the **pHMM scan** (V2R profile vs proteome) | fraction of the 823-column V2R profile this protein covers. **The per-transcript completeness measure** |
| 29 | `coverage_class` | **derived** from column 28 | `intact` ≥0.80 / `partial` 0.60–0.80 / `fragmented` <0.60 / `not_assessed` (no V2R hit) |
| 30–33 | `has_7tm`, `n_7tm_domains`, `tm7_len_aa`, `tm7_evalue` | **pHMM scan** (Pfam 7tm_3 vs `candidates.faa`) | the Class C seven-transmembrane domain |
| 34 | `eviann_evidence` | **EviAnn GFF** | `Evidence=` attribute, verbatim: `complete` (transcript **and** protein evidence) / `protein_only` / `transcript_only` |
| 35–36 | `start_codon`, `stop_codon` | **EviAnn GFF** | `StartCodon=`/`StopCodon=` attributes. Always present — EviAnn only emits complete ORFs |
| 37 | `gene_biotype` | **EviAnn GFF** | `gene_biotype=` attribute: `protein_coding` or `processed_pseudogene` (EviAnn's own call, 31 rows) |
| 38 | `eviann_label` | **EviAnn GFF** | the `Note=` attribute, verbatim — EviAnn's homology label |
| 39 | `eviann_label_status` | **derived** from column 38 | `annotated` / `unknown` (the literal string "function unknown") / `none` |
| 40 | `found_by_phmm_prot` | **derived**, from the **pHMM scan** target list | did the proteome pHMM hit this transcript directly |
| 41 | `tx_credited_to_locus` | **derived** — `01_run_curation.sh` §5 | `no` when the transcript straddles two loci and is credited to the other one. **Filter `!= no` for one row per transcript** |

## Call block — columns 42–51

| # | column | source | meaning |
|---|---|---|---|
| 42 | `classC_family` | **pHMM scan** (Class C panel `hmmscan` on `candidates.faa`) | the call: a family name, or `ambiguous` / `single_family_weak` / `below_threshold` / `not_scanned` |
| 43 | `best_bits` | **pHMM scan** | winning family's bitscore |
| 44 | `delta_bits` | **derived** from the pHMM scan — `make_tables.py::assign()` | margin over the runner-up family. **The confidence measure** |
| 45 | `evalue` | **pHMM scan** | winning hit's E-value |
| 46 | `confident` | **derived** from 44 | 1 when `delta_bits ≥ 20` (or ≥100 bits with no competitor) |
| 47 | `is_v2r` | **derived** — `make_tables.py::transcript_block()` | **the retention decision** — yes for 736 hmaj rows |
| 48 | `v2r_evidence` | **derived** | `panel` (731) or `eviann_rescue` (5: EviAnn-labelled, panel agrees, outside every locus) |
| 49 | `excluded_reason` | **derived** | `not_a_confident_V2R_call` where `is_v2r=no` |
| 50 | `eviann_product` | **EviAnn protein FASTA** | product name parsed from the `Name:"..."` field of the header |
| 51 | `product_agreement` | **derived**, comparing column 42 against column 50 | `agree` / `disagree` / `no_annotation` — our call vs EviAnn's label |

## eggNOG block — columns 52–63

All from **eggNOG** (`eggnog/run_emapper.sh`, eggNOG-mapper 5.0.2), run on the
retained V2R proteins only, *after* the panel call. `-` where eggNOG had
nothing; `NA` on rows never sent to it (non-V2R rows).

`eggnog_eggnog_ogs`, `eggnog_cog_category`, `eggnog_description`,
`eggnog_preferred_name`, `eggnog_gos`, `eggnog_ec`, `eggnog_kegg_ko`,
`eggnog_kegg_pathway`, `eggnog_kegg_module`, `eggnog_brite`, `eggnog_pfams`.

| # | column | source | meaning |
|---|---|---|---|
| 63 | `eggnog_arch_consistent` | **derived**, from `eggnog_pfams` (column 62) — `make_tables.py::transcript_block()` | `yes` when the Pfam architecture contains both `7tm_3` and `ANF_receptor`. **`no` = 7 rows, all GRIK1/GRIK2/PRRG1** — the panel's known blind spot |

## Repair and status — columns 64–68

Added by `02_model_repair.py` (64–67) and `03_pseudogene_calls.py` (68) — these
two scripts run *after* `01_run_curation.sh`/`make_tables.py` and rewrite
`v2r_master.tsv` with extra columns appended.

| # | column | source | meaning |
|---|---|---|---|
| 64 | `best_model_source` | **derived**, comparing **pHMM-scan** coverage of the EviAnn protein against **miniprot**'s own alignment and (where run) a getorf ORF, all rescored with the same V2R profile | which annotation gave the most complete model at this locus: `eviann` / `miniprot` / `orf` |
| 65 | `repaired_model_cov` | **derived** | that model's V2R-profile coverage |
| 66 | `repaired_coverage_class` | **derived** from 65 | `intact` / `partial` / `fragmented`, recomputed on the winning model |
| 67 | `merge_status` | **derived**, from miniprot alignment spans + EviAnn gene spans/evidence | `none`, `accepted`, or `rejected:<guard>` where a miniprot model claimed two EviAnn genes are one |
| 68 | `gene_status` | **derived**, from columns 66 (`repaired_coverage_class`) and 16–17 (`frameshifts`/`inframe_stops`) — `03_pseudogene_calls.py` using rule B from `v2rlib.py` (final spec: `intact` + no disruption, no 7TM requirement) | `functional` / `pseudogene` / `partial` |

> **Columns 64–67 are also per *locus*.** A short isoform of a locus whose best
> model is full length is labelled `intact`. Filtering
> `repaired_coverage_class == intact` selects **loci with an intact model**, not
> **transcripts that are intact**. For a strictly per-transcript measure use
> `v2r_model_cov` (column 28).

## A worked row

`LOC_00010372-mRNA-1` at `HMAJV2R_ch2_000252`:

```
coverage_class = partial   (v2r_model_cov 0.773, 623 aa)
has_7tm        = yes
frameshifts    = 3         <- locus-level, from 3 miniprot alignments
gene_status    = pseudogene
```

The transcript has a clean ORF and a 7TM domain. It is called `pseudogene`
because miniprot had to align through three frameshifts *somewhere in this
locus*. That is the intended behaviour — the disruption is genomic evidence about
the region, and EviAnn cannot report it because it only emits intact ORFs — but
it means `gene_status` is a statement about the locus, not proof this particular
transcript is broken.

---

# 2. Views of the master table

Same 68 columns, filtered.

| file | rows | filter |
|---|---|---|
| `v2r_genes.tsv` | 736 | `is_v2r == yes` — **the repertoire** |
| `classC_all_genes.tsv` | 761 | `confident == 1` — every confident Class C call, so it includes the 25 rejects (GRM, GABBR, CASR, GPRC6A). Use it to see what was excluded and why |

---

# 3. `v2r_catalogue.tsv` — one row per locus

1,073 rows, 15 columns. The tiering view, without transcripts. Written by
`make_tables.py::catalogue()`, same script/section as the master's locus block —
every column here has the identical provenance as its master-table counterpart.

| # | column | source | meaning |
|---|---|---|---|
| 1–6 | `locus_id`, `chrom`, `start`, `end`, `strand`, `length` | **derived** (see master cols 1–6) | the locus. `strand` is always `.` — strand is not resolved per locus |
| 7 | `tier` | **derived** | T2 607 / T1 453 / T0 13 |
| 8 | `family` | **pHMM scan** | strongest confident call, or `NA` |
| 9 | `delta_bits` | **pHMM scan**, derived | its margin |
| 10 | `n_exons` | **EviAnn GFF** (`Num_exons=` of the calling transcript) | exon count **of the calling transcript**, not of the locus |
| 11 | `n_tracks` | **derived** | primary tracks only (same as `n_primary_tracks` in the master) |
| 12 | `tracks` | **derived** | all supporting tracks |
| 13 | `eviann_gene_id` | **EviAnn GFF** | every EviAnn transcript and gene overlapping, comma-separated |
| 14 | `control_overlap` | **derived**, from genome2or | `genome2or` or empty |
| 15 | `flags` | **derived** | `unclassified` (452), `assigned_GRM` (8), `assigned_CASR` (2), `assigned_GABBR` (2), `control_overlap:genome2or` (1) |

---

# 4. `classC_assign.tsv` — the raw call

822 rows, one per scanned protein, straight from the hmmscan. No locus context.
Written by `make_tables.py::write_assignments()`.

| # | column | source | meaning |
|---|---|---|---|
| 1 | `query_id` | **pHMM scan** | transcript |
| 2–3 | `best_family`, `best_bits` | **pHMM scan** | winning family and its bitscore |
| 4–5 | `runner_up_family`, `runner_up_bits` | **pHMM scan** | second place — `NONE`/0 when only one family hit |
| 6 | `delta_bits` | **derived** | 2 − 4. When there is no runner-up this holds the raw bitscore instead |
| 7 | `evalue` | **pHMM scan** | winning hit |
| 8 | `confident` | **derived** | 1/0 |
| 9 | `assignment` | **derived** | the family name, or `ambiguous` (13) / `single_family_weak` (48) / `below_threshold` |

---

# 5. `chrom_summary.tsv` — per sequence

106 rows for hmaj: one per sequence carrying at least one candidate locus.
Written by `make_tables.py::chrom_summary()` — every value is a count over
`v2r_master.tsv` rows, so all of it is **derived**; nothing here is a raw field
from any upstream source.

`sequence`, `seq_type`, `n_loci`, `T2`, `T1`, `T0`, `n_v2r_transcripts`,
`n_v2r_genes`.

Its row count is itself informative — 106 sequences for hmaj and 126 for hcure
against 9–12 for the well-scaffolded genomes.

---

# 6. Intermediates

Inputs to the tables above, not results. All produced by `01_run_curation.sh`;
regenerated every run.

| file | rows | script section | columns / source |
|---|---|---|---|
| `loci.bed` | 1,073 | §3 | BED6: chrom, start, end, `locus_id`, 0, `.` — **derived** from `bedtools merge` of miniprot + pHMM-vs-proteome tracks |
| `locus_tracks.tsv` | 1,073 | §4 | `locus_id`, `n_primary_tracks`, `tracks`, `control_overlap` — all **derived** |
| `locus_track_pairs.tsv` | 4,357 | §4 | long form: `locus_id`, `track` — one row per pair, no header, **derived** |
| `locus_transcripts.tsv` | 916 | §5 | `locus_id`, `transcript_id`, `gene_id` (**EviAnn GFF**), `phmm_hit` (**derived** from the pHMM scan's target list). Every overlap, including transcripts credited elsewhere |
| `mrna_locus.tsv` | 908 | §5 | `transcript_id`, `locus_id` — **derived**; each transcript credited to exactly one locus (pHMM-supported wins, then lowest ID) |
| `locus_miniprot.tsv` | 695 | §5b | `locus_id`, `miniprot_alns`, `miniprot_frameshifts`, `miniprot_stops`, `miniprot_max_identity` — **miniprot**'s own GFF attributes, aggregated per locus. Source of master columns 15–18 |
| `v2r_domains.tsv` | 939 | §5b | `transcript_id`, `model_len` (823), `hmm_from`, `hmm_to` — **pHMM scan** (V2R profile vs proteome), one row per domain. Unioned into `v2r_model_cov` |
| `candidate_ids.txt` | 913 | §5 | transcripts sent to the panel, one per line — **derived** union of locus-overlapping transcripts and `eviann_v2r_labelled.txt` |
| `eviann_v2r_labelled.txt` | 624 | §1 | transcripts EviAnn itself called V2R (its `Note=` mentions vomeronasal/vmn2r/v2r), before this pipeline looked — **EviAnn GFF**, filtered |
| `locus_control.txt` | 1 | §4 | loci the olfactory-receptor control (genome2or) flags — **derived** |

---

# 7. Whole-annotation tables

**Not V2R-specific.** One row per mRNA in the entire annotation — the lookup
tables the pipeline joins against. All three come straight out of the **EviAnn
GFF** in a single pass, `01_run_curation.sh` §1 (one awk pass over the ~125 MB
file, since it's the same input read four ways).

| file | rows | source | columns |
|---|---|---|---|
| `eviann_mrna.bed` | 52,501 | **EviAnn GFF**, verbatim coordinates | chrom, start, end, `transcript_id`, `.`, strand, `gene_id` (7-column BED) |
| `mrna_attrs.tsv` | 52,501 | **EviAnn GFF**, verbatim attributes | `transcript_id`, `evidence`, `start_codon`, `stop_codon`, `num_exons`, `gene_biotype`, `note` |
| `exon_counts.tsv` | 52,500 | **EviAnn GFF**, `Num_exons=` attribute | `Parent=<transcript_id>`, exon count. **No header** — the `Parent=` prefix is legacy from the original input format |

`mrna_attrs.tsv` is the largest file in the directory at 6 MB.

---

# 8. Non-TSV outputs

`candidates.faa` (proteins scanned), `classC.tblout` (hmmscan output),
`pfam_7tm.domtbl` / `pfam_7tm.hmm` (7TM scan), `classC.err` / `pfam_7tm.err`
(HMMER stderr), `summary.txt` (parameters + headline numbers),
`v2r_genome_map.png` / `.svg`, `tracks/` (one BED per evidence track).

---

# Row counts across species

| file | hmaj | hcy | hcure | hcurw | horn |
|---|---|---|---|---|---|
| `v2r_master.tsv` | 1,294 | 924 | 870 | 1,084 | 864 |
| `v2r_genes.tsv` | **736** | **493** | **499** | **574** | **482** |
| `classC_all_genes.tsv` | 761 | 518 | 526 | 599 | 508 |
| `v2r_catalogue.tsv` | 1,073 | 715 | 673 | 895 | 668 |
| `classC_assign.tsv` | 822 | 552 | 558 | 636 | 528 |
| `chrom_summary.tsv` | 106 | 12 | 126 | 9 | 11 |
| `locus_transcripts.tsv` | 916 | 652 | 645 | 738 | 638 |
| `mrna_locus.tsv` | 908 | 642 | 617 | 726 | 623 |
| `locus_miniprot.tsv` | 695 | 528 | 527 | 652 | 492 |
| `v2r_domains.tsv` | 939 | 633 | 630 | 729 | 607 |
| `mrna_attrs.tsv` | 52,501 | 51,220 | 51,588 | 51,444 | 51,521 |

Only hmaj has floor tracks (ORF, 6-frame, tBLASTn) and the genome2or control, so
its `tracks` and `control_overlap` columns are richer. That cannot change locus
counts by design — floor tracks annotate a locus a primary engine already found.

---

# Cross-species tables

In `processing/`, not per species.

| file | what |
|---|---|
| `04_cross_species.tsv` | one row per species, every headline number; `.md` is the same as tables |
| `03_status_table.md` | functional/pseudogene/partial under all six rules, before and after repair |
| `02_repair_report_<sp>.md` | coverage before/after, winning model source, merges accepted |
| `02_merge_review_<sp>.tsv` | one row per proposed merge: `locus_id`, `miniprot_model`, `query`, `query_frac`, `n_genes_spanned`, `genes`, `eviann_cov`, `miniprot_cov`, `merge_status` |
| `01_control_report.md` | hcy scored against its independent benchmark |

---

# Common filters

```bash
M=results/hmaj/v2r_master.tsv

# the repertoire
awk -F'\t' 'NR==1||$47=="yes"' $M

# intact, full-length, undisrupted (no 7TM requirement)
awk -F'\t' 'NR==1||($47=="yes"&&$68=="functional")' $M

# high confidence, in-depth analysis set: functional + >=6 exons + 700-1100 aa
# (locus_tier==T2 is implied by gene_status==functional -- see README Method #8)
awk -F'\t' 'NR==1||($47=="yes"&&$68=="functional"&&$26>=6&&$27>=700&&$27<=1100)' $M

# drop the GRIK/PRRG1 contamination
awk -F'\t' 'NR==1||($47=="yes"&&$63!="no")' $M

# one row per transcript (drop the straddling duplicates)
awk -F'\t' 'NR==1||$41!="no"' $M

# count genes, never loci
awk -F'\t' 'NR>1&&$47=="yes"&&$21!="NA"{g[$21]}END{print length(g)}' $M
```

Column numbers are stable across species. Prefer looking the name up from the
header if you are scripting against this long-term.
