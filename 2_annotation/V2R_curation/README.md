# V2R curation — five *Hydrophis* sea snake genomes

Curates vomeronasal type-2 receptor (V2R) genes from EviAnn annotations by
reconciling every search engine run against each assembly, discriminating V2R
from the other Class C GPCR families with a competing profile-HMM panel,
repairing fragmented gene models, calling pseudogenes, and annotating the
result functionally.

---

# TLDR


| stage                                                                | hmaj            | hcy             | hcure           | hcurw           | horn            |
| -------------------------------------------------------------------- | --------------- | --------------- | --------------- | --------------- | --------------- |
| EviAnn original V2R-labelled (`Similar to Vmn2r*`)                   | 624             | 402             | 448             | 479             | 434             |
| T2 (loci, tier fixed pre-repair)                                     | 607             | 393             | 400             | 483             | 390             |
| T2 genes, post model repair                                          | 655             | 430             | 416             | 509             | 423             |
| partial (tx / genes)                                                 | 255 / 239       | 176 / 164       | 137 / 121       | 197 / 188       | 127 / 117       |
| pseudogene (tx / genes)                                              | 180 / 141       | 142 / 114       | 166 / 131       | 166 / 142       | 171 / 139       |
| functional (tx / genes)                                              | 301 / 281       | 175 / 155       | 196 / 166       | 211 / 183       | 184 / 169       |
| **high confidence** (tx / genes / loci) – intact + T2 + functional   | 301 / 281 / 268 | 175 / 155 / 143 | 196 / 166 / 163 | 211 / 183 / 180 | 184 / 169 / 165 |
| **strict** = + ≥6 exons + 700–1100 aa (tx / genes / loci)            | 181 / 176 / 172 | 112 / 106 / 103 | 124 / 122 / 121 | 117 / 110 / 110 | 118 / 118 / 118 |
| eggNOG annotated / canonical domain (`7tm_3`+`ANF_receptor`+`NCD3G`) | 736 / 729       | 493 / 487       | 499 / 494       | 574 / 568       | 482 / 476       |
| genes on chromosomes                                                 | 463             | 430             | 178             | 509             | 418             |
| genes on scaffolds (unplaced)                                        | 198 (30.0%)     | 3 (0.7%)        | 240 (57.4%)     | 4 (0.8%)        | 7 (1.6%)        |
| **final V2R protein count** (transcripts, `is_v2r`=yes)              | **736**         | **493**         | **499**         | **574**         | **482**         |
| **final V2R gene count**                                             | **661**         | **433**         | **418**         | **513**         | **425**         |




## Placement per chromosome (genes)


| chrom                | hmaj    | hcy     | hcure   | hcurw   | horn    |
| -------------------- | ------- | ------- | ------- | ------- | ------- |
| ch1                  | 40      | 25      | 34      | 56      | 35      |
| ch2                  | 191     | 188     | 75      | 181     | 205     |
| ch3                  | 4       | 0       | 2       | 0       | 0       |
| ch4                  | 112     | 111     | 1       | 102     | 56      |
| ch5                  | 4       | 21      | 1       | 4       | 1       |
| ch6                  | 0       | 0       | 0       | 0       | 0       |
| ch7                  | 2       | 2       | 4       | 2       | 2       |
| ch8                  | 0       | 1       | 2       | 0       | 0       |
| chZ                  | 110     | 82      | 59      | 164     | 119     |
| scaffolds (unplaced) | 198     | 3       | 240     | 4       | 7       |
| **total genes**      | **661** | **433** | **418** | **513** | **425** |


![[v2r position on genome]](./V2R_curation.png)


V2R curation workflow

---


| Script                              | Does                                                                  |
| ----------------------------------- | --------------------------------------------------------------------- |
| `v2rlib.py`                         | thresholds, shared readers, and the classification rules              |
| `01_run_curation.sh`                | evidence tracks → candidate loci → Class C panel → hmmscan → 7TM scan |
| `make_tables.py`                    | scans → assignments, catalogue, master table, summary                 |
| `plot_genome.py`                    | master table → genome map                                             |
| `processing/00_prepare_species.sh`  | the two primary evidence tracks for a new genome                      |
| `processing/01_positive_control.sh` | run on *H. cyanocinctus* and score against its benchmark              |
| `processing/02_model_repair.py`     | pick the most complete model per locus, guard any merge               |
| `processing/03_pseudogene_calls.py` | functional / pseudogene / partial                                     |
| `processing/04_collate.py`          | every species into one comparison                                     |
| `eggnog/`                           | database fetch and eggNOG-mapper run                                  |


Outputs land in `results/<species>/`; cross-species tables in `processing/`.

---



# METHOD

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



## Evidence tracks

Each search becomes a BED file and is given one of three roles:


| Track       | Role        | What it is                                                                 |
| ----------- | ----------- | -------------------------------------------------------------------------- |
| `miniprot`  | **primary** | spliced protein-to-genome alignment of 1,516 vomeronasal receptor proteins |
| `phmm_prot` | **primary** | V2R profile HMM vs the EviAnn proteome                                     |
| `genome2or` | control     | olfactory-receptor gene finder                                             |


**Only primary tracks can create a locus.** All five genomes now run on the same
two primary tracks (`miniprot` & `pHMM`)

## Candidate loci

```bash
cat tracks/miniprot.bed tracks/phmm_prot.bed \
  | sort -k1,1 -k2,2n | bedtools merge -d 1000 -i - | awk '$3-$2 >= 400'
```

- IDs are issued in genome order (`HMAJV2R_<chrom>_<n>`) and depend only on the input intervals, so they are stable across re-runs. 
- `merge -d 1000` is calibrated: on hmaj ... the locus count is 1,076 / 1,073 / 1,059 at gaps of 0 / 1,000 / 5,000 bp, and the **gene** count is identical at 0 and 1,000. 
- The oversized loci are therefore not a merge artefact —> they are chains of transitively overlapping miniprot alignments running the length of a tandem array.



## Class C panel

- each Class C family gets its own profile and assignment is argmax bitscore with a margin. 
- Decoys (376 non-V2R Class C proteins) are split to families by header keyword
- -> then `cd-hit -c 0.95` → `mafft --localpair` → `trimal -gt 0.5 -cons 60` → `hmmbuild`.


| family       | seeds    | after cd-hit | profile length |
| ------------ | -------- | ------------ | -------------- |
| GRM          | 161      | 25           | 888            |
| CASR         | 16       | 2            | 1089           |
| GPRC6A       | 8        | 4            | 926            |
| TAS1R        | 18       | 9            | 837            |
| GABBR        | 81       | 5            | 852            |
| GPRC5        | 19       | 7            | 400            |
| **V2R_OlfC** | prebuilt | —            | 823            |


- The panel is built once and **shared by all five species** (`panel/`), so bitscores are comparable across genomes. 
- The V2R profile is the pre-validated `v2r_full.hmm`

**What the engines were built from.** miniprot is not trained; its query is 1,516
vomeronasal receptor proteins (1,035 squamate UniProt, 43 rat, 41 *Xenopus*,
37 mouse, 36 squamate NCBI, 12 squamate Pfam). The V2R profile came from 581
curated seeds (332 mouse Vmn2r, 200 rat Vom2r, 49 *Varanus komodoensis*) expanded
and filtered to 1,683 validated sequences, reduced to 1,206 representatives,
aligned and trimmed to 823 columns. Both descend from vomeronasal collections —
independent in *method*, not in reference data, so their agreement is weaker
evidence than the EviAnn or eggNOG comparisons.

## Assignment

`hmmscan -E 1e-3` of every candidate protein against the whole panel:

```
delta_bits = bits(best family) - bits(runner-up family)
```


| Label                | Condition                                                 |
| -------------------- | --------------------------------------------------------- |
| `<family>`           | ≥2 families hit, `delta_bits >= 20` — confident           |
| `ambiguous`          | ≥2 families hit, `delta_bits < 20`                        |
| `single_family_weak` | only one family cleared the ceiling and scored < 100 bits |


- `single_family_weak` -> resolves if no runner-up + margin to test

---



## Tiers — per locus


| Tier   | Condition                                                       |
| ------ | --------------------------------------------------------------- |
| **T2** | the locus's best call is `V2R_OlfC` with `delta_bits >= 20`     |
| **T1** | a primary engine supports the locus, but no confident V2R call  |
| **T0** | confidently a *different* Class C family, or no primary support |


Overlap with the `genome2or` control demotes T2 → T1. There is deliberately no
structural gate inside the tier: exon count, coverage and 7TM presence are all
columns, so any integrity filter is a filter on the output.

---



## Model repair


| miniprot vs eviann coverage | intepretation          | model used |
| --------------------------- | ---------------------- | ---------- |
| ~equal                      | Eviann model = correct | `EviAnn`   |
| miniprot >> eviann          | Eviann fragmented gene | `miniprot` |
| miniprot << eviann          | miniprot mis-aligned   | `EviAnn`   |


- 45% of eviann v2r models cover < 60% of v2r profile
  - most = fragmented genes ... != short receptors
  - **@ the same loci: ==2 other model sources**==: `miniprot` alignments + `getorf ORFs`
    - for each locus = **most complete model wins** 
      - *pHMM scores proteins = doesnt build models*
- **miniprot model spanning >=2 eviann genes merges if:** 
  1. one query protein spans both fragments (≥80% of the query used);
  2. repaired coverage ≥ 0.90, not merely improved;
  3. no spanned gene has independent transcript support;
  4. the merged model's exons cover the fragments' spans.



## Completeness, before and after model repair


|                                        | hmaj           | hcy           | hcure         | hcurw         | horn          |
| -------------------------------------- | -------------- | ------------- | ------------- | ------------- | ------------- |
| intact, before → after                 | 348 → **459**  | 208 → **305** | 239 → **339** | 236 → **345** | 228 → **337** |
| partial                                | 56 → 27        | 34 → 18       | 30 → 9        | 46 → 22       | 38 → 18       |
| fragmented                             | 332 → **245**  | 251 → **167** | 230 → **149** | 292 → **203** | 216 → **125** |
| winning model: eviann / miniprot / orf | 500 / 209 / 22 | 327 / 163 / – | 313 / 184 / – | 364 / 206 / – | 299 / 181 / – |
| merges accepted                        | 6              | 4             | 2             | 4             | 7             |


- Repair moves 26–42% of fragmented models up a class in every genome. 
- miniprot gets 28–37% of loci —> consistent across five independent assemblies
  - which is what you would expect if the fragmentation is annotation breakage rather than short genes.

---



## Gene status

**Final specification:**

```
functional   repaired_coverage_class == intact
             AND inframe_stops == 0 AND frameshifts == 0
             (no 7TM requirement — see below)

pseudogene   inframe_stops > 0 OR frameshifts > 0

partial      everything else — truncated or too fragmented to judge.
             An assembly/annotation category, not a biological one.
```

- `inframe_stops` and `frameshifts` from **miniprot**, not from EviAnn. EviAnn cannot report a pseudogene: it only emits complete `ATG..stop` ORFs, so it declines to annotate the broken part rather than annotating it as broken. 
 `CHOSEN_RULE` **—** `intact`**, no disruption, no 7TM requirement => ****** `gene_status` 
- A locus's T2 tier (`classC_family == V2R_OlfC` at `delta_bits >= 20` over the runner-up family) and `gene_status == functional` turn out to be fully redundant in all five genomes: every `functional` transcript already sits at a T2 locus, in every species, with zero exceptions. 
- A locus's T2 tier (`classC_family == V2R_OlfC` at `delta_bits >= 20` over the runner-up family) and `gene_status == functional` turn out to be fully redundant in all five genomes: every `functional` transcript already sits at a T2 locus, in every species, with zero exceptions. 
  - => ***method validity***



## Functional / pseudogene / partial (final rule, after repair)

`repaired_coverage_class == intact AND inframe_stops == 0 AND frameshifts == 0`,
no 7TM requirement. Transcripts / genes.


|                | hmaj          | hcy           | hcure         | hcurw         | horn          |
| -------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| **functional** | 301 / **281** | 175 / **155** | 196 / **166** | 211 / **183** | 184 / **169** |
| **pseudogene** | 180 / **141** | 142 / **114** | 166 / **131** | 166 / **142** | 171 / **139** |
| **partial**    | 255 / 239     | 176 / 164     | 137 / 121     | 197 / 188     | 127 / 117     |




#### Disruption evidence


|                           | hmaj | hcy | hcure | hcurw | horn  |
| ------------------------- | ---- | --- | ----- | ----- | ----- |
| in-frame stops (miniprot) | 437  | 568 | 470   | 613   | 552   |
| frameshifts               | 845  | 804 | 851   | 1,012 | 1,040 |
| pseudogene genes          | 141  | 114 | 131   | 142   | 139   |


---



## Functional annotation (eggNOG)

- counts of retained V2R transcripts 
- **annotated** = V2R transcripts sent/returned by `eggNOG`
- **canonical** = How many have that expected Class C V2R Pfam set
- `ANF_receptor` = Venus Flytrap Domain


|                                        | hmaj | hcy | hcure | hcurw | horn |
| -------------------------------------- | ---- | --- | ----- | ----- | ---- |
| annotated                              | 736  | 493 | 499   | 574   | 482  |
| canonical `7tm_3, ANF_receptor, NCD3G` | 729  | 487 | 494   | 568   | 476  |


---



## Repertoire sizes


|                     | hmaj       | hcy        | hcure      | hcurw      | horn       |
| ------------------- | ---------- | ---------- | ---------- | ---------- | ---------- |
| assembly (Mb)       | 2,150      | 1,928      | 1,923      | 2,059      | 1,998      |
| annotated mRNAs     | 52,501     | 51,220     | 51,588     | 51,444     | 51,521     |
| candidate loci      | 1,073      | 715        | 673        | 895        | 668        |
| T2 / T1 / T0 loci   | 607/453/13 | 393/309/13 | 400/260/13 | 483/399/13 | 390/265/13 |
| **V2R transcripts** | **736**    | **493**    | **499**    | **574**    | **482**    |
| **V2R genes**       | **661**    | **433**    | **418**    | **513**    | **425**    |
| EviAnn-rescued      | 5          | 3          | 2          | 4          | 2          |
| 7TM present         | 592        | 388        | 411        | 459        | 416        |




## Placement in the assembly


| genes          | hmaj | hcy  | hcure   | hcurw | horn |
| -------------- | ---- | ---- | ------- | ----- | ---- |
| on chromosomes | 463  | 430  | **178** | 509   | 418  |
| on scaffolds   | 198  | 3    | **240** | 4     | 7    |
| % unplaced     | 30%  | 0.7% | **57%** | 0.8%  | 1.6% |


---



## High-confidence set

- **only** for specific/in-depth analysis (*e.g., for predictive protein folding*)
- V2Rs typically have 6 exons & aa length ~700-1100


|                                  | hmaj                | hcy                 | hcure               | hcurw               | horn                |
| -------------------------------- | ------------------- | ------------------- | ------------------- | ------------------- | ------------------- |
| **high confidence**              | 301 / 281 / 268     | 175 / 155 / 143     | 196 / 166 / 163     | 211 / 183 / 180     | 184 / 169 / 165     |
| + ≥6 exons + 700–1100 aa protein | **181 / 176 / 172** | **112 / 106 / 103** | **124 / 122 / 121** | **117 / 110 / 110** | **118 / 118 / 118** |


---



## Positive control

- Pipeline run on original Li et al. (2021) *Hydrophis cyanocinctus* -> `funannotate` proteome


| Metric                    | Benchmark | This pipeline | Difference |
| ------------------------- | --------- | ------------- | ---------- |
| V2R-family proteins       | 458       | 493           | +7.6%      |
| confident V2R transcripts | 426       | 493           | +15.7%     |
| **V2R genes**             | **~408**  | **433**       | **+6.1%**  |


---



#### Why does *H. major* have >> V2Rs?

- ***Theory* = Uncollapsed heterozygosity**
  - hmaj has the **largest assembly (2150 Mb, 90–230 Mb more than the others) and the most sequences (1111)**. Both are classic signatures of a diploid assembly that retained haplotypes instead of collapsing them. If a V2R array is represented on both a chromosome and an unplaced scaffold, every gene in it is counted twice

**TEST** = cluster each species' V2R proteins @ 99% identity & count scaffold genes that fall in a cluster with a chromosome gene

- **63 of hmaj's 198 scaffold genes are near-identical duplicates of a chromosome copy**
  - No other species shows this (*inc. curtus (east) which has more scaffold genes (240) but only 6 duplicates*)
    - *Therefore hmaj estimate -> 661 - 63 = 598 genes*


|                                               | hmaj    | hcy | hcure    | hcurw | horn |
| --------------------------------------------- | ------- | --- | -------- | ----- | ---- |
| genes on scaffolds                            | 198     | 3   | **240**  | 4     | 7    |
| of those, ≥99% identical to a chromosome copy | **63**  | 1   | **6**    | 1     | 3    |
| as % of scaffold genes                        | **32%** | –   | **2.5%** | –     | –    |




# FINAL


| stage                                                                | hmaj            | hcy             | hcure           | hcurw           | horn            |
| -------------------------------------------------------------------- | --------------- | --------------- | --------------- | --------------- | --------------- |
| EviAnn original V2R-labelled (`Similar to Vmn2r*`)                   | 624             | 402             | 448             | 479             | 434             |
| T2 (loci, tier fixed pre-repair)                                     | 607             | 393             | 400             | 483             | 390             |
| T2 genes, post model repair                                          | 655             | 430             | 416             | 509             | 423             |
| partial (tx / genes)                                                 | 255 / 239       | 176 / 164       | 137 / 121       | 197 / 188       | 127 / 117       |
| pseudogene (tx / genes)                                              | 180 / 141       | 142 / 114       | 166 / 131       | 166 / 142       | 171 / 139       |
| functional (tx / genes)                                              | 301 / 281       | 175 / 155       | 196 / 166       | 211 / 183       | 184 / 169       |
| **high confidence** (tx / genes / loci) – intact + T2 + functional   | 301 / 281 / 268 | 175 / 155 / 143 | 196 / 166 / 163 | 211 / 183 / 180 | 184 / 169 / 165 |
| **strict** = + ≥6 exons + 700–1100 aa (tx / genes / loci)            | 181 / 176 / 172 | 112 / 106 / 103 | 124 / 122 / 121 | 117 / 110 / 110 | 118 / 118 / 118 |
| eggNOG annotated / canonical domain (`7tm_3`+`ANF_receptor`+`NCD3G`) | 736 / 729       | 493 / 487       | 499 / 494       | 574 / 568       | 482 / 476       |
| genes on chromosomes                                                 | 463             | 430             | 178             | 509             | 418             |
| genes on scaffolds (unplaced)                                        | 198 (30.0%)     | 3 (0.7%)        | 240 (57.4%)     | 4 (0.8%)        | 7 (1.6%)        |
| **final V2R protein count** (transcripts, `is_v2r`=yes)              | **736**         | **493**         | **499**         | **574**         | **482**         |
| **final V2R gene count**                                             | **661**         | **433**         | **418**         | **513**         | **425**         |




## Placement per chromosome (genes)


| chrom                | hmaj    | hcy     | hcure   | hcurw   | horn    |
| -------------------- | ------- | ------- | ------- | ------- | ------- |
| ch1                  | 40      | 25      | 34      | 56      | 35      |
| ch2                  | 191     | 188     | 75      | 181     | 205     |
| ch3                  | 4       | 0       | 2       | 0       | 0       |
| ch4                  | 112     | 111     | 1       | 102     | 56      |
| ch5                  | 4       | 21      | 1       | 4       | 1       |
| ch6                  | 0       | 0       | 0       | 0       | 0       |
| ch7                  | 2       | 2       | 4       | 2       | 2       |
| ch8                  | 0       | 1       | 2       | 0       | 0       |
| chZ                  | 110     | 82      | 59      | 164     | 119     |
| scaffolds (unplaced) | 198     | 3       | 240     | 4       | 7       |
| **total genes**      | **661** | **433** | **418** | **513** | **425** |


