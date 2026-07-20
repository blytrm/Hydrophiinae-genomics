# Association of Repetitive Elements at Chemosensory Receptor Genes in *Hydrophis* Genomes
### *`Are LTRs Enriched at Candidate V2R Loci? Could This Help Explain Gene Family Expansion & Copy-Number Variation?`*
**BIOL 4001 — Advanced Biological Sciences · Technical Report**  
---

[PDF download](./BILLY_TRIM-TECH_REPORT.pdf) · [Analysis scripts](https://github.com/blytrm/Hydrophiinae-genomics/tree/main/3_repeats/LTR-density_atV2Rs)


# `TO DO` = **_Repo map_**

---

## Project outline

*Hydrophis* sea snakes shifted from air- to water-borne chemosensation after diverging from terrestrial Australian elapids (~9–18 MYA). Squamates expand vomeronasal type-2 receptor (**V2R**) gene families in tandem arrays; long terminal repeat retrotransposons (**LTRs**) can supply substrates for non-allelic homologous recombination and thus copy-number change. This report tests whether LTR density at V2R loci exceeds spatial null expectations on chromosomes **C2** and **CZ** of a reassembled *H. ornatus* genome, and whether local LTR density predicts larger paralog arrays.

### Aims

| ID | Hypothesis | Approach |
|----|------------|----------|
| **H1** | LTR coverage at V2R loci/clusters exceeds a chromosome-aware spatial null | Permutation enrichment (`regioneR`), cluster-level randomisation, cyclic LTR-track rotation |
| **H2** | V2R cluster paralog count increases with local LTR density | Negative binomial GLM + BH multiple-testing correction |

### Workflow

1. -> ![Reassembled](../../1_assembly) *H. ornatus* genome
2. Repeat Annotation (![EDTA + RepeatMasker](../../3_repeats))
3. Proxy V2R Candidate Loci:
  -> ![tBLASTn domain scan](../../2_annotation/7tm-proxy-anno)
  -> ![profile Hidden Markov Model via HMMER](../../2_annotation/pHMM-proxy-anno)
6. Permutation Enrichment Tests
7. Cluster-Level Enrichment + Cyclic Permutation Enrichment Tests
8. Dose-Response Regression Prediction

### Null hierarchy (H1)

1. **Null A** — genome-wide shuffle *(biologically invalid foil; V2Rs are chromosome-structured)*
2. **Null B** — shuffle within C2 / CZ only
3. **Null C** — GC-matched resampling (`matchRanges` logic)
4. **Null D** — exclude low-mappability / coverage-anomaly regions
5. **Cyclic permutation** — rotate LTR track; keep V2Rs fixed *(most conservative)*

### Headline results

- Genome-wide Null A **inflates** enrichment (~3.6×) relative to chromosome-matched nulls (~1.4–1.6×).
- Cluster-level enrichment survives on **C2** (all nulls significant); **CZ** is method-sensitive (all nulls except cyclic permutation significant).
- **No robust monotonic dose–response** between LTR density and cluster size after BH correction.
- Enrichment is concentrated in short `LTR/unknown` fragments — compatible with degraded remnants or annotation artefacts, not necessarily intact drivers of duplication.
- Strong co-localisation of V2Rs with low-mappability sequence is a major assembly confound (esp. CZ).

### Repo map (analysis code)

```text
Hydrophiinae-genomics/
└── 3_repeats/LTR-density_atV2Rs/
    ├── annotation/          # EDTA, RepeatMasker, V2R calls
    ├── tech-report-pages/   # Document pages, figures, LaTeX document code, pdf
    ├── dose_response/       # NB-GLM
    ├── figures/             # 
    └── README               # method notes
```

---

![Page 1](./tech-report-pages/page-01.png)
<details>
  <summary><h2>->Click to view the full technical report (15 pages)</h2></summary>

![Page 2](./tech-report-pages/page-02.png)

![Page 3](./tech-report-pages/page-03.png)

![Page 4](./tech-report-pages/page-04.png)

![Page 5](./tech-report-pages/page-05.png)

![Page 6](./tech-report-pages/page-06.png)

![Page 7](./tech-report-pages/page-07.png)

![Page 8](./tech-report-pages/page-08.png)

![Page 9](./tech-report-pages/page-09.png)

![Page 10](./tech-report-pages/page-10.png)

![Page 11](./tech-report-pages/page-11.png)

![Page 12](./tech-report-pages/page-12.png)

![Page 13](./tech-report-pages/page-13.png)

![Page 14](./tech-report-pages/page-14.png)

![Page 15](./tech-report-pages/page-15.png)

</details>

---

## License/citation

Course technical report (BIOL 4001). Word count: 2,850.
