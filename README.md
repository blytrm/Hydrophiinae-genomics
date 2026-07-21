# Hydrophiinae Genomics

**Computational resolution of sea snake chemoreceptomes: phylogenomic reconstruction of vomeronasal type-2 receptor (V2R) gene-family dynamics in _Hydrophiinae_.**

Honours thesis · Sanders Laboratory, School of Biological Sciences, University of Adelaide
Billy Trim (`a1864358`) · Supervisors: **Dr Alastair Ludington** & **Prof. Kate Sanders**

---

### The problem & Preliminary Work (2025)
Published V2R repertoire sizes for a **single clade** span more than an order of magnitude:

| Study | Assembly / species | V2R loci | Reading |
|---|---|---:|---|
| Kishida et al. (2019) | early sea-snake assemblies | ~137 | marine **contraction** |
| Li et al. (2021) | _H. cyanocinctus_ (HCya_v2) | >1,411 | **expansion** |
| Policarpo et al. (2024) | marine tetrapods (pre-2021 assemblies) | — | convergent **contraction** across ALL marine tetropods |
| **[Trim](https://medium.com/@billyt13/parameter-driven-refinement-of-a-sea-snake-genome-assembly-enhances-structural-accuracy-and-enables-715b5c2ecf95?sharedUserId=billyt13) (2026), _preliminary_** | **_H. major_, reassembled** | **~2,370** | highest estimate to date |

These numbers cannot all be biology. **The central thesis of this project is that assembly-quality artefacts — collapsed tandem arrays, Z-chromosome miscounting, lineage-specific sequence invisible to a linear reference — are the dominant driver of this variation, not genuine copy-number differences between species.**

The strongest evidence is a controlled, same-genome comparison. Reassembling _Hydrophis major_ during 2025 held the biology fixed and varied only assembly quality:

- BUSCO completeness **96.0% → 99.24%**
- duplication rate roughly **halved**; error rate **reduced**
- chromosome-2 V2R loci **340 → 570** (genome-wide ≈ **2,370**, tBLASTn hits)

That gives a rare labelled pair — *same genome, two assembly qualities, a quantified annotation difference* — which anchors the rest of the thesis and scales it across the _Hydrophis_ panel.

---

### Project Aims
* **Improved reference assemblies** across the species panel (_H. major_, _H. curtus_, _H. ornatus_, _H_cyanocinctus_).
* **A layered V2R annotation catalogue** 
* **Phylogenetic + evolutionary dynamics**
* **Role of transposable elements in gene family expansion**

## Hypotheses

- **H1 — Assembly quality confounds V2R counts.** Higher-quality reassembly recovers more V2R loci; BUSCO completeness correlates positively with locus recovery, and the original assemblies undercounts.
- **H2 — LTR co-expansion.** LTR-retrotransposon density in the neighbourhood of V2R loci predicts copy number after controlling for assembly quality (negative-binomial GLM; spatial-autocorrelation checked).
- **H3 — Lineage-specific birth rates.** Marine _Hydrophis_ lineages show elevated V2R birth rates relative to terrestrial elapid outgroups (CAFE5), contrary to a blanket marine-contraction narrative.
- **H4 — Positive selection on ligand binding.** Retained functional V2Rs in marine lineages carry episodic positive selection (ω > 1) concentrated in the ligand-binding domain, consistent with ligand-driven diversification rather than drift.

---
Genomic resources, scripts, plots and comparative analyses for Hydrophiinae sea snake research, for my honours thesis
>*Billy Trim · Supervisors: Dr Alastair Ludington & A/Prof. Kate Sanders*

## Project Description 
```text
abstract
```

### Project Aims
```text
- aims
- hypotheses
```

---
### Outline
# [Genome Assembly](./1_assembly)
# Genome Annotation
* ## [Repeat Annotation](./3_repeats)
    * ### [LTR Density & V2R Gene Array Association (*Technical Report Investigation*)](./3_repeats/LTR-density_atV2Rs)
* ## [Gene Annotation](./2_annotation)
* ## [V2R Annotation](./4_v2rs)

#### Comparative Genomics

#### Phylogenetics & Evolution

#### Vomeronasal Type-2 Receptors

---
---


# Progress
```
    ☑️ = done
    ✅ = uploaded to repo
```

### Genome Assembly
- outline
- steps 1-10 ✅
- qc✅
- plots☑️

### Genome Annotation
- tblastn☑️
- phmm☑️
- eviann
- repeats☑️

### Comparative Genomics
- improvement☑️
- across species
- plots
- synteny
- conserved / diverged regions
- SVs
- v2rs

### Vomeronasal type-2 Receptor Genes & Repetitive Sequence Association
- permutation tests☑️
- regression☑️
- on eviann v2rs + all species

### Phylogenetics & Evolution
- selection testing
- birth death modelling
- ortho/paralogous

### Vomeronasal type-2 Receptors
- rna seq -> tissue specificity
- cnv vs expression
- protein prediction
- subtype/family clustering
- extracellular codons
- ligand binding pocket
- comparison w terrestrials


---

## Background & Relevant Biology


...

---


<table>
  <tr>
    <td colspan="3" align="center"><b>Conference Poster — GSA 2026</b><br><img src="./resources/GSA26-BTRIMa.png" width="400"></td>
  </tr>
  <tr>
    <td colspan="3" align="center"><b>Literature Review &amp; Proposal</b></td>
  </tr>
  <tr>
    <td align="center"><sub>Page 1</sub><br><img src="./resources/pages/lit_rev-01.png" width="320"></td>
    <td align="center"><sub>Page 11</sub><br><img src="./resources/pages/lit_rev-11.png" width="320"></td>
    <td align="center"><sub>Page 20</sub><br><img src="./resources/pages/lit_rev-20.png" width="320"></td>
  </tr>
  <tr>
    <td colspan="3" align="center"><b>Technical Report</b></td>
  </tr>
  <tr>
    <td align="center"><sub>Page 1</sub><br><img src="./resources/pages/page-01.png" width="320"></td>
    <td align="center"><sub>Page 4</sub><br><img src="./resources/pages/page-04.png" width="320"></td>
    <td align="center"><sub>Page 8</sub><br><img src="./resources/pages/page-08.png" width="320"></td>
  </tr>
</table>

---

<details>
  <summary><h2>Click to view the full literature review & project proposal report</h2></summary>

![Page 1](./resources/pages/lit_rev-01.png)

![Page 2](./resources/pages/lit_rev-02.png)

![Page 3](./resources/pages/lit_rev-03.png)

![Page 4](./resources/pages/lit_rev-04.png)

![Page 5](./resources/pages/lit_rev-05.png)

![Page 6](./resources/pages/lit_rev-06.png)

![Page 7](./resources/pages/lit_rev-07.png)

![Page 8](./resources/pages/lit_rev-08.png)

![Page 9](./resources/pages/lit_rev-09.png)

![Page 10](./resources/pages/lit_rev-10.png)

![Page 11](./resources/pages/lit_rev-11.png)

![Page 12](./resources/pages/lit_rev-12.png)

![Page 13](./resources/pages/lit_rev-13.png)

![Page 14](./resources/pages/lit_rev-14.png)

![Page 15](./resources/pages/lit_rev-15.png)

![Page 16](./resources/pages/lit_rev-16.png)

![Page 17](./resources/pages/lit_rev-17.png)

![Page 18](./resources/pages/lit_rev-18.png)

![Page 19](./resources/pages/lit_rev-19.png)

![Page 20](./resources/pages/lit_rev-20.png)

![Page 21](./resources/pages/lit_rev-21.png)

![Page 22](./resources/pages/lit_rev-22.png)

![Page 23](./resources/pages/lit_rev-23.png)

![Page 24](./resources/pages/lit_rev-24.png)

![Page 25](./resources/pages/lit_rev-25.png)

![Page 26](./resources/pages/lit_rev-26.png)

![Page 27](./resources/pages/lit_rev-27.png)

![Page 28](./resources/pages/lit_rev-28.png)

![Page 29](./resources/pages/lit_rev-29.png)

![Page 30](./resources/pages/lit_rev-30.png)

![Page 31](./resources/pages/lit_rev-31.png)

![Page 32](./resources/pages/lit_rev-32.png)

![Page 33](./resources/pages/lit_rev-33.png)

![Page 34](./resources/pages/lit_rev-34.png)

[PDF download](./resources/BTrim_Literature_Review_Proposal.pdf)

</details>

---

<details>
  <summary><h2>Click to view the full technical report (15 pages)</h2></summary>

![Page 1](./resources/pages/page-01.png)

![Page 2](./resources/pages/page-02.png)

![Page 3](./resources/pages/page-03.png)

![Page 4](./resources/pages/page-04.png)

![Page 5](./resources/pages/page-05.png)

![Page 6](./resources/pages/page-06.png)

![Page 7](./resources/pages/page-07.png)

![Page 8](./resources/pages/page-08.png)

![Page 9](./resources/pages/page-09.png)

![Page 10](./resources/pages/page-10.png)

![Page 11](./resources/pages/page-11.png)

![Page 12](./resources/pages/page-12.png)

![Page 13](./resources/pages/page-13.png)

![Page 14](./resources/pages/page-14.png)

![Page 15](./resources/pages/page-15.png)

[PDF download](./resources/Billy%20Trim%20Technical%20Report.pdf) · [Analysis scripts](https://github.com/blytrm/Hydrophiinae-genomics/tree/main/3_repeats/LTR-density_atV2Rs)

</details>

---

