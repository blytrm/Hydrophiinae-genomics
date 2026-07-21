# Hydrophiinae Genomics

**Resolution of sea snake chemoreceptomes: _In silico_ reconstruction of vomeronasal type-2 receptor (V2R) gene-family dynamics in _Hydrophiinae_.**

Honours thesis · Sanders Laboratory, School of Biological Sciences, University of Adelaide
Billy Trim (`a1864358`) · Supervisors: **Dr Alastair Ludington** & **Prof. Kate Sanders**

---

  | `Index` |
  | --- |
   - [Assembly](./1_assembly)
        - [QC](./1_assembly/QC)
        - [Pipeline](./1_assembly/pipeline)
        - [Results](./1_assembly/results)
        - [Sequence Resolution/Investigation](./1_assembly/sequence_resolution)
             - [Coverage](./1_assembly/sequence_resolution/coverage)
   - [Annotation](./2_annotation)
        - [V2R Structure/Homology-based Proxy Search (Conserved Protein Domain)](./2_annotation/7tm-proxy-anno)
        - [V2R Probabalistic/Homology-based Proxy Search (Hidden Markov Model)](./2_annotation/pHMM-proxy-anno)
   - [Repeat Annotation](./3_repeats)
        - [EDTA_RM](./3_repeats/EDTA_RM)
        - [LTR & V2R Tandem Array Association](./3_repeats/LTR-density_atV2Rs) <- **`technical report`**
             - [code_&scripts](./3_repeats/LTR-density_atV2Rs/code_&scripts)
             - [tech-report-pages](./3_repeats/LTR-density_atV2Rs/tech-report-pages)
        - [RROC-AUC_array-prediction](./3_repeats/RROC-AUC_array-prediction)
   - [V2R Curation + Investigation](./4_v2rs)

---

### The problem // Preliminary Work (2025)
Published V2R repertoire sizes for a **single clade** span more than an order of magnitude:

| Study | Assembly / species | V2R loci | Reading |
|---|---|---:|---|
| Kishida et al. (2019) | early sea-snake assemblies | ~137 | marine **contraction** |
| Li et al. (2021) | _H. cyanocinctus_ (HCya_v2) | >1,411 | **expansion** |
| Policarpo et al. (2024) | marine tetrapods (pre-2021 assemblies) | — | convergent **contraction** across ALL marine tetropods |
| **Trim 2026 {[_report_](https://medium.com/@billyt13/parameter-driven-refinement-of-a-sea-snake-genome-assembly-enhances-structural-accuracy-and-enables-715b5c2ecf95?sharedUserId=billyt13)} _preliminary_** | **_H. major_, reassembled** | **~2,370** | highest estimate |

These numbers cannot all be biology. **The central thesis of this project is that assembly-quality artefacts — collapsed tandem arrays, Z-chromosome miscounting, lineage-specific sequence invisible to a linear reference —> the dominant driver of this variation, not genuine copy-number differences between recently (1MYA) diverged species**

[Results showed:](https://github.com/blytrm/Hydrophiinae-chromosome-assembly-investigation)
- BUSCO completeness **96.0% → 99.24%**
- duplication rate roughly **halved**; error rate **reduced**
- chromosome-2 V2R loci **340 → 570** (genome-wide ≈ **2,370**, tBLASTn hits)

That gives a rare labelled pair —> *same genome, two assembly qualities, a quantified annotation difference* — which anchors the rest of the thesis and scales it across the _Hydrophis_ panel.

---

## Research Aims
   1. Generate improved chromosome-scale genome assemblies and annotations for six *Hydrophis* species.
        * [Genome Assembly](./1_assembly)
        * [Gene Annotation](./2_annotation)
        * [Repeat Annotation](./3_repeats)
   2. Produce a high-confidence, standardised V2R gene catalogue across the full species panel using a consistent layered annotation framework.
        * [V2R Annotation](./4_v2rs) 
   3. Model V2R gene family dynamics across the radiation and marine transition, and test structural genomic predictors of copy-number variation.

---

## Project Hypotheses

* **H1:** Improved [assemblies](./1_assembly/results) will reveal V2R repertoires larger than published—confirming that assembly artefacts systematically deflate existing counts, and hence V2Rs are underestimated in the literature.
* **H2:** *Hydrophis* species show elevated rates of V2R gene gain relative to terrestrial elapid outgroups and to neutral birth-death expectation, as an adaptive response to the marine environment.
* **H3:** Expanded V2R gene lineages in *Hydrophis* exhibit signatures of positive selection consistent with diversifying selection maintaining receptor specificity breadth across the radiation.
* **H4:** Local TE density in V2R-clusters is positively associated with V2R copy number across *Hydrophis* species, consistent with TE-mediating recombination, regulation and driving tandem array expansion.
     * >[LTR Density & V2R Gene Array Association (*Technical Report Investigation*)](./3_repeats/LTR-density_atV2Rs)

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
---
---
---


# Progress + notes
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

