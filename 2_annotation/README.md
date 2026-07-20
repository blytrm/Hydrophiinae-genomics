# ANNOTATION
---
From [Assembly Pipeline](../1_assembly):
```text
--> 5 re-assembled genomes
  --> All w improved quality
```
---

## Structural & Functional Annotation
* [EviAnn Pipeline](./eviann-anno)

---
### Initial Proxy V2R Annotation
  * Used for initial/broad V2R count estimations, to identify tandem arrays, array complexity, and in [*'Technical Report'*](../3_repeats/LTR-density_atV2Rs)
1. [Conserved Domain Genome Scan](./7tm-proxy-anno)
   * **The 7-transmembrane domain of V2Rs was utilised to scan genomes for count estimates**:
       >The 7-transmembrane domain is an ideal genomic search target because it is the most highly conserved structural element across all vomeronasal type-2 receptors and related Class C GPCRs. This domain spans the cell membrane and is under strong functional constraint—it must maintain the precise topology and chemistry required for signal transduction across the lipid bilayer. Because the 7TM architecture is essential for the basic mechanism of G-protein coupled receptor signalling and is shared across all functional V2Rs regardless of species or ligand specificity*

2. [Homology-based Search](./pHMM-proxy-anno)
   * **_Profile Hidden Markov Model_** **<-- curated from consensus MSA from high-quality vertebrate V2R protein databases**
       >Profile HMMs capture position-specific residue preferences and indel patterns across a curated V2R alignment, making them more sensitive to remote homologs than pairwise BLAST while remaining specific to the V2R family. Seeds were drawn from reviewed vertebrate V2R databases (mouse Vmn2r, rat Vom2r, Swiss-Prot, squamate and Xenopus outgroups), then validated against a Class C GPCR decoy panel (TAS1R, CaSR, GRM, GABBR, GPRC5/6) so that gathering thresholds sit above the strongest false-positive scores before genome search.*
---    
