## Hi-C Scaffolding 

Hi-C scaffolding of contigs with `HapHiC`
**input** -> polished de novo contig assembly

---

**param guide**
definition of `HapHiC` parameters used:

parameter space exploration and optimisation was done previously with _Hydrophis major_ genome.

| Parameter | Value | Purpose | Rationale |
|-----------|-------|---------|-----------|
| **`--correct_nrounds`** | **`4`** | Rounds of Hi-C–based misjoin (chimeric contig) detection and breaking, applied before clustering | Targets inherited chimeric joins — notably chr2<->Z fusions in V2R-dense regions — that the fragmentation-permissive [hifiasm](1_assembly/pipeline/01_denovo-asm/hifiasm.sh): `no-post-join`, `lowQ50`) input deliberately left unresolved for downstream correction. Aggressive for a HiFi assembly; retained here because chimerism is the central assembly problem under investigation |
| **`--mutprob`** | **`0.15`** | Mutation probability in the ALLHiC genetic-algorithm optimiser (contig ordering/orientation) | Sustains exploration of the ordering solution space; mitigates premature convergence to local optima without destabilising the search |
| **`--npop`** | **`200`** | Population size of the GA optimiser | Larger population broadens the search over candidate orderings, improving robustness for a high-contig-count (n=1,896) input at modest runtime cost |
| **`--N_CHROM`** | **`14`** | Expected chromosome number; target for inflation selection and final pseudomolecule count | Karyotype provisional (prior estimate ~18 pairs; current evidence suggests fewer). Set liberally at 14: over-specification is tolerated, as the reassignment step concatenates surplus groups, whereas under-specification risks fusing distinct chromosomes |
| **`--min_inflation`** | **`1.8`** | Lower bound of the MCL inflation sweep | Coarse end of the search; permits recovery of larger, contiguous chromosomal groups |
| **`--max_inflation`** | **`10.0`** | Upper bound of the MCL inflation sweep | Raised from the 3.0 default. Higher inflation yields finer clustering; extending the ceiling lets the sweep reach inflations that resolve distinct chromosomes under the elevated inter-chromosomal Hi-C background expected in repeat-dense, structurally complex termini |
| **`--topN`** | **`25`** | Number of top Hi-C–linked neighbours retained per contig during ordering | Bounds the search space, keeping ordering tractable for the fragmented input while preserving the strongest contact signal |
| **`--bin_size`** | **`200`** | Bin resolution for the correction-step link-density scan used to call breakpoints | Coarser bins reduce noise sensitivity in breakpoint detection, appropriate for moderate Hi-C depth |

---

#### Why These Specific Values for Your Sea Snake Assembly:

- **`N_CHROM=14`** — Helps guide contig scaffolding (although suspected under-estimated)
- **`correct_nrounds=4`** — Iterative refinement improves layout consistency
- **`bin_size=200`** — Balances resolution vs. Hi-C coverage; coarser bins are more robust to sequencing noise
- **`npop=200` + `mutprob=0.15`** — Moderate genetic algorithm parameters; aggressive enough to escape local optima but not so chaotic as to be slow
- **Inflation bounds (1.8–10.0)** — Prevents both under-clustering and over-splitting of Hi-C contacts into spurious scaffolds

---

#### Hi-C scaffolding parameter exploration + optimisation:
Report -> https://github.com/blytrm/Hydrophiinae-chromosome-assembly-investigation
      ---> https://medium.com/@billyt13/parameter-driven-refinement-of-a-sea-snake-genome-assembly-enhances-structural-accuracy-and-enables-715b5c2ecf95

[Download/View PDF](https://raw.githubusercontent.com/blytrm/Hydrophiinae-chromosome-assembly-investigation/main/BTrim_Parameter-Driven-Refinement-of-a-Sea-Snake-Genome-Assembly-Enhances-Structural-Accuracy-and-Enables-Improved-Gene-Quantification.pdf)
