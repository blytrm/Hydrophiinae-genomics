# Genome Assembly Methods & Protocol

Reference-quality, chromosome-scale assembly of Hydrophis sea snake genomes (H. major, H. cyanocinctus, H. ornatus, and the two H. curtus lineages). This document records the assembly workflow step by step, the tool and parameter choices, and the justification for each decision. Each numbered step corresponds to a script in [...].
All steps are run under SLURM on an HPC cluster. Polishing precedes scaffolding, and scaffolding precedes gap-filling and a final polish, because each downstream tool assumes its input is locally correct at the base level (Rhie et al., 2021). QC metrics (aiming to cover contiguity, completeness and assembly correctness; scripts available here) are recomputed after every step to confirm monotonic improvement rather than regression (see table X for these results).

[flow chart]
[add in versions]

### De novo Contig Assembly
The genome assembly pipeline begins from raw sequencing reads being assembled into contigs, whereby the assembler tool is read-type dependent [something about POC/exploration validation & flye experimentation].

HiFi reads (PacBio) → Hifiasm
ONT & CLR reads → NextDenovo
XXX → Flye
H. ornatus was assembled with Flye to quantify the effect of assembler choice on final quality.

### Long-read Misassembly Detection & Correction with Inspector

Raw long reads are first subsampled to 40× coverage against the contig set with rasusa (Hall, 2022), then Inspector (Chen et al., 2021) maps the reads back to the assembly, calls structural and small-scale errors from read–assembly disagreement, and emits a corrected FASTA via inspector-correct.py.
**Justification**. Inspector evaluates an assembly against the reads that built it and corrects read-inconsistent contigs, removing errors that would otherwise propagate into scaffolding. Subsampling to 40× before mapping caps runtime and memory, and avoids depth-driven bias from over-represented regions without sacrificing the signal needed for reliable error calls.

### Structural Error Breaking with CRAQ
Short reads (Illumina) and long reads are aligned to the Inspector-corrected assembly, and CRAQ (Li et al., 2023) breaks contigs at positions where regional and clip-based error signals jointly indicate a structural error.

**Justification**. CRAQ quantifies clip-based regional (CRH) and clip-based assembly (CRE) error signals at single-nucleotide resolution and breaks contigs at conflict breakpoints, correcting structural errors that Inspector's read-consistency model can miss. Using both short and long reads gives orthogonal evidence: short reads localise base-level conflict precisely, long reads span and resolve repeat-induced misjoins. samblaster marks PCR/optical duplicates and `-F 3340` discards unmapped, secondary, supplementary, and duplicate alignments so only primary, properly-mapped reads inform the error model.

### Pre-Scaffolding Polish with NextPolish

NextPolish (Hu et al., 2020) polishes the structurally corrected contigs using Illumina short reads (`sgs_option`, `task = best`, `-max_depth 100 -bwa`).
Justification. Base-level correction is performed before scaffolding because Hi-C scaffold-breaking and curation tools assume the input sequence is locally correct (Rhie et al., 2021); indels and base errors left in place would corrupt restriction-site detection and contact-density estimation downstream. Short reads provide the highest per-base accuracy for resolving residual indels and substitutions. `task = best` runs NextPolish's recommended multi-round short-read protocol; `-max_depth 100` caps pile-up depth to bound memory.

### Hi-C Scaffolding with HapHiC

Hi-C reads are aligned with bwa mem -5SP, filtered, and passed to HapHiC (Zeng et al., 2024) to cluster, order, and orient contigs into the expected number of pseudo-chromosomes (N_CHROM).
Justification. Hi-C contact frequency reflects 3D physical proximity, which decays with genomic distance and so encodes long-range ordering information that long reads alone cannot supply. HapHiC is reference-independent and allele-aware, avoiding bias from a (non-existent, for Hydrophis) reference and handling heterozygous, repeat-rich genomes robustly. `bwa mem -5SP` is the standard Hi-C alignment mode (skips mate rescue/pairing so chimeric Hi-C reads are reported as split alignments); filter_bam --nm 3 removes alignments with >3 mismatches; the `--RE` motif set is the Arima two-/multi-enzyme signature used to score restriction-site density during scaffolding.

Hi-C reads are re-aligned `(bwa mem -5SP | samblaster | samtools view -F 3340)` to the HapHiC scaffolds to produce the BAM used for error detection in Step 6.
**Justification**. AutoHiC's error model operates on a contact map built from the scaffolded coordinate system, so reads must be re-mapped to the scaffolds rather than the pre-scaffold contigs. The identical FASTA used here must be reused in Step 6 to keep coordinates consistent.

### Deep-Learning based Scaffold Correction with AutoHiC’s MoreHiC module

The scaffold AGP/FASTA is converted to JuiceBox .assembly format, the BAM is converted to a .hic contact matrix (via matlock and the 3D-DNA visualiser), and AutoHiC's onehic module (Jiang et al., 2024) detects and corrects scaffold-level errors directly from the contact-map image, then writes a corrected FASTA.

```bash
agp2assembly.py / makeAgpFromFasta.py     # FASTA/AGP -> .assembly
matlock bam2 juicer  +  run-assembly-visualizer.sh   # BAM -> .hic
onehic.py -hic <hic> -asy <assembly> -p models/error_model.pth   # detect+correct
juicebox_assembly_converter.py            # adjusted .assembly -> FASTA
```

**Justification**. Scaffold-level errors: inversions, misorientations, and misjoins, manifest as recognisable patterns in Hi-C contact heatmaps. AutoHiC applies a trained deep-learning detector to these images, reporting >90% correction accuracy, which removes the most error-prone class of mistakes before a human curator looks at the map and reduces the manual workload to genuine edge cases.

```bash
Agp_tools
AGP2fa ow whatever
```

### Manual Curation with PretextView

Hi-C reads are aligned to the AutoHiC-corrected assembly and rendered as a PretextMap contact map; the assembly is then manually inspected, broken, inverted, and re-oriented in PretextView (Harry, 2022).

```bash
bwa mem -5SP ... | samblaster --ignoreUnmated \
  | samtools view -F 2316 | PretextMap -o <prefix>.pretext \
        --sortby length --sortorder descend --mapq 10
```

**Justification**. Automated tools do not catch every error. Howe et al. (2021) showed manual curation recovers ~6% additional placed sequence and corrects a median of ~4 misjoins per vertebrate assembly, errors largely invisible to N50 and BUSCO yet directly damaging to syntenic inference across tandem gene arrays, which is the central biological question here. `-F 2316` and `--mapq 10` retain only confidently, uniquely mapped read pairs so the curator views a clean, high-confidence contact signal.

### Gap Filling with TGS-GapCloser

The curated assembly's residual scaffold gaps (N-runs introduced during scaffolding) are closed with … 
Justification. Scaffolding joins contigs across gaps of estimated size; filling those gaps with real long-read sequence converts ambiguous N-runs into called bases, improving contiguity (scaftig N50) and completeness without altering the validated scaffold order/orientation. A consolidated alternative, gap-pol.sh, performs the equivalent step with TGS-GapCloser (Xu et al., 2020), which is purpose-built for closing gaps in large genomes from low-coverage error-prone long reads, followed immediately by the Step 9 polish.

### Final Polish with NextPolish/NextPolish2

A final NextPolish round is run on the curated, gap-filled assembly. This round uses only short reads (sgs_option), with minimap2 -x map-pb/map-ont for the long-read alignment.
Justification. Gap-filling inserts new long-read sequence at gap junctions, and curation creates new contig joins; both introduce fresh, unpolished bases and junction indels. A final polish corrects these, with short reads supplying per-base accuracy and long reads disambiguating repeat-region context, yielding the best achievable consensus QV before release.

Alternative Pipeline for Hydrophis elegans
Ragtag scaffolded off cyanocinctus

## Quality Control

QC follows the "3C" framework: contiguity, completeness and correctness. Key metrics are recomputed after each step to verify gradual, monotonic improvement:
- **Contiguity**: total size, contig/scaffold counts, N50, and auN (area under the N-curve; a length-weighted contiguity average), using QUAST (Gurevich et al., 2013).
- **Completeness**: single-copy ortholog recovery against squamata_odb12 via compleasm (Huang & Li, 2023), plus k-mer completeness and consensus quality value (QV) from MerquryFK against FastK k-mer databases.
- **Correctness**: MerquryFK QV and mosdepth coverage statistics (e.g. breadth and uniformity of coverage; Pedersen & Quinlan, 2018), confirming that more reads map, and map cleanly (without excessive clipping or multi-mapping), to the assembly.

**Justification**. N50 alone rewards over-joining and can rise while correctness falls; pairing contiguity with read-based completeness and correctness metrics detects regressions (e.g. a curation edit that inflates N50 but creates a chimera) that any single metric would hide.
 Other tools / metrics experimented with …


#### References
