# Genome Annotation Evidence Summary
**Species:** *Hydrophis cyanocinctus* (hcy)  
**Genome:** `10-hcy-final.orignames.fa` (ch1–ch14, contig1–contig34; 48 sequences)  
**Annotation tool:** EviAnn (job `10840840`, submitted 2026-05-11)  
**Pipeline:** `/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/pipeline/`

---

## 1. RNA-seq Evidence

### 1.1 Summary

| | Count |
|---|---|
| Raw samples (input) | 31 |
| Samples passing QC | 30 |
| Samples excluded | 1 (`350843-vitFol` — persistent OOM during alignment) |
| Aligned BAMs | 30 |
| Total BAM size | ~123 GB |

### 1.2 QC Pipeline (step 02)

Tool versions from conda env `rnaqc`:

| Step | Tool | Purpose |
|---|---|---|
| 1 | FastQC | Pre-QC quality assessment |
| 2 | fastp | Adapter trimming, quality filtering |
| 3 | BBDuk | rRNA depletion |
| 4 | Kraken2 | Contamination screening |
| 5 | FastQC | Post-QC quality assessment |
| 6 | MultiQC | Aggregated QC report |

**fastp parameters:**

| Parameter | Value |
|---|---|
| Min read length | 50 bp |
| Average quality threshold | 15 |
| Qualified quality threshold | 15 |
| Max N bases | 5 |

**Databases:**

| Database | Path |
|---|---|
| Adapter sequences | `/hpcfs/.../bbmap/resources/adapters.fa` |
| rRNA reference | `/scratchdata1/.../dbs/rrna_refs/rrna_combined.fa` (16S/18S/23S/5S) |
| Kraken2 DB | `/scratchdata1/.../dbs/kraken2/k2_db` |

Cleaned reads: `/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/results/rna-qc/cleaned/`

### 1.3 Alignment (step 04)

**Tool:** HISAT2 (16 threads) → samtools sort  
**Index:** built from `10-hcy-final.orignames.fa`  

| Parameter | Value |
|---|---|
| Min intron length | 20 bp |
| Max intron length | 1,388,389 bp |

BAMs: `/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/results/rna-align/bam/`

### 1.4 Samples Used

#### In-house — *H. cyanocinctus* (target species)

| Sample ID | Tissue |
|---|---|
| `hcy-brain` | Brain |
| `hcy-parietal` | Parietal eye |
| `hcy-pineal` | Pineal gland |

#### In-house — mixed Hydrophiinae

| Sample ID | Species (inferred) | Tissue |
|---|---|---|
| `350840-liver` | — | Liver |
| `350841-skin` | — | Skin |
| `350842-brain` | — | Brain |
| `ala-KLS0468-vmo` | *Acalyptophis* sp. | Vomeronasal organ |
| `C9FMBANXX-2-B223_B_RHINOCEROS_VNO_TGACCA_L002` | *B. rhinoceros* | VNO |
| `C9FMBANXX-2-NS313_N_SCUTATUS_VNO__F_ATCACG_L002` | *Notechis scutatus* | VNO |
| `hda-eye` | *H. darwinensis* | Eye |
| `hja-eye` | *H. jaegeri* | Eye |
| `hma-eye-L2` | *H. major* | Eye |
| `hma-KLS0460-heart` | *H. major* | Heart |
| `hma-KLS0460-liver` | *H. major* | Liver |
| `hma-KLS0460-testis` | *H. major* | Testis |

#### Public — *Hydrophis curtus* (SRA)

| Accession | Tissue |
|---|---|
| SRR11659669 | — |
| SRR11659670 | — |
| SRR11659671 | — |
| SRR12882912 | — |
| SRR13325964 | — |

#### Public — *Hydrophis cyanocinctus* (SRA)

| Accession | Tissue |
|---|---|
| SRR11659657 | — |
| SRR11659658 | — |
| SRR11659659 | — |

#### Public — *Hydrophis melanocephalus* (DDBJ)

| Accession | Tissue |
|---|---|
| DRR147123 | — |
| DRR147124 | — |
| DRR147125 | — |
| DRR147126 | — |

#### Public — *Hydrophis platurus* (SRA)

| Accession | Tissue |
|---|---|
| SRR3955060 | — |
| SRR3955061 | — |

#### Public — *Hydrophis schistosus* (SRA)

| Accession | Tissue |
|---|---|
| SRR16480669 | — |

#### Excluded

| Sample ID | Reason |
|---|---|
| `350843-vitFol` | Repeated OOM failure during HISAT2 alignment (>128 GB); excluded from annotation |

---

## 2. Protein Evidence

### 2.1 Summary

| | Value |
|---|---|
| Input files | 3 |
| Total proteins | 891,380 |
| Total amino-acid residues | 476,490,871 |
| Empty records skipped | 0 |
| Duplicate IDs renamed | 0 |
| Output size | 527 MB |
| Min protein length | 1 aa |
| Max protein length | 39,205 aa |

Output: `/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/results/protein-qc/proteins.clean.faa`

### 2.2 Input Files (step 03)

| File | Records | Description |
|---|---|---|
| `txid8570_refseq_proteins.faa` | 890,908 | Squamata RefSeq proteins (NCBI taxid 8570) |
| `squamata_v2rs.faa` | 455 | Squamata V2R chemoreceptor proteins |
| `mus_vmn2r_specific.faa` | 17 | *Mus musculus* VMN2R-specific proteins |

Cleaning: non-standard amino acids replaced with X; gap/stop characters stripped; duplicate FASTA IDs de-duplicated with source suffix.

---

## 3. EviAnn Run

| Parameter | Value |
|---|---|
| Job ID | `10840840` |
| Partition | icelake |
| Threads | 16 |
| Memory | 60 GB |
| Time limit | 48 h |
| Genome | `10-hcy-final.orignames.fa` |
| RNA BAMs | 30 sorted BAMs (`bam` tag) |
| Proteins `-p` | `proteins.clean.faa` |
| UniProt `-s` | not supplied |
| Work dir | `/hpcfs/.../annotation/results/eviann/hcy/` |

---

## 4. File Locations

| Item | Path |
|---|---|
| Genome (annotation) | `/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/files/hcy.backup0.backup0/10-hcy-final.orignames.fa` |
| Genome (renamed chr) | `/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/files/hcy.backup0.backup0/10-hcy-final.fa` |
| Cleaned RNA reads | `results/rna-qc/cleaned/` |
| Aligned BAMs | `results/rna-align/bam/` |
| HISAT2 index | `results/rna-align/index/hcy` |
| Protein FASTA | `results/protein-qc/proteins.clean.faa` |
| Protein manifest | `results/protein-qc/proteins.clean.manifest.tsv` |
| EviAnn output | `results/eviann/hcy/` |
| Pipeline scripts | `pipeline/` |
| SLURM logs | `pipeline/joblogs/` |

All relative paths from: `/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/`
