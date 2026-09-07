# AutoHiC Genome Assembly Workflow

-> modular SLURM scripts for Hi-C scaffold correction, using deep-learning-assisted curation (using Hi-C contact heatmaps). 

---

### 5–7. AutoHiC Error Detection & Curation Prep — `05`–`07`

| Step | Script | Function |
|------|--------|----------|
| 5/3 | `05-autohic-hic-realign.sh` | Realign Hi-C to HapHiC scaffolds |
| 6/3 | `06-autohic-morehic.sh` | Convert BAM → `.hic`; run AutoHiC `onehic.py` (deep learning error classifier); output corrected `.assembly` + FASTA |
| 7/3 | `07-autohic-pretext.sh` | Align Hi-C to corrected assembly; generate `.pretext` contact map |

AutoHiC uses a ResNet-based model trained on Hi-C contact maps to detect misassembly breakpoints without requiring a reference genome.

> Li *et al.* (2023) AutoHiC: a deep-learning method for automatic and accurate chromosome-scale genome assembly. *Brief. Bioinform.* **24**, bbad310.
