#!/usr/bin/env bash
# Central config 

# Conda
CONDA_PROFILE="/hpcfs/users/a1864358/miniconda/miniconda3/etc/profile.d/conda.sh"
RNAQC_ENV="rnaqc"           # fastqc, bbduk, kraken2, fastp, multiqc, kraken2-build
ALIGN_ENV="rnaqc"           # hisat2, minimap2, samtools, ufasta  (TODO: change if separate)
EVIANN_ENV="eviann"         # eviann v2.0.4 installed via bioconda
HISAT2_MODULE=""   # hisat2 2.2.2 installed in rnaqc conda env

# Project layout
ANNO_ROOT="/hpcfs/users/a1864358/sanders_lab/asm/files/annotation"
DATA="${ANNO_ROOT}/data"
RNA_DIR="${DATA}/rna-seq"
PROT_DIR="${DATA}/prot-seq"
PIPE="${ANNO_ROOT}/pipeline"
FOFN_DIR="${PIPE}/fofn"

RESULTS="${ANNO_ROOT}/results"
OUT_RNAQC="${RESULTS}/rna-qc"
OUT_PROTQC="${RESULTS}/protein-qc"
OUT_ALIGN="${RESULTS}/rna-align"
OUT_EVIANN="${RESULTS}/eviann"

GENOME_LABEL="hcy"
GENOME="/hpcfs/users/a1864358/sanders_lab/asm/files/hcy/files/hcy.backup0.backup0/10-hcy-final.orignames.fa"
HISAT2_IDX="${OUT_ALIGN}/index/${GENOME_LABEL}"   # built by 04_align_rna.sh on first run

# Databases
KRAKEN2_DB="/scratchdata1/users/a1864358/dbs/kraken2/k2_db"
ADAPTER_REF="/hpcfs/users/a1864358/miniconda/envs/general/share/bbmap/resources/adapters.fa"
RRNA_REF="/scratchdata1/users/a1864358/dbs/rrna_refs/rrna_combined.fa"
UNIPROT_REF=""              # optional: path to uniprot_sprot.fasta for EviAnn -s flag

# QC parameters
MINLEN=50
FASTP_AVG_QUAL=15
FASTP_QUALIFIED_QUAL=15
FASTP_N_BASE_LIMIT=5

# Alignment parameters
MIN_INTRON=20
MAX_INTRON=1388389

# SLURM (override per-script with #SBATCH directives where useful)
SLURM_PARTITION="icelake"
SLURM_MAIL="13billy.trim13@gmail.com"

mkdir -p "${RESULTS}" "${OUT_RNAQC}" "${OUT_PROTQC}" "${OUT_ALIGN}" "${OUT_EVIANN}" "${FOFN_DIR}" "${PIPE}/joblogs"
