# De novo genome assembly
## Using [NextDenovo](https://github.com/Nextomics/NextDenovo)

#### Configuration
```bash
ls SRR16961044_1.fastq.gz SRR16961045_1.fastq.gz SRR16961047_1.fastq.gz SRR16961048_1.fastq.gz SRR17056026_1.fastq.gz SRR17056027_1.fastq.gz > heleg.fofn

cat > run.cfg << 'EOF'
[General]
job_type = slurm
submit = sbatch --cpus-per-task={cpu} --mem 250GB -o {out} -e {err} {script}
job_prefix = hor
task = all
rewrite = yes
deltmp = yes
parallel_jobs = 20
input_type = raw
read_type = clr
input_fofn = heleg.fofn
workdir = /home/a1645424/al-biohub/billy_trim/honours/01_genome_assembly/hydrophis_curtus/contig/assembly/workdir_heleg

[correct_option]
read_cutoff = 5k
seed_cutoff = 10k
genome_size = 2g
blocksize = 15g
sort_options = -m 100g -t 16
minimap2_options_raw = -t 16
pa_correction = 10
correction_options = -p 4

[assemble_option]
minimap2_options_cns = -t 10 -k17 -w17
minimap2_options_map = -t 10
nextgraph_options = -a 1
EOF
```

#### Running NextDenovo
```bash

#!/bin/bash
#SBATCH --job-name=nextdenovo_heleg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=2500G
#SBATCH --time=7-00:00:00
#SBATCH -p icelake
#SBATCH --output=nextdenovo_%j.out
#SBATCH --error=nextdenovo_%j.err
#SBATCH --mail-user=13billy.trim13@gmail.com
#SBATCH --mail-type=ALL

module purge

NEXTDENOVO_PATH="/hpcfs/users/a1864358/sanders_lab/asm/files/hcur-e/NextDenovo"
cd $NEXTDENOVO_PATH
export PATH=$NEXTDENOVO_PATH:$PATH

WORKDIR="/home/a1645424/al-biohub/billy_trim/honours/01_genome_assembly/hydrophis_curtus/contig/assembly/workdir_heleg"
cd $WORKDIR || mkdir -p $WORKDIR && cd $WORKDIR

which nextDenovo || { echo "ERROR: nextDenovo not found in PATH"; exit 1; }

nextDenovo run.cfg
```
