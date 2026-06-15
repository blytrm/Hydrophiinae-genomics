## NextDenovo Config Scripts

Submitting/running:
```bash
nextDenovo <in.cfg>
```

**_Hydrophis curtus_ (west)**
`hcurw.cfg`
```bash
[General]
job_type = slurm
submit = sbatch --cpus-per-task={cpu} --mem 250GB -o {out} -e {err} {script}
job_prefix = hcu_west
task = all
rewrite = yes
deltmp = yes
parallel_jobs = 20
input_type = raw
read_type = ont
input_fofn = hcu_west.fofn
workdir = /home/a1645424/al-biohub/billy_trim/honours/01_genome_assembly/hydrophis_curtus-west/contig/assembly/workdir_hcu-west

[correct_option]
read_cutoff = 5k
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
```
`hcurw.fofn`
```bash
/home/path/to/data/sra-data/hydrophis_curtus-SAMN23040034/nanopore_wgs/hcu_west.fastq.gz
```
---

**_Hydrophis curtus_ (east)**
`hcure.cfg`
```bash
[General]
job_type = slurm
submit = sbatch --cpus-per-task={cpu} --mem 250GB -o {out} -e {err} {script}
job_prefix = hcu
task = all
rewrite = yes
deltmp = yes
parallel_jobs = 20
input_type = raw
read_type = clr
input_fofn = hcu.fofn
workdir = /home/a1645424/al-biohub/billy_trim/honours/01_genome_assembly/hydrophis_curtus/contig/assembly/workdir_hcu

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
```
`hcure.fofn`
```bash
/home/path/to/data/sra-data/hydrophis_curtus-PRJNA616080/pacbio_wgs-SRX8032201/SRR11454605_subreads.fastq.gz
```
---

**_Hydrophis cyanocinctus_**
`hcy.cfg`
```bash
[General]
job_type = slurm
submit = sbatch --cpus-per-task={cpu} --mem 250GB -o {out} -e {err} {script}
job_prefix = hcy
task = all
rewrite = yes
deltmp = yes
parallel_jobs = 20
input_type = raw
read_type = clr
input_fofn = hcy.fofn
workdir = /home/a1645424/al-biohub/billy_trim/honours/01_genome_assembly/workdir_hcy

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
```
`hcy.fofn`
```bash
/home/path/to/data/sra-data/hydrophis_cyanocinctus-PRJNA573877/pacbio_wgs-SRX6956830/hydrophis_cyanocinctus-subreads.fastq.gz
```
---
