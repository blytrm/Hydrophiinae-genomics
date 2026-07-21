# Old/New Assembly Synteny Plot Pipeline 

- tBLASTn query is `annotation/v2r/hma_7tm.fa`; alignment orientation is **ref = old, query = new**.

### _Hydrophis cyanocinctus old vs new_

![Hydrophis synteny plots (old vs new assembly)](synteny_facet.png)

- Left = all chromosomes; right = ch2 + chz subset (500 bp alignment filter).
- Red band = new assembly, blue band = old assembly, orange ticks = V2R genes.
- Ribbons = nucmer synteny links

## Run order

```bash
cd /hpcfs/users/a1864358/sanders_lab/asm/files/synteny
bash prep-asms.sh        # 1. prep assemblies
sbatch nucmer-all.sh     # 2. whole-genome alignment
sbatch v2r-tblastn.sh    # 3. tBLASTn V2R overlay
bash plot-all.sh         # 4. plot
```

- Step 4 needs steps 1, 2 and 3 done.

---



## 1. `prep-asms.sh` — prep assemblies

- Keeps the top-N largest sequences of each new/old assembly.
- Renames them to `ch1..chN` (and `chz` lowercase for the sex chromosome).
- Writes per-species `new.fa`, `old.fa`, `new-ref.csv`, `old-ref.csv`, `ch_rename.tsv`.

```bash
#!/usr/bin/env bash
# Per-species rename + subset for synteny.
set -eu

SYNT=/hpcfs/users/a1864358/sanders_lab/asm/files/synteny
BIN=/hpcfs/users/a1864358/miniconda/envs/general/bin
PATH=$BIN:$PATH

# new (final): take top-N largest seqs, rename to ch1..chN.
rename_top_to_ch() {
  local infa=$1 outfa=$2 N=$3
  [ -f ${infa}.fai ] || samtools faidx $infa
  sort -k2 -n -r ${infa}.fai | head -$N | awk '{print $1"\tch"NR}' > ${outfa}.map.tsv
  awk '{print $1}' ${outfa}.map.tsv > ${outfa}.regions
  samtools faidx -r ${outfa}.regions $infa > ${outfa}.tmp
  awk -v map=${outfa}.map.tsv '
    BEGIN { while ((getline l<map)>0) { split(l,a,"\t"); m[a[1]]=a[2] } }
    /^>/  { n=substr($0,2); sub(/ .*/,"",n); if (n in m) print ">"m[n]; else print $0; next }
    { print }
  ' ${outfa}.tmp > $outfa
  rm ${outfa}.tmp ${outfa}.regions
  samtools faidx $outfa
}

# old (orig): chr1..chrN + chrZ -> ch1..chN + chz, dropping empty seqs first.
rename_chr_to_ch() {
  local infa=$1 outfa=$2 N=$3
  seqkit seq -m 1 $infa > ${outfa}.cln 2>/dev/null
  samtools faidx ${outfa}.cln
  sort -k2 -n -r ${outfa}.cln.fai | head -$N | awk '
    { n=$1; if (n=="chrZ") nn="chz"; else { sub(/^chr/,"ch",n); nn=n } print $1"\t"nn }
  ' > ${outfa}.map.tsv
  awk '{print $1}' ${outfa}.map.tsv > ${outfa}.regions
  samtools faidx -r ${outfa}.regions ${outfa}.cln > ${outfa}.tmp
  awk -v map=${outfa}.map.tsv '
    BEGIN { while ((getline l<map)>0) { split(l,a,"\t"); m[a[1]]=a[2] } }
    /^>/  { n=substr($0,2); sub(/ .*/,"",n); if (n in m) print ">"m[n]; else print $0; next }
    { print }
  ' ${outfa}.tmp > $outfa
  rm ${outfa}.tmp ${outfa}.regions ${outfa}.cln ${outfa}.cln.fai
  samtools faidx $outfa
}

# Write chromosome-size CSVs and an identity rename table per species.
write_refs() {
  local sp=$1
  awk 'BEGIN{print "chromosome,size"} {print $1","$2}' $SYNT/$sp/new.fa.fai > $SYNT/$sp/new-ref.csv
  awk 'BEGIN{print "chromosome,size"} {print $1","$2}' $SYNT/$sp/old.fa.fai > $SYNT/$sp/old-ref.csv
  printf "original\tnew\n" > $SYNT/$sp/ch_rename.tsv
}

# Run per species (new = top-8/9 final chrs, old = original chr set).
sp=hcurw; cd $SYNT/$sp
rename_top_to_ch /hpcfs/users/a1864358/sanders_lab/asm/files/final-asms/10-hcurw-final.fa new.fa 8
rename_chr_to_ch /hpcfs/users/a1864358/sanders_lab/asm/files/hcur-w/hydrophis_curtus-west.fa old.fa 16
write_refs $sp

sp=hmaj; cd $SYNT/$sp
rename_top_to_ch /hpcfs/users/a1864358/sanders_lab/asm/files/final-asms/10-hmaj-final.fa new.fa 8
rename_chr_to_ch /hpcfs/users/a1864358/sanders_lab/asm/files/hmaj/hydrophis_major.fa old.fa 16
write_refs $sp

sp=horn; cd $SYNT/$sp
rename_top_to_ch /hpcfs/users/a1864358/sanders_lab/asm/files/final-asms/10-horn-final.fa new.fa 8
rename_top_to_ch /hpcfs/users/a1864358/sanders_lab/asm/files/horn/horn-og.fa old.fa 20
write_refs $sp

sp=hcy; cd $SYNT/$sp   # hcy uses pre-made top9/top18 subsets
[ -f new.fa ] || { ln -sf $SYNT/hcy-synt-top9.fasta new.fa; ln -sf $SYNT/hcy-synt-top9.fasta.fai new.fa.fai; }
[ -f old.fa ] || { awk -v map=$SYNT/ch_rename.tsv '
    BEGIN { while ((getline l<map)>0) { split(l,a,"\t"); m[a[1]]=a[2] } }
    /^>/  { n=substr($0,2); sub(/ .*/,"",n); if (n in m) print ">"m[n]; else print ">"n; next }
    { print }' $SYNT/hcy-synt-org-top18.fasta > old.fa; samtools faidx old.fa; }
write_refs $sp
```

---



## 2. `nucmer-all.sh` (SLURM) —> whole-genome alignment

- Runs `nucmer old.fa new.fa` then `show-coords -lTH` per species.
- Converts coords to `algn.csv` with the fixed 11-column header.
- Filters to `filtered-500-algn.csv` (ref & query alignment length ≥ 500 bp).

```bash
#!/usr/bin/env bash
#SBATCH --job-name=nucmer-all
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=48
#SBATCH --time=24:00:00
#SBATCH --mem=128GB
#SBATCH -o %x_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com

set -eu
PATH=/hpcfs/users/a1864358/miniconda/envs/general/bin:$PATH

SYNT=/hpcfs/users/a1864358/sanders_lab/asm/files/synteny
SPECIES=(hcurw hmaj horn hcy)
HEADER="reference_start,reference_end,query_start,query_end,reference_alignment_length,query_alignment_length,percent_identity,reference_length,query_length,reference_chromosome,query_chromosome"

for sp in "${SPECIES[@]}"; do
  cd $SYNT/$sp
  # ref = old, query = new (must match plot-synteny-run.R orientation).
  nucmer --prefix algn old.fa new.fa -t ${SLURM_CPUS_PER_TASK}
  show-coords -lTH algn.delta > algn.coords
  (echo "${HEADER}" && awk '{$1=$1}1' OFS="," algn.coords) > algn.csv
  awk -F',' 'NR==1 || ($5>=500 && $6>=500)' algn.csv > filtered-500-algn.csv
done
```

---



## 3. `v2r-tblastn.sh` (SLURM) —> tBLASTn V2R overlay

- Builds a BLAST db from each `new.fa` / `old.fa`.
- Runs `tblastn -outfmt 6 -evalue 1e-10` with the V2R protein query.
- Converts HSP coords to `v2r.new.bed` / `v2r.old.bed` (feeds per-chromosome V2R counts).

```bash
#!/usr/bin/env bash
#SBATCH --job-name=v2r-tblastn
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --mem=32GB
#SBATCH -o %x_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=13billy.trim13@gmail.com

set -eu
PATH=/hpcfs/users/a1864358/miniconda/envs/general/bin:$PATH

SYNT=/hpcfs/users/a1864358/sanders_lab/asm/files/synteny
QUERY=/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/v2r/hma_7tm.fa
SPECIES=(hcurw hmaj horn hcy)
THREADS=4

# tBLASTn each new/old assembly in parallel.
for sp in "${SPECIES[@]}"; do
  for which in new old; do
    DB=$SYNT/$sp/${which}.fa
    OUT=$SYNT/$sp/v2r.${which}.tblastn.outfmt6
    [ -f ${DB}.nhr ] || makeblastdb -input_type fasta -in $DB -dbtype nucl -parse_seqids
    tblastn -query $QUERY -db $DB -out $OUT -outfmt 6 -num_threads $THREADS -evalue 1e-10 &
  done
done
wait

# Convert hits to BED3 (chrom, start-1, end), coordinates sorted.
for sp in "${SPECIES[@]}"; do
  for which in new old; do
    OUT=$SYNT/$sp/v2r.${which}.tblastn.outfmt6
    BED=$SYNT/$sp/v2r.${which}.bed
    awk 'BEGIN{OFS="\t"} { s=$9; e=$10; if (s>e){t=s;s=e;e=t} print $2, s-1, e }' $OUT \
      | sort -k1,1 -k2,2n > $BED
  done
done
```

---



## 4. `plot-all.sh` → `plot-synteny-run.R` —> plot

- Checks all required inputs exist per species before plotting.
- Runs `Rscript plot-synteny-run.R <synteny/sp> 500`.
- Outputs `synteny_plot_500.png` (full) and `synteny_plot_500_chr2_Z.png` (ch2+chz subset).

```bash
#!/usr/bin/env bash
set -eu
SYNT=/hpcfs/users/a1864358/sanders_lab/asm/files/synteny
R=/hpcfs/users/a1864358/miniconda/envs/general/bin/Rscript
TAG=500

for sp in hcurw hmaj horn hcy; do
  WD=$SYNT/$sp
  # Skip species missing any required input.
  for f in new.fa.fai old.fa.fai new-ref.csv old-ref.csv ch_rename.tsv v2r.new.bed v2r.old.bed filtered-${TAG}-algn.csv; do
    if [ ! -e $WD/$f ]; then echo "[!] missing $WD/$f"; continue 2; fi
  done
  $R $SYNT/plot-synteny-run.R $WD $TAG
done
```

- `plot-synteny-run.R` draws the circlize plot (ribbons = nucmer links, palette = `ltc` "heatmap3", V2R beds merged within 5 kb give the `(n=K)` sector labels).

```r
# Wrapper around plot-synteny.R: parametrised inputs/outputs, bug fixes.
# Usage: Rscript plot-synteny-run.R <workdir> <tag>
#   workdir contains: new-ref.csv, old-ref.csv, ch_rename.tsv,
#                     v2r.new.bed, v2r.old.bed, filtered-<tag>-algn.csv
#   produces synteny_plot_<tag>.png, synteny_plot_<tag>_redribbons.png,
#   and chr2/chz subset variants in workdir.

suppressPackageStartupMessages({
  library(circlize)
  library(dplyr)
  library(ltc)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript plot-synteny-run.R <workdir> <tag>")
workdir <- args[1]
tag     <- args[2]

base_dir   <- workdir
algn_path  <- file.path(workdir, sprintf("filtered-%s-algn.csv", tag))
rename_tsv <- file.path(workdir, "ch_rename.tsv")
v2r_new_path <- file.path(workdir, "v2r.new.bed")
v2r_old_path <- file.path(workdir, "v2r.old.bed")

merge_bed <- function(df, gap = 0) {
  df_unique <- df %>% distinct(chr, start, end) %>% arrange(chr, start, end)
  if (nrow(df_unique) == 0) {
    return(data.frame(chr = character(), start = numeric(), end = numeric(),
                      stringsAsFactors = FALSE))
  }
  result <- data.frame(chr = character(), start = numeric(), end = numeric(),
                       stringsAsFactors = FALSE)
  current_chr   <- df_unique$chr[1]
  current_start <- df_unique$start[1]
  current_end   <- df_unique$end[1]
  if (nrow(df_unique) >= 2) {
    for (i in 2:nrow(df_unique)) {
      row <- df_unique[i, ]
      if (row$chr == current_chr && row$start <= current_end + gap) {
        current_end <- max(current_end, row$end)
      } else {
        result <- rbind(result, data.frame(chr = current_chr,
                                           start = current_start,
                                           end = current_end,
                                           stringsAsFactors = FALSE))
        current_chr   <- row$chr
        current_start <- row$start
        current_end   <- row$end
      }
    }
  }
  rbind(result, data.frame(chr = current_chr, start = current_start,
                           end = current_end, stringsAsFactors = FALSE))
}

check_v2r_overlap <- function(chr, start, end, v2r_data) {
  chr_v2r <- v2r_data %>% filter(chr == !!chr)
  if (nrow(chr_v2r) == 0) return(FALSE)
  any((start <= chr_v2r$end) & (end >= chr_v2r$start))
}

# Refs 
old_ref <- read.csv(file.path(base_dir, "old-ref.csv"),
                    stringsAsFactors = FALSE, strip.white = TRUE)
new_ref <- read.csv(file.path(base_dir, "new-ref.csv"),
                    stringsAsFactors = FALSE, strip.white = TRUE)
names(old_ref) <- trimws(names(old_ref))
names(new_ref) <- trimws(names(new_ref))
old_ref$assembly <- "old"
new_ref$assembly <- "new"
ref_data <- rbind(old_ref, new_ref)

old_chr <- ref_data %>% filter(assembly == "old") %>%
  select(chromosome, size) %>% mutate(chr = paste0("old_", chromosome))
new_chr <- ref_data %>% filter(assembly == "new") %>%
  select(chromosome, size) %>% mutate(chr = paste0("new_", chromosome))

all_chr <- rbind(
  data.frame(chr = old_chr$chr, start = 0, end = old_chr$size),
  data.frame(chr = new_chr$chr, start = 0, end = new_chr$size)
)
chr_sizes <- setNames(all_chr$end, all_chr$chr)

# Alignment 
cat(sprintf("[%s] reading alignment %s\n", tag, algn_path))
alignment <- read.csv(algn_path, stringsAsFactors = FALSE)

links <- alignment %>%
  mutate(chr1 = paste0("old_", reference_chromosome),
         start1 = reference_start, end1 = reference_end,
         chr2 = paste0("new_", query_chromosome),
         start2 = query_start, end2 = query_end,
         identity = percent_identity) %>%
  select(chr1, start1, end1, chr2, start2, end2, identity)

# Order start<end (nucmer can flip on reverse strand)
links <- links %>%
  mutate(s1 = pmin(start1, end1), e1 = pmax(start1, end1),
         s2 = pmin(start2, end2), e2 = pmax(start2, end2)) %>%
  mutate(start1 = s1, end1 = e1, start2 = s2, end2 = e2) %>%
  select(-s1, -e1, -s2, -e2)

links <- links %>%
  filter(identity >= 90) %>%
  filter((end1 - start1) >= 1000 | (end2 - start2) >= 1000) %>%
  filter(chr1 %in% names(chr_sizes) & chr2 %in% names(chr_sizes)) %>%
  mutate(chr1_size = chr_sizes[chr1], chr2_size = chr_sizes[chr2]) %>%
  filter(start1 >= 0 & end1 <= chr1_size &
         start2 >= 0 & end2 <= chr2_size &
         start1 < end1 & start2 < end2) %>%
  select(-chr1_size, -chr2_size)

cat(sprintf("[%s] links after filter: %d\n", tag, nrow(links)))
if (nrow(links) > 200000) {
  links <- links %>% arrange(desc(identity)) %>% slice_head(n = 200000)
  cat(sprintf("[%s] capped to 200000 by identity\n", tag))
}

# V2R bed (tBLASTn HSPs; merged below to count loci) 
v2r_new_raw <- tryCatch(read.table(v2r_new_path, sep = "\t", header = FALSE,
                                   stringsAsFactors = FALSE,
                                   col.names = c("chr", "start", "end")),
                        error = function(e) data.frame(chr = character(),
                                                       start = integer(),
                                                       end = integer()))
v2r_old_raw <- tryCatch(read.table(v2r_old_path, sep = "\t", header = FALSE,
                                   stringsAsFactors = FALSE,
                                   col.names = c("chr", "start", "end")),
                        error = function(e) data.frame(chr = character(),
                                                       start = integer(),
                                                       end = integer()))

# Old bed uses CM accessions; rename to chN
if (file.exists(rename_tsv) && nrow(v2r_old_raw) > 0) {
  rn <- read.table(rename_tsv, sep = "\t", header = FALSE,
                   stringsAsFactors = FALSE, col.names = c("from", "to"))
  map <- setNames(rn$to, rn$from)
  v2r_old_raw$chr <- ifelse(v2r_old_raw$chr %in% names(map),
                            map[v2r_old_raw$chr], v2r_old_raw$chr)
}

# 5kb gap-merge collapses HSP fragments from same V2R locus into one row.
v2r_new <- merge_bed(v2r_new_raw, gap = 5000)
v2r_old <- merge_bed(v2r_old_raw, gap = 5000)

v2r_new_plot <- v2r_new %>% mutate(chr = paste0("new_", chr)) %>%
  filter(chr %in% all_chr$chr) %>%
  mutate(chr_size = chr_sizes[chr]) %>%
  filter(start >= 0 & end <= chr_size & start < end) %>%
  select(chr, start, end)
v2r_old_plot <- v2r_old %>% mutate(chr = paste0("old_", chr)) %>%
  filter(chr %in% all_chr$chr) %>%
  mutate(chr_size = chr_sizes[chr]) %>%
  filter(start >= 0 & end <= chr_size & start < end) %>%
  select(chr, start, end)
v2r_all <- rbind(v2r_new_plot, v2r_old_plot)

# V2R loci per sector (used in chr labels)
v2r_counts <- v2r_all %>% count(chr, name = "n")
v2r_count_map <- setNames(v2r_counts$n, v2r_counts$chr)
v2r_count_for <- function(sector) {
  n <- v2r_count_map[sector]
  if (is.na(n)) 0L else as.integer(n)
}
cat(sprintf("[%s] V2R loci per sector:\n", tag))
print(v2r_counts)

# Colours / V2R link tagging 
links <- links %>% mutate(chr1_name = gsub("^(old_|new_)", "", chr1),
                          chr2_name = gsub("^(old_|new_)", "", chr2))
all_chr_names <- sort(unique(c(links$chr1_name, links$chr2_name)))
chr_pair_colours <- setNames(
  ltc("heatmap3", length(all_chr_names), "continuous"),
  all_chr_names
)

links <- links %>%
  mutate(link_colour_base = ifelse(chr1_name %in% names(chr_pair_colours),
                                   chr_pair_colours[chr1_name], "#808080"),
         overlaps_v2r_old = mapply(function(ch, st, en) check_v2r_overlap(ch, st, en, v2r_old_plot),
                                   chr1, start1, end1),
         overlaps_v2r_new = mapply(function(ch, st, en) check_v2r_overlap(ch, st, en, v2r_new_plot),
                                   chr2, start2, end2),
         overlaps_v2r = overlaps_v2r_new | overlaps_v2r_old)

NEW_COL <- "#cd5c5c"
OLD_COL <- "#4169E1"
V2R_ORANGE <- "#FF8C00"
V2R_BLUE <- "#4169E1"

chr_label <- function(chr) gsub("old_|new_|chr|ch", "", chr)

draw_circos_plot <- function(out_path, all_chr_ord, n_new, n_old, plot_links,
                             plot_v2r, link_mode = c("default", "red_faint"),
                             v2r_col = V2R_ORANGE, label_cex = 0.6,
                             legend_cex = 0.6) {
  link_mode <- match.arg(link_mode)
  gap_vec <- c(if (n_new > 1) rep(2, n_new - 1), 10,
               if (n_old > 1) rep(2, n_old - 1), 10)
  chr_bar_colours <- c(rep(NEW_COL, n_new), rep(OLD_COL, n_old))

  png(out_path, width = 3000, height = 3000, res = 300)
  on.exit({
    circos.clear()
    dev.off()
  }, add = TRUE)

  circos.clear()
  circos.par(start.degree = 90, gap.degree = gap_vec,
             cell.padding = c(0.01, 0, 0.01, 0))
  circos.genomicInitialize(all_chr_ord, plotType = NULL)

  circos.track(ylim = c(0, 1),
               panel.fun = function(x, y) {
                 xlim <- CELL_META$xlim
                 circos.rect(xlim[1], 0, xlim[2], 1,
                             col = chr_bar_colours[CELL_META$sector.numeric.index],
                             border = NA)
                 circos.text(mean(xlim), 2.6, chr_label(CELL_META$sector.index),
                             cex = label_cex, facing = "bending.inside",
                             niceFacing = TRUE, font = 2, col = "#06233D")
               }, bg.border = NA, track.height = 0.06)

  if (nrow(plot_v2r) > 0) {
    v2r_v <- plot_v2r %>% mutate(value = 0.5)
    circos.genomicTrack(v2r_v, ylim = c(0, 1),
                        panel.fun = function(region, value, ...) {
                          circos.genomicRect(region, value, ytop = 1, ybottom = 0,
                                             col = v2r_col, border = v2r_col, ...)
                        }, track.height = 0.05, bg.border = NA)
  }

  region1 <- plot_links %>% select(chr = chr1, start = start1, end = end1)
  region2 <- plot_links %>% select(chr = chr2, start = start2, end = end2)
  if (link_mode == "red_faint") {
    link_colour_vec <- adjustcolor(NEW_COL, alpha.f = 0.12)
  } else {
    link_colour_vec <- ifelse(plot_links$overlaps_v2r,
                              adjustcolor(V2R_ORANGE, alpha.f = 0.85),
                              adjustcolor(plot_links$link_colour_base, alpha.f = 0.6))
  }
  circos.genomicLink(region1, region2, col = link_colour_vec, border = NA)

  if (link_mode == "red_faint") {
    legend("topright",
           legend = c("New assembly", "Old assembly", "Synteny", "V2R"),
           fill = c(NEW_COL, OLD_COL,
                    adjustcolor(NEW_COL, alpha.f = 0.12), V2R_BLUE),
           border = NA, bty = "n", cex = legend_cex, ncol = 2)
  } else {
    legend("topright",
           legend = c("New assembly", "Old assembly", "V2R genes", "V2R links"),
           fill = c(NEW_COL, OLD_COL, V2R_ORANGE, V2R_ORANGE),
           border = NA, bty = "n", cex = legend_cex, ncol = 2)
  }
}

all_chr_ord <- rbind(
  data.frame(chr = new_chr$chr, start = 0, end = new_chr$size),
  data.frame(chr = old_chr$chr, start = 0, end = old_chr$size)
)

# Main plots 
out_main <- file.path(base_dir, sprintf("synteny_plot_%s.png", tag))
cat(sprintf("[%s] writing %s\n", tag, out_main))
draw_circos_plot(out_main, all_chr_ord, nrow(new_chr), nrow(old_chr),
                 links, v2r_all, link_mode = "default", v2r_col = V2R_ORANGE)

out_red <- file.path(base_dir, sprintf("synteny_plot_%s_redribbons.png", tag))
cat(sprintf("[%s] writing %s\n", tag, out_red))
draw_circos_plot(out_red, all_chr_ord, nrow(new_chr), nrow(old_chr),
                 links, v2r_all, link_mode = "red_faint", v2r_col = V2R_BLUE)

# Subset plots (ch2 / chz) 
chr_subset <- c("ch2", "chz")
links_subset <- links %>% filter(chr1_name %in% chr_subset & chr2_name %in% chr_subset)
v2r_subset <- v2r_all %>% filter(gsub("old_|new_", "", chr) %in% chr_subset)

if (nrow(links_subset) > 0) {
  new_sub <- new_chr %>% filter(chromosome %in% chr_subset)
  old_sub <- old_chr %>% filter(chromosome %in% chr_subset)
  all_chr_subset <- rbind(
    data.frame(chr = new_sub$chr, start = 0, end = new_sub$size),
    data.frame(chr = old_sub$chr, start = 0, end = old_sub$size)
  )
  n_new_s <- nrow(new_sub)
  n_old_s <- nrow(old_sub)

  out_sub <- file.path(base_dir, sprintf("synteny_plot_%s_chr2_Z.png", tag))
  cat(sprintf("[%s] writing %s (%d links)\n", tag, out_sub, nrow(links_subset)))
  draw_circos_plot(out_sub, all_chr_subset, n_new_s, n_old_s,
                   links_subset, v2r_subset, link_mode = "default",
                   v2r_col = V2R_ORANGE, label_cex = 0.8, legend_cex = 0.8)

  out_sub_red <- file.path(base_dir, sprintf("synteny_plot_%s_chr2_Z_redribbons.png", tag))
  cat(sprintf("[%s] writing %s\n", tag, out_sub_red))
  draw_circos_plot(out_sub_red, all_chr_subset, n_new_s, n_old_s,
                   links_subset, v2r_subset, link_mode = "red_faint",
                   v2r_col = V2R_BLUE, label_cex = 0.8, legend_cex = 0.8)
} else {
  cat(sprintf("[%s] no ch2/chz links — subset plot skipped\n", tag))
}

cat(sprintf("[%s] done\n", tag))
```



---

