#!/usr/bin/env Rscript
# Usage: plot_cpg.R <genome.fa> <out.png> <out.tsv>
suppressMessages({
  library(Biostrings); library(GenomicRanges); library(GenomeInfoDb)
  library(regioneR);   library(karyoploteR)
})

args    <- commandArgs(trailingOnly = TRUE)
fa      <- args[1]; out_png <- args[2]; out_tsv <- args[3]
WIN <- 200L; STEP <- 100L                 # island-detection window / step
DENS_WIN <- 1e6L                          # density smoothing window (1 Mb)
CHR_RE  <- "^ch([0-9]+|[ZWXY])$"          # assembled chromosomes only

RED <- "#C0392B"; BLACK <- "#000000"

seqs <- readDNAStringSet(fa)
names(seqs) <- sub("\\s.*", "", names(seqs))   # first token of header

# --- chromosomes only, ordered by number then Z/W -------------------------
chr_all <- names(seqs)[grepl(CHR_RE, names(seqs))]
if (!length(chr_all)) stop("no chromosome-named sequences (", CHR_RE, ") found in ", fa)
num <- suppressWarnings(as.integer(sub("^ch", "", chr_all)))     # NA for chZ/chW
ord <- order(is.na(num), num, chr_all)                          # numeric first, then Z/W
keepch <- chr_all[ord]

lens   <- width(seqs); names(lens) <- names(seqs)
custom <- toGRanges(data.frame(chr = keepch, start = 1, end = as.numeric(lens[keepch])))

# --- de-novo CpG-island detection (Views: one row per window, not per base) --
find_islands <- function(seq, chrom) {
  L <- length(seq)
  if (L < WIN) return(NULL)
  starts <- seq.int(1L, L - WIN + 1L, by = STEP)
  v      <- Views(seq, start = starts, width = WIN)
  cg     <- letterFrequency(v, c("C", "G"))
  gc_fr  <- rowSums(cg) / WIN
  exp    <- (cg[, "C"] * cg[, "G"]) / WIN
  obs    <- vcountPattern("CG", v)
  keep   <- gc_fr > 0.5 & exp > 0 & (obs / exp) > 0.6
  if (!any(keep)) return(NULL)
  GRanges(chrom, IRanges(starts[keep], width = WIN))
}

hits    <- Filter(Negate(is.null), Map(find_islands, as.list(seqs[keepch]), keepch))
islands <- if (length(hits)) reduce(do.call(c, unname(hits))) else GRanges()  # merge adjacent
if (length(islands))
    islands <- keepSeqlevels(islands, keepch, pruning.mode = "coarse")

# --- karyoplot ------------------------------------------------------------
pp <- getDefaultPlotParams(plot.type = 1)
pp$data1height     <- 300
pp$ideogramheight  <- 30
pp$topmargin       <- 60
pp$bottommargin    <- 40
pp$leftmargin      <- 0.14
pp$rightmargin     <- 0.06

png(out_png, width = 1800, height = 120 * length(keepch) + 300, res = 150)
kp <- plotKaryotype(genome = custom, chromosomes = "all", plot.type = 1,
                    plot.params = pp, cex = 1.1,
                    main = "CpG-island density across chromosomes")
kpAddBaseNumbers(kp, tick.dist = 50e6, tick.len = 6, cex = 0.55,
                 tick.col = BLACK, add.units = TRUE)

if (length(islands)) {
  # black rug of individual island positions at the base of the data panel
  kpPlotRegions(kp, data = islands, col = BLACK, border = NA,
                r0 = 0, r1 = 0.06, data.panel = 1, avoid.overlapping = FALSE)
  # red density curve above it
  kd <- kpPlotDensity(kp, data = islands, window.size = DENS_WIN,
                      col = RED, border = RED, r0 = 0.10, r1 = 1,
                      data.panel = 1)
  ymax <- kd$latest.plot$computed.values$max.density
  # axis on the right so it never collides with chromosome names on the left
  kpAxis(kp, ymin = 0, ymax = ymax, r0 = 0.10, r1 = 1, numticks = 3,
         side = 2, cex = 0.55, col = BLACK, data.panel = 1)
  kpAddLabels(kp, labels = "CpG / Mb", srt = 90, pos = 3, cex = 0.7,
              col = BLACK, data.panel = 1, label.margin = 0.055)
}
invisible(dev.off())

# --- table ----------------------------------------------------------------
tab <- as.data.frame(table(factor(as.character(seqnames(islands)), levels = keepch)))
colnames(tab) <- c("chromosome", "cpg_islands")
write.table(tab, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
