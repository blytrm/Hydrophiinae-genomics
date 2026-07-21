#!/usr/bin/env python3
"""Visualise CRAQ assembly-QC reports across genomes (hcy, hcur-w, horn).

AQI = Assembly Quality Indicator, 0-100, higher = better.
  R-AQI (from Avg.CRE) : regional / structural accuracy
  S-AQI (from Avg.CSE) : small-scale (base/local) accuracy
Parenthetical value in the report = the AQI; bare value = raw error count.
"""

import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid")
OUT = "/scratchdata1/users/a1864358/sanders_lab/craq_plots"
import os; os.makedirs(OUT, exist_ok=True)

BASE = "/hpcfs/users/a1864358/sanders_lab/asm/files"
SHORT = {
    "hcy":    f"{BASE}/hcy/files/hcy.backup0.backup0/3-dnqc/craq.Report",
    "hcur-w": f"{BASE}/hcur-w/2-dnqc/craq/out_final.Report",
    "horn":   f"{BASE}/horn/craq-report.Report",
}
REGIONAL = {
    "hcy":    f"{BASE}/hcy/files/hcy.backup0.backup0/3-dnqc/craq-regional.report",
    "hcur-w": f"{BASE}/hcur-w/2-dnqc/craq/out_regional.Report",
}

paren = re.compile(r"\(([-\d.eE+]+)\)")

def parse_short(path):
    """Return (genome_summary_dict, per_contig_DataFrame)."""
    rows, genome = [], None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("Short Report") or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 7:
                continue
            try:
                cov   = float(f[1]); lowc = float(f[2])
                raqi  = float(paren.search(f[5]).group(1))
                saqi  = float(paren.search(f[6]).group(1))
            except (ValueError, AttributeError):
                continue
            rec = dict(chr=f[0], covered=cov, lowconf=lowc, raqi=raqi, saqi=saqi)
            if f[0] == "Genome":
                genome = rec
            else:
                rows.append(rec)
    return genome, pd.DataFrame(rows)

def parse_regional(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    return df.rename(columns={"AQI": "aqi"})

summary, contigs = {}, {}
for g, p in SHORT.items():
    summary[g], contigs[g] = parse_short(p)
    print(f"{g}: {len(contigs[g])} contigs, genome R-AQI={summary[g]['raqi']:.2f} S-AQI={summary[g]['saqi']:.2f}")

genomes = list(SHORT)
colors = dict(zip(genomes, sns.color_palette("Set2", len(genomes))))

# genome-level summary bars (R-AQI, S-AQI, covered, low-conf)
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
x = np.arange(len(genomes)); w = 0.38
r = [summary[g]["raqi"] for g in genomes]
s = [summary[g]["saqi"] for g in genomes]
ax[0].bar(x - w/2, r, w, label="R-AQI (structural)", color="#4878CF")
ax[0].bar(x + w/2, s, w, label="S-AQI (small-scale)", color="#EE854A")
for i, (rv, sv) in enumerate(zip(r, s)):
    ax[0].text(i - w/2, rv + 0.5, f"{rv:.1f}", ha="center", fontsize=8)
    ax[0].text(i + w/2, sv + 0.5, f"{sv:.1f}", ha="center", fontsize=8)
ax[0].set(title="Genome-level AQI", ylabel="AQI (0-100, higher=better)", ylim=(0, 105))
ax[0].set_xticks(x); ax[0].set_xticklabels(genomes); ax[0].legend(fontsize=8)

cov = [summary[g]["covered"] * 100 for g in genomes]
ax[1].bar(x, cov, color=[colors[g] for g in genomes])
for i, v in enumerate(cov):
    ax[1].text(i, v + 0.005, f"{v:.3f}%", ha="center", fontsize=8)
ax[1].set(title="Mapped/Covered rate", ylabel="% genome covered", ylim=(99.5, 100.02))
ax[1].set_xticks(x); ax[1].set_xticklabels(genomes)

low = [summary[g]["lowconf"] * 100 for g in genomes]
ax[2].bar(x, low, color=[colors[g] for g in genomes])
for i, v in enumerate(low):
    ax[2].text(i, v, f"{v:.2e}", ha="center", va="bottom", fontsize=8)
ax[2].set(title="Low-confidence rate", ylabel="% low-confidence")
ax[2].set_xticks(x); ax[2].set_xticklabels(genomes)
fig.suptitle("CRAQ assembly QC — genome summary", fontweight="bold")
fig.tight_layout(); fig.savefig(f"{OUT}/01_genome_summary.png", dpi=150)
plt.close(fig)

#  per-contig AQI distribution across genomes
long = pd.concat(
    [contigs[g].assign(genome=g) for g in genomes], ignore_index=True
)
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for j, (col, title) in enumerate([("raqi", "Per-contig R-AQI"),
                                  ("saqi", "Per-contig S-AQI")]):
    sns.violinplot(data=long, x="genome", y=col, ax=ax[j], hue="genome",
                   palette=colors, cut=0, legend=False, inner="quartile")
    sns.stripplot(data=long, x="genome", y=col, ax=ax[j], color="k",
                  size=2, alpha=0.25)
    ax[j].set(title=title, ylabel="AQI", xlabel="")
fig.suptitle("CRAQ — per-contig AQI distribution", fontweight="bold")
fig.tight_layout(); fig.savefig(f"{OUT}/02_per_contig_aqi.png", dpi=150)
plt.close(fig)

# regional AQI track distribution (genomes with regional data) 
reg = {g: parse_regional(p) for g, p in REGIONAL.items()}
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for g, d in reg.items():
    ax[0].hist(d["aqi"], bins=40, alpha=0.55, label=g, color=colors[g],
               density=True)
ax[0].set(title="Regional (windowed) AQI distribution",
          xlabel="window AQI", ylabel="density"); ax[0].legend()

# fraction of windows below quality thresholds
thr = [50, 70, 90, 95, 99]
xb = np.arange(len(thr)); wb = 0.38
for k, g in enumerate(reg):
    frac = [(reg[g]["aqi"] < t).mean() * 100 for t in thr]
    ax[1].bar(xb + (k - 0.5) * wb, frac, wb, label=g, color=colors[g])
ax[1].set(title="Windows below AQI threshold", ylabel="% of windows",
          xlabel="AQI threshold")
ax[1].set_xticks(xb); ax[1].set_xticklabels(thr); ax[1].legend()
fig.suptitle("CRAQ — regional quality (hcy, hcur-w)", fontweight="bold")
fig.tight_layout(); fig.savefig(f"{OUT}/03_regional_aqi.png", dpi=150)
plt.close(fig)

print("Saved plots to", OUT)
