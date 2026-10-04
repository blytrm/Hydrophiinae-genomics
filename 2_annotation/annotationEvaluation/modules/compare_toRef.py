#!/usr/bin/env python3
"""Overlay new vs reference gene-length distributions + stats table.
Usage: compare_toRef.py <new.gff3> <ref.gff3> <out_prefix>"""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

new_gff, ref_gff, prefix = sys.argv[1], sys.argv[2], sys.argv[3]
cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]


def gene_lengths(path, label):
    df = pd.read_csv(path, sep="\t", comment="#", names=cols)
    g = df[df["type"] == "gene"].copy()
    g["gene_length"] = g["end"] - g["start"] + 1
    g["set"] = label
    return g[["gene_length", "set"]]


new = gene_lengths(new_gff, "new")
ref = gene_lengths(ref_gff, "reference")
both = pd.concat([new, ref], ignore_index=True)

plt.figure(figsize=(10, 6))
sns.histplot(data=both, x="gene_length", hue="set", bins=50, log_scale=True,
             element="step", common_norm=False, stat="density")
plt.title("Gene Length Distribution: new vs reference", fontsize=12)
plt.xlabel("Gene Length (bp; log scale)", fontsize=10)
plt.tight_layout()
plt.savefig(f"{prefix}_genelength.png", dpi=150)

stats = both.groupby("set")["gene_length"].agg(["count", "min", "max", "mean", "median"])
stats.to_csv(f"{prefix}_stats.tsv", sep="\t")
