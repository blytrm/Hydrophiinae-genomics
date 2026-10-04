#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

gff_path, out_png = sys.argv[1], sys.argv[2]
cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
df = pd.read_csv(gff_path, sep="\t", comment="#", names=cols)

genes = df[df["type"] == "gene"].copy()
genes["gene_length"] = genes["end"] - genes["start"] + 1

plt.figure(figsize=(10, 6))
sns.histplot(data=genes, x="gene_length", bins=50, log_scale=True,
             color="darkred", edgecolor="black")
plt.title("Gene Length Distribution", fontsize=12)
plt.xlabel("Gene Length (bp; log scale)", fontsize=10)
plt.ylabel("Frequency", fontsize=10)
plt.tight_layout()
plt.savefig(out_png, dpi=150)
