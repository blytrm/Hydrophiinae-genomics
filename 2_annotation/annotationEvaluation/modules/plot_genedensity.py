#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

gff_path, out_png = sys.argv[1], sys.argv[2]
cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
df = pd.read_csv(gff_path, sep="\t", comment="#", names=cols)

counts = df[df["type"] == "gene"]["seqid"].value_counts().head(30)

plt.figure(figsize=(12, 6))
counts.plot(kind="bar", color="steelblue", edgecolor="black")
plt.title("Gene Count per Sequence (top 30)", fontsize=12)
plt.xlabel("Sequence", fontsize=10)
plt.ylabel("Gene count", fontsize=10)
plt.xticks(rotation=90, fontsize=7)
plt.tight_layout()
plt.savefig(out_png, dpi=150)
