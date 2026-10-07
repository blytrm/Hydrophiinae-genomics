#!/usr/bin/env python3
"""v2r merged-locus count table (post bedtools merge). ch1-9 by size + scaffolds,
old vs new. Two totals: [scaffs]=all, [no scaffs]=ch1-9 only."""
import os
from collections import defaultdict

OUTDIR = "/hpcfs/users/a1864358/sanders_lab/asm/files/synteny/v2r-full"
SPECIES = ["hcurw", "hmaj", "horn", "hcy"]
WHICH = ["old", "new"]
CH = [f"ch{i}" for i in range(1, 10)]
ROWS = CH + ["scaffolds"]

def size_map(fai):
    seqs = []
    with open(fai) as f:
        for line in f:
            p = line.split("\t"); seqs.append((p[0], int(p[1])))
    seqs.sort(key=lambda x: -x[1])
    return {sid: (f"ch{i+1}" if i < 9 else "scaffolds") for i, (sid, _) in enumerate(seqs)}

counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
for sp in SPECIES:
    for w in WHICH:
        m = size_map(f"{OUTDIR}/{sp}.{w}.fai")
        with open(f"{OUTDIR}/{sp}.{w}.filt.merged.bed") as f:
            for line in f:
                sid = line.split("\t", 1)[0]
                counts[sp][w][m.get(sid, "scaffolds")] += 1

cols = [(sp, w) for sp in SPECIES for w in WHICH]
hdr = ["chromosome"] + [f"{sp}_{w}" for sp, w in cols]
lines = ["\t".join(hdr)]
for r in ROWS:
    lines.append("\t".join([r] + [str(counts[sp][w].get(r, 0)) for sp, w in cols]))
lines.append("\t".join(["TOTAL [scaffs]"] + [str(sum(counts[sp][w].values())) for sp, w in cols]))
lines.append("\t".join(["TOTAL [no scaffs]"] + [str(sum(counts[sp][w].get(r, 0) for r in CH)) for sp, w in cols]))

tsv = f"{OUTDIR}/v2r-loci-merged.tsv"
open(tsv, "w").write("\n".join(lines) + "\n")
md = ["| " + " | ".join(hdr) + " |", "|" + "|".join(["---"]*len(hdr)) + "|"]
for ln in lines[1:]:
    md.append("| " + " | ".join(ln.split("\t")) + " |")
open(f"{OUTDIR}/v2r-loci-merged.md", "w").write("\n".join(md) + "\n")
print(f"wrote {tsv}\nwrote {OUTDIR}/v2r-loci-merged.md\n")
print("\n".join(md))
