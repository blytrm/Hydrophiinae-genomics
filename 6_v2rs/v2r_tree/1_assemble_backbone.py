#!/usr/bin/env python3
"""Assemble the intact-only V2R backbone FASTA for the 12-taxon tree.

Reads each taxon's curated FASTA, keeps only class=intact records, re-headers
to <sptag>_V2R#### leaf name, concatenates, and writes the gene->species map for GeneRax
"""
import os
import sys

CUR = "/hpcfs/users/a1864358/sanders_lab/asm/files/annotation/eviann_results/hmaj/curation"
V2R12 = "/scratchdata1/users/a1864358/sanders_lab/phylogenomics/orthofinder/v2r_all12"
OUTD = "/scratchdata1/users/a1864358/sanders_lab/phylogenomics/v2r_tree/workspace/data"

# (sptag, source fasta) — headers of every source contain "class=intact" for intact genes
SOURCES = [
    ("hmaj",    f"{V2R12}/hmaj.fa"),
    ("hcy",     f"{V2R12}/hcy.fa"),
    ("hcure",   f"{V2R12}/hcure.fa"),
    ("hcurw",   f"{V2R12}/hcurw.fa"),
    ("horn",    f"{V2R12}/horn.fa"),
    ("nscut",   f"{CUR}/results/nscut_denovo/v2r_nscut_denovo.fasta"),
    ("ptext",   f"{CUR}/results/ptext_denovo/v2r_ptext_denovo.fasta"),
    ("nhel",    f"{CUR}/results/nhel/v2r_nhel.fasta"),
    ("tele",    f"{CUR}/results/tele/v2r_tele.fasta"),
    ("vberus",  f"{CUR}/results/vberus/v2r_vberus.fasta"),
    ("caspera", f"{CUR}/results/caspera/v2r_caspera.fasta"),
    ("vkomodo", f"{CUR}/results/vkomodo/v2r_vkomodo.fasta"),
]

os.makedirs(OUTD, exist_ok=True)
bb_path = os.path.join(OUTD, "backbone.faa")
map_path = os.path.join(OUTD, "gene_species_map.tsv")

counts = {}
total = 0
with open(bb_path, "w") as bb, open(map_path, "w") as mp:
    for tag, path in SOURCES:
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            sys.exit(f"ERROR: missing/empty source for {tag}: {path}")
        n = 0
        keep = False
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    keep = "class=intact" in line
                    if keep:
                        n += 1
                        leaf = f"{tag}_V2R{n:04d}"
                        bb.write(f">{leaf}\n")
                        mp.write(f"{leaf}\t{tag}\n")
                elif keep:
                    bb.write(line if line.endswith("\n") else line + "\n")
        counts[tag] = n
        total += n
        print(f"  {tag:<8} intact tips: {n}")

print(f"\nTOTAL backbone tips: {total}")
print(f"backbone: {bb_path}")
print(f"map:      {map_path}")
if total < 100:
    sys.exit(f"ERROR: backbone unexpectedly small ({total}); aborting")
