#!/usr/bin/env python3
"""One protein per gene (longest isoform), using an NCBI GFF for protein->gene.

    python3 longest_isoform_gff.py <protein.faa> <genomic.gff> <out.fa>

longest_isoform.py splits headers on "-mRNA-" (EviAnn). NCBI proteomes carry no
gene grouping in the header -- isoforms map to a gene only through the GFF:
    CDS   protein_id=XP_...  Parent=rna-XM_...
    mRNA  ID=rna-XM_...      Parent=gene-LOC...
 protein_id -> mRNA -> gene, then keep the longest protein per gene.
Proteins with no GFF gene (rare: no CDS record) fall back to their own id as gene.
"""

import re
import sys


def attr(field, key):
    m = re.search(key + r"=([^;]+)", field)
    return m.group(1) if m else None


def read_fasta(path):
    seqs, hdr = {}, None
    for line in open(path):
        if line.startswith(">"):
            hdr = line[1:].split()[0]
            seqs[hdr] = []
        elif hdr:
            seqs[hdr].append(line.strip())
    return {h: "".join(v) for h, v in seqs.items()}


def protein_to_gene(gff_path):
    rna_gene, prot_rna = {}, {}
    for line in open(gff_path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        if f[2] == "mRNA":
            rid, parent = attr(f[8], "ID"), attr(f[8], "Parent")
            if rid and parent:
                rna_gene[rid] = parent
        elif f[2] == "CDS":
            pid, parent = attr(f[8], "protein_id"), attr(f[8], "Parent")
            if pid and parent and pid not in prot_rna:
                prot_rna[pid] = parent
    return {pid: rna_gene.get(rid, pid) for pid, rid in prot_rna.items()}


def main():
    faa, gff, out = sys.argv[1], sys.argv[2], sys.argv[3]
    seqs = read_fasta(faa)
    p2g = protein_to_gene(gff)

    best = {}
    for h, seq in seqs.items():
        g = p2g.get(h, h)            # fall back to protein id if unmapped
        if g not in best or len(seq) > len(seqs[best[g]]):
            best[g] = h

    with open(out, "w") as fh:
        for g in sorted(best):
            h, seq = best[g], seqs[best[g]]
            fh.write(f">{h}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    print(f"{faa}: {len(seqs)} isoforms, {len(best)} genes -> {out}")


if __name__ == "__main__":
    main()
