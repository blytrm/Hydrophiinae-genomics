#!/usr/bin/env python3
"""One protein per gene, longest isoform wins.

    python3 scripts/longest_isoform.py <proteome.fa> <out.fa>

EviAnn's *.proteins.fasta has every isoform (LOC_00004150-mRNA-1, -mRNA-2, ...)
-- OrthoFinder wants one representative per gene, or splice variants inflate
orthogroups. Gene ID is just the header up to "-mRNA-", no GFF needed.
"""

import sys


def read_fasta(path):
    seqs, hdr = {}, None
    for line in open(path):
        if line.startswith(">"):
            hdr = line[1:].split()[0]
            seqs[hdr] = []
        elif hdr:
            seqs[hdr].append(line.strip())
    return {h: "".join(v) for h, v in seqs.items()}


def gene_of(header):
    return header.rsplit("-mRNA-", 1)[0] if "-mRNA-" in header else header


def main():
    in_path, out_path = sys.argv[1], sys.argv[2]
    seqs = read_fasta(in_path)

    best = {}
    for h, seq in seqs.items():
        g = gene_of(h)
        if g not in best or len(seq) > len(seqs[best[g]]):
            best[g] = h

    with open(out_path, "w") as fh:
        for g in sorted(best):
            h, seq = best[g], seqs[best[g]]
            fh.write(f">{h}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    print(f"{in_path}: {len(seqs)} isoforms, {len(best)} genes -> {out_path}")


if __name__ == "__main__":
    main()
