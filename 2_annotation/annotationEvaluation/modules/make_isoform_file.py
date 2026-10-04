#!/usr/bin/env python3
"""Build an OMArk --isoform_file from a protein FASTA.


Usage: make_isoform_file.py <proteins.fasta> <out.splice>
"""
import re
import sys
from collections import OrderedDict

ISO_RE = re.compile(r"-mRNA-\d+$")


def main(fasta: str, out: str) -> int:
    genes: "OrderedDict[str, list]" = OrderedDict()
    n_prot = 0
    with open(fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            pid = line[1:].split()[0]
            n_prot += 1
            gene = ISO_RE.sub("", pid)
            genes.setdefault(gene, []).append(pid)

    with open(out, "w") as oh:
        for prots in genes.values():
            oh.write(";".join(prots) + "\n")

    multi = sum(1 for p in genes.values() if len(p) > 1)
    sys.stderr.write(
        f"[make_isoform_file] {n_prot} proteins -> {len(genes)} genes "
        f"({multi} multi-isoform)\n"
    )
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: make_isoform_file.py <proteins.fasta> <out.splice>")
    sys.exit(main(sys.argv[1], sys.argv[2]))
