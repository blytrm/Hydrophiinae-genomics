#!/usr/bin/env python3
"""Extract spliced CDS from an EviAnn transcripts FASTA.
"""
import re
import sys

CDS_RE = re.compile(r"CDS=(\d+)-(\d+)")


def iter_fasta(path):
    hdr, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(seq)
                hdr, seq = line[1:], []
            else:
                seq.append(line)
    if hdr is not None:
        yield hdr, "".join(seq)


def main(src: str, out: str) -> int:
    n_in = n_out = n_skip = 0
    with open(out, "w") as oh:
        for hdr, seq in iter_fasta(src):
            n_in += 1
            m = CDS_RE.search(hdr)
            if not m:
                n_skip += 1
                continue
            a, b = int(m.group(1)), int(m.group(2))
            if a < 1 or b > len(seq) or a > b:
                n_skip += 1
                continue
            pid = hdr.split()[0]
            oh.write(f">{pid}\n{seq[a - 1:b]}\n")
            n_out += 1
    sys.stderr.write(
        f"[extract_cds] {n_in} transcripts -> {n_out} CDS ({n_skip} skipped)\n"
    )
    return 0 if n_out else 1


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: extract_cds.py <transcripts.fasta> <out_cds.fasta>")
    sys.exit(main(sys.argv[1], sys.argv[2]))
