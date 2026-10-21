#!/usr/bin/env python3
"""Add validated family-C anchors (mouse + anole + corn snake) to the backbone,
keeping anchors in the INGROUP and TAS1R as the only outgroup.
"""
FAMC = "/scratchdata1/users/a1864358/sanders_lab/phylogenomics/v2r_tree/workspace/data/familyC"
A = f"{FAMC}/anchors"


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:], []
        else:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


old = read_fasta(f"{FAMC}/anchored.clean.faa")
mouse = list(read_fasta(f"{A}/mouse_vmn2r1.faa").values())[0].rstrip("*")
anole = list(read_fasta(f"{A}/anole_anchor.faa").values())[0].rstrip("*")
corn = list(read_fasta(f"{A}/corn_anchor.faa").values())[0].rstrip("*")

anchors = [
    (">ANCHOR_Vmn2r1_MOUSE", mouse),
    (">ANCHOR_V2R1_ANOLE", anole),
    (">ANCHOR_V2R1like_CORNSNAKE", corn),
]
backbone = [(">" + n, s) for n, s in old.items() if not n.startswith("ANCHOR")]

with open(f"{FAMC}/anchored2.clean.faa", "w") as fh:
    for name, seq in anchors + backbone:
        fh.write(name + "\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")

with open(f"{FAMC}/outgroup2.txt", "w") as fh:
    tas1r = [n for n in old if n.startswith("TAS1R")]
    fh.write(",".join(tas1r))

print(f"anchors={len(anchors)} backbone={len(backbone)} tas1r_outgroup={len(tas1r)}")
