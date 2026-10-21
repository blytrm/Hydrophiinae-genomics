#!/usr/bin/env python3
"""
Map each HyPhy site (trimmed codon-alignment column) to:
  (a) protein domain  VFT / linker / 7TM/9Cys  (majority vote across all reps)
  (b) reference residue, exon and genomic coordinate

Reads domain boundaries from split.tsv produced by build_domain_hmms.sh.
Emits selected_site_annotation.tsv to <out_dir>.


Args:
  codon.fasta      untrimmed codon alignment 
  codon.trim.fasta trimmed codon alignment 
  prot.aln         protein MSA 
  split.tsv        output of build_domain_hmms.sh -> split_domains.py
  ref_id           reference gene ID in prot.aln (e.g. hcure_V2R0001)
  ref.gff          miniprot GFF for the reference gene (CDS rows with Target=)
  parsed_dir       directory with meme_sites.tsv / fel_sites.tsv / cfel_sites.tsv
  out_dir          where to write selected_site_annotation.tsv
"""

import sys, os, re
from collections import Counter
from Bio import SeqIO

codon_fa, trim_fa, prot_fa, split_tsv, ref_id, gff, parsed_dir, out_dir = sys.argv[1:9]
os.makedirs(out_dir, exist_ok=True)

# --- recover kept codon columns (same >50%-gap rule as the trim step) ---
cod = list(SeqIO.parse(codon_fa, "fasta"))
n = len(cod)
ncodon = len(cod[0].seq) // 3
seqs = [str(r.seq) for r in cod]
kept = [c for c in range(ncodon)
        if sum(1 for s in seqs if "-" in s[c*3:c*3+3] or s[c*3:c*3+3] == "NNN") / n <= 0.5]
trim_n = len(next(SeqIO.parse(trim_fa, "fasta")).seq) // 3
site2orig = {i + 1: kept[i] for i in range(min(len(kept), trim_n))}

# --- domain boundaries per seq (from build_domain_hmms.sh -> split_domains.py) ---
dom = {}
with open(split_tsv) as f:
    f.readline()  # header
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if p[6] == "ok":
            dom[p[0]] = (int(p[3]), int(p[4]))  # (vft_end, core_start)

prot = {r.id: str(r.seq) for r in SeqIO.parse(prot_fa, "fasta")}

def residue_of(seq, col):
    """1-based ungapped residue index at alignment column col (None if gap)."""
    if seq[col] == "-":
        return None
    return sum(1 for a in seq[:col + 1] if a != "-")

# --- per-column consensus domain (majority vote across all reps) ---
col_dom = {}
for i, orig in site2orig.items():
    votes = Counter()
    for sid, seq in prot.items():
        if sid not in dom:
            continue
        res = residue_of(seq, orig)
        if res is None:
            continue
        vft_end, core_start = dom[sid]
        votes["VFT" if res <= vft_end else "7TM/9Cys" if res >= core_start else "linker"] += 1
    if votes:
        top, nvote = votes.most_common(1)[0]
        col_dom[i] = (top, nvote / sum(votes.values()))
    else:
        col_dom[i] = ("NA", 0)

# --- reference residue -> exon + genomic coordinate (miniprot GFF Target= ranges) ---
exons = []
for ln in open(gff):
    if "\tCDS\t" not in ln:
        continue
    f = ln.rstrip("\n").split("\t")
    m = re.search(r"Target=(\S+)\s+(\d+)\s+(\d+)", f[8])
    if m and m.group(1) == ref_id:
        exons.append((int(m.group(2)), int(m.group(3)), int(f[3]), int(f[4]), f[6]))
exons.sort()

def ref_exon_coord(res):
    for k, (ps, pe, gs, ge, strand) in enumerate(exons, 1):
        if ps <= res <= pe:
            off = res - ps
            gcoord = gs + off * 3 if strand == "+" else ge - off * 3
            return k, len(exons), gcoord, strand
    return None, len(exons), None, None

ref_seq = prot.get(ref_id)

# --- load selected site sets ---
def load_sites(tsv, filt=None):
    s = set()
    if not os.path.exists(tsv):
        return s
    with open(tsv) as f:
        hd = f.readline().rstrip("\n").split("\t")
        ci = hd.index("site")
        fi = hd.index(filt[0]) if filt else None
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if filt and p[fi] != filt[1]:
                continue
            s.add(int(p[ci]))
    return s

meme   = load_sites(f"{parsed_dir}/meme_sites.tsv")
felpos = load_sites(f"{parsed_dir}/fel_sites.tsv", filt=("class", "pos"))
cfel   = load_sites(f"{parsed_dir}/cfel_sites.tsv")

# --- write selected_site_annotation.tsv ---
out_tsv = f"{out_dir}/selected_site_annotation.tsv"
with open(out_tsv, "w") as o:
    o.write("hyphy_site\tdomain\tdomain_frac\tref_residue\tref_exon\tn_exons\tgenomic_coord\tstrand\tmeme\tfelpos\tcfel\n")
    for i in sorted(meme | felpos | cfel):
        d, fr = col_dom.get(i, ("NA", 0))
        res = residue_of(ref_seq, site2orig[i]) if ref_seq else None
        ex, ne, gc, st = ref_exon_coord(res) if res else (None, len(exons), None, None)
        o.write(f"{i}\t{d}\t{fr:.2f}\t{res or ''}\t{ex or ''}\t{ne}\t{gc or ''}\t{st or ''}\t"
                f"{int(i in meme)}\t{int(i in felpos)}\t{int(i in cfel)}\n")

# --- domain composition summary ---
bg = Counter(col_dom[i][0] for i in site2orig)
btot = sum(bg.values())
print("ALL columns: " + ", ".join(f"{k} {v} ({100*v/btot:.0f}%)" for k, v in bg.most_common()))
for name, sel in [("MEME episodic", meme), ("FEL positive", felpos), ("contrast-FEL", cfel)]:
    c = Counter(col_dom[i][0] for i in sel if i in col_dom)
    tot = sum(c.values())
    if tot:
        print(f"{name} (n={tot}): " + ", ".join(f"{k} {v} ({100*v/tot:.0f}%)" for k, v in c.most_common()))
print(f"reference: {ref_id}, {len(exons)} exons")
print(f"wrote {out_tsv}")
