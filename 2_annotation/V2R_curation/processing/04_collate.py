#!/usr/bin/env python3
"""
collate every species into one comparison.
"""

import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from v2rlib import (RULES, SPECIES, coverage_of, gene_status, n_genes,
                    read_dicts, results_dir)

HERE = Path(__file__).resolve().parent.parent
FAI = "/hpcfs/users/a1864358/sanders_lab/asm/files/final-asms/10-{}-final.renamed.fa.fai"


def tally(rows, key):
    c = defaultdict(int)
    for r in rows:
        c[r[key]] += 1
    return c


data = {}
for sp in (sys.argv[1:] or SPECIES):
    path = results_dir(sp, HERE) / "v2r_master.tsv"
    if not path.exists():
        print(f"skipping {sp}: no {path}", file=sys.stderr)
        continue
    rows = list(read_dicts(path))
    v2r = [r for r in rows if r["is_v2r"] == "yes"]
    loci = {r["locus_id"] for r in rows if r["locus_id"] != "UNPLACED"}
    per_locus = {r["locus_id"]: r for r in rows}
    repaired = "repaired_coverage_class" in rows[0]
    cov_raw = lambda r: r["coverage_class"]
    cov_rep = coverage_of if repaired else cov_raw

    asm = {}
    try:
        for line in open(FAI.format(sp)):
            f = line.split("\t")
            asm[f[0]] = int(f[1])
    except OSError:
        pass

    d = dict(
        species=sp,
        assembly_mb=round(sum(asm.values()) / 1e6, 1) if asm else "-",
        n_chrom=sum(1 for c in asm if c.startswith("ch")),
        n_mrna_annotated=sum(1 for _ in open(results_dir(sp, HERE) / "eviann_mrna.bed")),
        loci=len(loci),
        v2r_tx=len(v2r), v2r_genes=n_genes(v2r),
        v2r_loci=len({r["locus_id"] for r in v2r if r["locus_id"] != "UNPLACED"}),
        rescued=sum(1 for r in v2r if r["v2r_evidence"] == "eviann_rescue"),
        tracks=",".join(sorted({t for r in rows for t in r["tracks"].split(",") if t})),
    )
    tiers = tally(list(per_locus.values()), "locus_tier")
    d.update({f"tier_{t}": tiers.get(t, 0) for t in ("T2", "T1", "T0")})

    for tag, cov in (("raw", cov_raw), ("rep", cov_rep)):
        c = defaultdict(int)
        for r in v2r:
            c[cov(r)] += 1
        for k in ("intact", "partial", "fragmented"):
            d[f"{k}_{tag}"] = c.get(k, 0)

    src = tally(v2r, "best_model_source") if repaired else {}
    for k in ("eviann", "miniprot", "orf"):
        d[f"src_{k}"] = src.get(k, "-" if not repaired else 0)
    d["merges_accepted"] = sum(1 for r in per_locus.values()
                               if r.get("merge_status") == "accepted") if repaired else "-"

    for rule in RULES:
        c = defaultdict(lambda: [0, set()])
        for r in v2r:
            k = gene_status(r, rule, cov_rep(r))
            c[k][0] += 1
            if r["gene_id"] != "NA":
                c[k][1].add(r["gene_id"])
        for k in ("functional", "pseudogene", "partial"):
            d[f"{rule}_{k}_tx"] = c[k][0]
            d[f"{rule}_{k}_genes"] = len(c[k][1])

    hc = [r for r in v2r if r["has_7tm"] == "yes" and int(r["n_exons"]) >= 6
          and 700 <= int(r["protein_len_aa"]) <= 1100]
    d["high_conf_tx"], d["high_conf_genes"] = len(hc), n_genes(hc)
    d["has_7tm"] = sum(1 for r in v2r if r["has_7tm"] == "yes")
    d["arch_inconsistent"] = sum(1 for r in v2r if r["eggnog_arch_consistent"] == "no")
    d["on_chrom_genes"] = n_genes([r for r in v2r if r["seq_type"] == "chromosome"])
    d["on_scaffold_genes"] = n_genes([r for r in v2r if r["seq_type"] == "scaffold"])
    d["inframe_stops"] = sum(int(r["inframe_stops"] or 0) for r in per_locus.values())
    d["frameshifts"] = sum(int(r["frameshifts"] or 0) for r in per_locus.values())
    ev = tally(v2r, "eviann_evidence")
    d["evidence_complete"] = ev.get("complete", 0)
    d["evidence_protein_only"] = ev.get("protein_only", 0)
    data[sp] = d

cols = list(next(iter(data.values())).keys())
with open(HERE / "processing" / "04_cross_species.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for sp in data:
        fh.write("\t".join(str(data[sp][c]) for c in cols) + "\n")


def table(title, rows_spec):
    out = [f"### {title}", "", "| | " + " | ".join(data) + " |",
           "|---|" + "---|" * len(data)]
    for label, key in rows_spec:
        out.append(f"| {label} | " + " | ".join(str(data[sp].get(key, "-")) for sp in data) + " |")
    return out + [""]


md = ["# 04 — cross-species collation", "",
      f"Species: {', '.join(data)}. All numbers computed identically.", ""]
md += table("Assembly and annotation", [
    ("assembly size (Mb)", "assembly_mb"), ("chromosomes", "n_chrom"),
    ("annotated mRNAs", "n_mrna_annotated"), ("evidence tracks used", "tracks")])
md += table("Candidate loci and tiers", [
    ("candidate loci", "loci"), ("T2", "tier_T2"), ("T1", "tier_T1"), ("T0", "tier_T0")])
md += table("V2R set", [
    ("transcripts", "v2r_tx"), ("genes", "v2r_genes"), ("loci", "v2r_loci"),
    ("EviAnn-rescued", "rescued"), ("7TM present", "has_7tm"),
    ("eggNOG architecture inconsistent", "arch_inconsistent")])
md += table("Completeness before repair", [
    ("intact", "intact_raw"), ("partial", "partial_raw"), ("fragmented", "fragmented_raw")])
md += table("Completeness after repair", [
    ("intact", "intact_rep"), ("partial", "partial_rep"), ("fragmented", "fragmented_rep"),
    ("winning model: eviann", "src_eviann"), ("winning model: miniprot", "src_miniprot"),
    ("winning model: orf", "src_orf"), ("merges accepted", "merges_accepted")])
md += table("Disruption evidence (miniprot)", [
    ("in-frame stops", "inframe_stops"), ("frameshifts", "frameshifts")])
md += table("Placement", [
    ("genes on chromosomes", "on_chrom_genes"), ("genes on scaffolds", "on_scaffold_genes")])
md += table("High confidence (7TM + >=6 exons + 700-1100 aa)", [
    ("transcripts", "high_conf_tx"), ("genes", "high_conf_genes")])

for rule, (desc, _) in RULES.items():
    md += table(f"Rule {rule} — {desc}", [
        ("functional (tx)", f"{rule}_functional_tx"),
        ("functional (genes)", f"{rule}_functional_genes"),
        ("pseudogene (tx)", f"{rule}_pseudogene_tx"),
        ("pseudogene (genes)", f"{rule}_pseudogene_genes"),
        ("partial (tx)", f"{rule}_partial_tx"),
        ("partial (genes)", f"{rule}_partial_genes")])

(HERE / "processing" / "04_cross_species.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
