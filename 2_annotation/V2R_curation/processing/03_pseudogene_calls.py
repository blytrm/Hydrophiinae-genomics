#!/usr/bin/env python3
"""
functional / pseudogene / partial calls.

EviAnn only emits complete ATG..stop ORFs, so its own models are stop-free by
construction and cannot tell you a locus is a pseudogene. 

miniprot, which aligns a known-good V2R protein
through the locus regardless of what it finds, reports the in-frame stops and
frameshifts it had to align through. Those are the disruption signal.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from v2rlib import (CHOSEN_RULE, RULES, SPECIES, coverage_of, gene_status,
                    n_genes, read_dicts, results_dir)

HERE = Path(__file__).resolve().parent.parent


def counts(rows, rule, cov):
    """class -> (transcripts, genes) under one rule."""
    out = {}
    for r in rows:
        k = gene_status(r, rule, cov(r))
        tx, genes = out.get(k, (0, set()))
        out[k] = (tx + 1, genes | ({r["gene_id"]} - {"NA"}))
    return {k: (tx, len(g)) for k, (tx, g) in out.items()}


def cell(c, key):
    tx, genes = c.get(key, (0, 0))
    return f"{tx}/{genes}"


def main():
    species = sys.argv[1:] or SPECIES
    report = ["# 03 — functional / pseudogene / partial", ""]

    for sp in species:
        out = results_dir(sp, HERE)
        rows = list(read_dicts(out / "v2r_master.tsv"))
        v2r = [r for r in rows if r["is_v2r"] == "yes"]
        raw = lambda r: r["coverage_class"]

        report += [f"## {sp} — {len(v2r)} V2R transcripts, {n_genes(v2r)} genes", "",
                   "Transcripts/genes, **before** → **after** the step-02 repair.", "",
                   "| rule for `functional` | functional | pseudogene | partial |",
                   "|---|---|---|---|"]
        for rule, (desc, _) in RULES.items():
            before, after = counts(v2r, rule, raw), counts(v2r, rule, coverage_of)
            mark = " **←**" if rule == CHOSEN_RULE else ""
            report.append(
                f"| {rule} {desc}{mark} "
                f"| {cell(before,'functional')} → **{cell(after,'functional')}** "
                f"| {cell(before,'pseudogene')} → {cell(after,'pseudogene')} "
                f"| {cell(before,'partial')} → {cell(after,'partial')} |")

        # write the chosen call back onto every row of both tables
        header = [c for c in rows[0] if c != "gene_status"] + ["gene_status"]
        for r in rows:
            r["gene_status"] = gene_status(r) if r["is_v2r"] == "yes" else "NA"
        for name, keep in (("v2r_master.tsv", lambda r: True),
                           ("v2r_genes.tsv", lambda r: r["is_v2r"] == "yes")):
            with open(out / name, "w") as fh:
                fh.write("\t".join(header) + "\n")
                for r in rows:
                    if keep(r):
                        fh.write("\t".join(str(r.get(c, "NA")) for c in header) + "\n")

        final = counts(v2r, CHOSEN_RULE, coverage_of)
        report += ["", f"Written as `gene_status` using rule {CHOSEN_RULE}: "
                   + ", ".join(f"**{k} {tx} tx / {g} genes**"
                               for k, (tx, g) in sorted(final.items())), ""]

    text = "\n".join(report) + "\n"
    (HERE / "processing" / "03_status_table.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
