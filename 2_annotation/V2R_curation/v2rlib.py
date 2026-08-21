"""
Shared helpers and thresholds for the V2R curation.
"""

from collections import defaultdict
from pathlib import Path

SPECIES = ["hmaj", "hcy", "hcure", "hcurw", "horn"]

TARGET = "V2R_OlfC"       # the family being curated
MIN_DELTA = 20.0          # bitscore margin over the runner-up for a confident call
EVALUE_CEIL = 1e-3        # a family only competes on hits it actually supports
MIN_SINGLE_BITS = 100.0   # absolute score needed when no other family competes
COV_INTACT = 0.80         # fraction of the V2R profile covered -> intact
COV_PARTIAL = 0.60        # 0.60-0.80 -> partial; below -> fragmented


def results_dir(species, here=None):
    here = Path(here or Path(__file__).resolve().parent)
    return here / "results" / species


def read_tsv(path, header=True):
    """Yield each row as a list of fields."""
    with open(path) as fh:
        if header:
            next(fh, None)
        for line in fh:
            if line.strip():
                yield line.rstrip("\n").split("\t")


def read_dicts(path):
    """Yield each row as a dict keyed by the header line."""
    with open(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if line.strip():
                yield dict(zip(cols, line.rstrip("\n").split("\t")))


def write_tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def union_len(spans):
    """Length covered by a set of closed intervals, overlaps counted once."""
    spans = sorted(spans)
    total, cur_s, cur_e = 0, *spans[0]
    for s, e in spans[1:]:
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
    return total + cur_e - cur_s + 1


def read_domtbl(path, hmm_coords=False):
    """target -> (spans, model_len, best i-Evalue) from an hmmsearch domtblout.

    hmm_coords picks the model columns the domain covers (fields 16/17) instead
    of the protein columns (18/19). HMMER forbids whitespace in names, so a
    plain split is exact.
    """
    lo, hi = (15, 16) if hmm_coords else (17, 18)
    spans, model_len, evalue = defaultdict(list), {}, {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        spans[f[0]].append((int(f[lo]), int(f[hi])))
        model_len[f[0]] = int(f[5])
        e = float(f[12])
        if f[0] not in evalue or e < evalue[f[0]]:
            evalue[f[0]] = e
    return spans, model_len, evalue


def profile_coverage(path):
    """target -> fraction of the V2R model its domains cover."""
    spans, model_len, _ = read_domtbl(path, hmm_coords=True)
    return {t: union_len(v) / model_len[t] for t, v in spans.items()}


def coverage_class(cov):
    if cov is None:
        return "not_assessed"
    return ("intact" if cov >= COV_INTACT else
            "partial" if cov >= COV_PARTIAL else "fragmented")


# Functional / pseudogene / partial.
# `pseudogene` depends only on whether miniprot had to align through a broken reading frame
RULES = {
    "A": ("intact + 7TM + no disruption",
          lambda r, cls: cls == "intact" and r["has_7tm"] == "yes"),
    "B": ("intact + no disruption",
          lambda r, cls: cls == "intact"),
    "C": ("intact|partial + 7TM + no disruption",
          lambda r, cls: cls in ("intact", "partial") and r["has_7tm"] == "yes"),
    "D": ("intact|partial + no disruption",
          lambda r, cls: cls in ("intact", "partial")),
    "E": ("7TM + no disruption",
          lambda r, cls: r["has_7tm"] == "yes"),
    "F": ("intact + 7TM, ignores disruption",
          lambda r, cls: cls == "intact" and r["has_7tm"] == "yes"),
}
# Final spec: functional = repaired_coverage_class == intact AND inframe_stops == 0 AND frameshifts == 0.
CHOSEN_RULE = "B"


def disrupted(row):
    return int(row["inframe_stops"] or 0) > 0 or int(row["frameshifts"] or 0) > 0


def coverage_of(row):
    """The repaired class where step 02 has run, the raw one otherwise."""
    return row.get("repaired_coverage_class") or row["coverage_class"]


def gene_status(row, rule=CHOSEN_RULE, cls=None):
    """functional / pseudogene / partial for one transcript row."""
    cls = cls or coverage_of(row)
    ok = RULES[rule][1](row, cls)
    if rule == "F":                      # F deliberately ignores disruption
        return "functional" if ok else ("pseudogene" if disrupted(row) else "partial")
    if ok and not disrupted(row):
        return "functional"
    return "pseudogene" if disrupted(row) else "partial"


def count_by(rows, key):
    c = defaultdict(int)
    for r in rows:
        c[r[key] if isinstance(r, dict) else r] += 1
    return c


def n_genes(rows):
    return len({r["gene_id"] for r in rows if r["gene_id"] != "NA"})
