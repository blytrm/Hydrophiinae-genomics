#!/usr/bin/env python3
"""02 — gene-model repair.

45% of EviAnn V2R models cover under 60% of the 823-column V2R profile. Most of
those are not short receptors: they are one gene chopped into pieces by the
annotator. Two other sources of sequence models exist at the same loci:

    miniprot   spliced alignments of full-length V2R proteins, one per locus
    orf        getorf ORFs already scanned with the V2R profile

---> so for every locus we can ask which of the three models is the most complete,
using V2R-profile coverage. 
The pHMM itself is not a source: it scores proteins -> does not build models (scores)

A miniprot model that spans two or more EviAnn genes -> merged if suitable
"""

import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from v2rlib import COV_INTACT, COV_PARTIAL, coverage_class, profile_coverage, read_dicts

SP = sys.argv[1] if len(sys.argv) > 1 else "hmaj"
HERE = Path(__file__).resolve().parent.parent
OUT = HERE / "results" / SP
ROOT = Path("/scratchdata1/users/a1864358/sanders_lab/asm/files/annotation/eviann_results")
BASE = ROOT / SP
GENOME = Path(f"/hpcfs/users/a1864358/sanders_lab/asm/files/final-asms/10-{SP}-final.renamed.fa")
V2R_HMM = Path("/hpcfs/users/a1864358/sanders_lab/v2r_hmm/hmm/v2r_full.hmm")
WORK = HERE / "processing" / f"02_work_{SP}"
WORK.mkdir(parents=True, exist_ok=True)

GUARD_MIN_COV = 0.90         # guard 2: a merged model must be near-complete
GUARD_MIN_QUERY_FRAC = 0.80  # guard 1: the query protein must be spanned, not clipped


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)


def attr(s, k):
    m = re.search(k + r"=([^;]+)", s)
    return m.group(1) if m else None



# miniprot models: translate them, then score with the same V2R profile the  pipeline uses, so the coverage numbers are directly comparable.
mp_faa, mp_domtbl = WORK / "miniprot_prot.faa", WORK / "miniprot_v2r.domtbl"
if not mp_domtbl.exists():
    # miniprot writes ##PAF comment lines that contain literal tabs; gffread
    # parses them as feature rows and bails, so strip them first.
    clean = WORK / "miniprot.clean.gff"
    sh(f"grep -v '^##PAF' {BASE}/miniprot.v2r.gff > {clean}")
    sh(f'gffread -y {mp_faa} -g {GENOME} {clean} 2>/dev/null')
    sh(f'hmmsearch --cpu 16 -E 1e-3 --domtblout {mp_domtbl} {V2R_HMM} {mp_faa} > /dev/null')

mp_cov = profile_coverage(mp_domtbl)

# miniprot model coordinates, exon blocks, and how much of its query it used
mp_model, mp_cds = {}, defaultdict(list)
for line in open(BASE / "miniprot.v2r.gff"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9:
        continue
    tgt = attr(f[8], "Target")
    if f[2] == "mRNA" and tgt:
        q = tgt.split()
        mp_model[attr(f[8], "ID")] = dict(
            chrom=f[0], start=int(f[3]) - 1, end=int(f[4]), strand=f[6],
            query=q[0], qs=int(q[1]), qe=int(q[2]),
            stops=int(attr(f[8], "StopCodon") or 0),
            fs=int(attr(f[8], "Frameshift") or 0))
    elif f[2] == "CDS":
        mp_cds[attr(f[8], "Parent")].append((int(f[3]) - 1, int(f[4])))

# query protein lengths, for guard 1
qlen = {}
name = None
for line in open(ROOT / "hmaj" / "query_nr.faa"):
    if line.startswith(">"):
        name = line[1:].split()[0]
        qlen[name] = 0
    elif name:
        qlen[name] += len(line.strip())

# ORF models, if the getorf scan was run for this species. Coordinates are in the domtbl description
orf_cov, orf_pos = {}, {}
orf_domtbl = BASE / f"{SP}.orf.domtbl"
if orf_domtbl.exists():
    orf_cov = profile_coverage(orf_domtbl)
    for line in open(orf_domtbl):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split(maxsplit=22)
        m = re.search(r"\[\s*(\d+)\s*-\s*(\d+)\s*\]", f[22] if len(f) > 22 else "")
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            chrom = re.sub(r"_\d+$", "", f[0])
            orf_pos[f[0]] = (chrom, min(a, b) - 1, max(a, b))

# Which models sit in which locus.
def write_bed(path, rows):
    with open(path, "w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

loci_bed = WORK / "loci.bed"
sh(f"cut -f1-4 {OUT}/loci.bed | sort -k1,1 -k2,2n > {loci_bed}")

def locus_members(bed_rows, tag):
    src = WORK / f"{tag}.bed"
    write_bed(src, sorted(bed_rows, key=lambda r: (r[0], r[1])))
    sh(f"sort -k1,1 -k2,2n -o {src} {src}")
    out = subprocess.run(f"bedtools intersect -wa -wb -a {loci_bed} -b {src}",
                         shell=True, capture_output=True, text=True).stdout
    hit = defaultdict(list)
    for line in out.strip().split("\n"):
        if line:
            f = line.split("\t")
            hit[f[3]].append(f[7])
    return hit

mp_in_locus = locus_members(
    [(v["chrom"], v["start"], v["end"], k) for k, v in mp_model.items()], "miniprot")
orf_in_locus = locus_members(
    [(p[0], p[1], p[2], k) for k, p in orf_pos.items()], "orf") if orf_pos else {}

# The master table, and the EviAnn side of the comparison.
rows = list(read_dicts(OUT / "v2r_master.tsv"))
header = list(rows[0].keys())

ev_cov, ev_genes, ev_support, ev_span = defaultdict(float), defaultdict(set), {}, {}
for r in rows:
    if r["transcript_id"] == "NA":
        continue
    c = r["v2r_model_cov"]
    if c not in ("NA", ""):
        ev_cov[r["locus_id"]] = max(ev_cov[r["locus_id"]], float(c))
    if r["gene_id"] != "NA":
        ev_genes[r["locus_id"]].add(r["gene_id"])
        ev_support[r["gene_id"]] = r["eviann_evidence"]
        lo, hi = int(r["tx_start"]), int(r["tx_end"])
        s = ev_span.get(r["gene_id"])
        ev_span[r["gene_id"]] = (min(lo, s[0]), max(hi, s[1])) if s else (lo, hi)

# Per locus: pick the most complete model, then guard any merge it implies.
review, decisions = [], {}
for locus in {r["locus_id"] for r in rows if r["locus_id"] != "UNPLACED"}:
    e = ev_cov.get(locus, 0.0)
    best_mp = max((mp_cov.get(m, 0.0), m) for m in mp_in_locus.get(locus, [])) \
        if mp_in_locus.get(locus) else (0.0, None)
    best_orf = max((orf_cov.get(o, 0.0), o) for o in orf_in_locus.get(locus, [])) \
        if orf_in_locus.get(locus) else (0.0, None)

    source, cov, model = "eviann", e, None
    if best_mp[0] > cov:
        source, cov, model = "miniprot", best_mp[0], best_mp[1]
    if best_orf[0] > cov:
        source, cov, model = "orf", best_orf[0], best_orf[1]

    merge_status = "none"
    if source == "miniprot" and model:
        m = mp_model[model]
        # genes wholly or partly inside this miniprot model's span
        spanned = [g for g in ev_genes.get(locus, ())
                   if g in ev_span and ev_span[g][0] < m["end"] and ev_span[g][1] > m["start"]]
        if len(spanned) > 1:
            fails = []
            # guard 1 — one query protein spanning the whole thing, not a clipped alignment
            qfrac = (m["qe"] - m["qs"]) / max(1, qlen.get(m["query"], 1))
            if qfrac < GUARD_MIN_QUERY_FRAC:
                fails.append(f"g1_query_frac={qfrac:.2f}")
            # guard 2 — the repaired model must be near-complete, not merely better
            if cov < GUARD_MIN_COV:
                fails.append(f"g2_cov={cov:.2f}")
            # guard 3 — no spanned gene may have independent transcript evidence
            supported = [g for g in spanned
                         if ev_support.get(g) in ("complete", "transcript_only")]
            if supported:
                fails.append("g3_transcript_support=" + ",".join(sorted(supported)))
            # guard 4 — the merged model's exons must cover the fragments' spans
            cds = mp_cds.get(model, [])
            uncovered = [g for g in spanned
                         if not any(a < ev_span[g][1] and b > ev_span[g][0] for a, b in cds)]
            if uncovered:
                fails.append("g4_exons_missing=" + ",".join(sorted(uncovered)))

            merge_status = "accepted" if not fails else "rejected:" + ";".join(fails)
            if fails:                       # rejected merges keep the EviAnn model
                source, cov = "eviann", e
            review.append([locus, model, m["query"], f"{qfrac:.2f}", len(spanned),
                           ",".join(sorted(spanned)), f"{e:.3f}", f"{best_mp[0]:.3f}",
                           merge_status])

    cls = coverage_class(cov) if cov else "not_assessed"
    decisions[locus] = (source, f"{cov:.3f}", cls, merge_status)

# Write the columns back, and the review file.
new_cols = ["best_model_source", "repaired_model_cov", "repaired_coverage_class",
            "merge_status"]
for c in new_cols:
    if c in header:
        header.remove(c)
header += new_cols
with open(OUT / "v2r_master.tsv", "w") as fh:
    fh.write("\t".join(header) + "\n")
    for r in rows:
        d = decisions.get(r["locus_id"], ("NA", "NA", "NA", "none"))
        r.update(dict(zip(new_cols, d)))
        fh.write("\t".join(str(r.get(c, "NA")) for c in header) + "\n")

with open(HERE / "processing" / f"02_merge_review_{SP}.tsv", "w") as fh:
    fh.write("locus_id\tminiprot_model\tquery\tquery_frac\tn_genes_spanned\t"
             "genes\teviann_cov\tminiprot_cov\tmerge_status\n")
    for r in sorted(review):
        fh.write("\t".join(str(x) for x in r) + "\n")

# Report.
v2r = [r for r in rows if r["is_v2r"] == "yes"]
def tally(rs, key):
    c = defaultdict(int)
    for r in rs:
        c[r[key]] += 1
    return c

before, after = tally(v2r, "coverage_class"), tally(v2r, "repaired_coverage_class")
src = tally(v2r, "best_model_source")
acc = sum(1 for r in review if r[8] == "accepted")
rej = len(review) - acc

lines = [f"# 02 — gene-model repair: {SP}", "",
         "Coverage of the winning model per locus, EviAnn vs miniprot vs ORF,",
         "measured against the same 823-column V2R profile.", "",
         "| coverage class | before | after repair |", "|---|---|---|"]
for c in ("intact", "partial", "fragmented", "not_assessed"):
    if before[c] or after[c]:
        lines.append(f"| {c} | {before[c]} | {after[c]} |")
lines += ["", "| winning model source | transcripts |", "|---|---|"]
for k, v in sorted(src.items(), key=lambda kv: -kv[1]):
    lines.append(f"| {k} | {v} |")
lines += ["", f"Merges proposed: {len(review)} — **{acc} accepted**, {rej} rejected by a guard.",
          f"Detail in `02_merge_review_{SP}.tsv`.", ""]
if rej:
    gfail = defaultdict(int)
    for r in review:
        if r[8] != "accepted":
            for f in r[8].split(":", 1)[1].split(";"):
                gfail[f.split("=")[0]] += 1
    lines += ["| guard failed | merges |", "|---|---|"]
    for k, v in sorted(gfail.items()):
        lines.append(f"| {k} | {v} |")
report = "\n".join(lines) + "\n"
(HERE / "processing" / f"02_repair_report_{SP}.md").write_text(report)
print(report)
