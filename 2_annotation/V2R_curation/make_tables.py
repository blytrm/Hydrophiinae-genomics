#!/usr/bin/env python3
"""Turn the Class C pHMM scan into the curation tables.

Writes, into the results directory:

    classC_assign.tsv     one row per scanned protein: family call + margin
    v2r_catalogue.tsv     one row per locus, tiered
    v2r_master.tsv        every locus x overlapping transcript  <- filter this
    v2r_genes.tsv         the is_v2r rows of the master table
    classC_all_genes.tsv  every confident Class C call, rejects included
    chrom_summary.tsv     loci and genes per chromosome / scaffold
    summary.txt           parameters and every headline number
"""

import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from v2rlib import (COV_INTACT, COV_PARTIAL, EVALUE_CEIL, MIN_DELTA,
                    MIN_SINGLE_BITS, TARGET, coverage_class, read_domtbl,
                    read_tsv, union_len, write_tsv)

EGGNOG_COLS = ["eggNOG_OGs", "COG_category", "Description", "Preferred_name",
               "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "BRITE",
               "PFAMs"]

# A product label counts as agreeing with a family call when it names that
# family in the annotation's own vocabulary.
PRODUCT_HINTS = {
    "V2R_OlfC": ("vomeronasal", "vmn2r", "v2r"),
    "TAS1R": ("taste receptor type 1", "tas1r"),
    "GRM": ("metabotropic glutamate", "grm"),
    "CASR": ("calcium-sensing", "casr"),
    "GABBR": ("gamma-aminobutyric", "gabbr", "gaba"),
    "GPRC6A": ("group 6 member a", "gprc6"),
    "GPRC5": ("group 5 member", "gprc5", "retinoic acid-induced"),
}

CATALOGUE_HEADER = ["locus_id", "chrom", "start", "end", "strand", "length",
                    "tier", "family", "delta_bits", "n_exons", "n_tracks",
                    "tracks", "eviann_gene_id", "control_overlap", "flags"]

LOCUS_COLS = ["locus_id", "chrom", "seq_type", "locus_start", "locus_end",
              "locus_len_bp", "locus_tier", "locus_family", "locus_delta_bits",
              "n_primary_tracks", "tracks", "control_overlap",
              "n_transcripts_in_locus", "n_genes_in_locus", "miniprot_alns",
              "frameshifts", "inframe_stops", "miniprot_max_identity",
              "locus_flags"]

TX_COLS = ["transcript_id", "gene_id", "tx_start", "tx_end", "strand",
           "tx_span_bp", "n_exons", "protein_len_aa", "v2r_model_cov",
           "coverage_class", "has_7tm", "n_7tm_domains", "tm7_len_aa",
           "tm7_evalue",
           "tblastn_7tm_hit", "tblastn_7tm_coords",
           "pfam_7tm_hit", "pfam_7tm_coords",
           "eviann_evidence", "start_codon", "stop_codon",
           "gene_biotype", "eviann_label", "eviann_label_status",
           "found_by_phmm_prot", "tx_credited_to_locus", "classC_family",
           "best_bits", "delta_bits", "evalue", "confident", "is_v2r",
           "v2r_evidence", "excluded_reason", "eviann_product",
           "product_agreement"] \
    + ["eggnog_" + c.lower() for c in EGGNOG_COLS] + ["eggnog_arch_consistent"]

MASTER_HEADER = LOCUS_COLS + TX_COLS
COL = {name: i for i, name in enumerate(MASTER_HEADER)}


# Family assignment
def parse_tblout(path):
    """query -> [(family, evalue, bits)].  hmmscan --tblout is 18 whitespace
    fields then a free-text description: 0 target (family), 2 query, 4 E-value,
    5 bitscore."""
    hits = defaultdict(list)
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split(maxsplit=18)
        if len(f) >= 18:
            hits[f[2]].append((f[0], float(f[4]), float(f[5])))
    return hits


def assign(query_hits):
    """(family, bits, runner_up, runner_bits, delta, evalue, confident, call)."""
    kept = [h for h in query_hits if h[1] <= EVALUE_CEIL]
    if not kept:
        # Nothing cleared the ceiling: report the least-bad hit, call nothing.
        fam, ev, bits = min(query_hits, key=lambda h: h[1])
        return (fam, bits, "NONE", 0.0, 0.0, ev, False, "below_threshold")

    best = {}
    for fam, ev, bits in kept:
        if fam not in best or bits > best[fam][1]:
            best[fam] = (ev, bits)
    ranked = sorted(best.items(), key=lambda kv: kv[1][1], reverse=True)
    fam, (ev, bits) = ranked[0]

    if len(ranked) > 1:
        runner, (_, runner_bits) = ranked[1]
        delta = bits - runner_bits
        confident = delta >= MIN_DELTA
        return (fam, bits, runner, runner_bits, delta, ev, confident,
                fam if confident else "ambiguous")

    # One family, so there is no runner-up and no margin to test. Treating the
    # raw bitscore as the margin would mark every such hit maximally confident,
    # so gate on an absolute score and label the call to keep it auditable.
    confident = bits >= MIN_SINGLE_BITS
    return (fam, bits, "NONE", 0.0, bits, ev, confident,
            fam if confident else "single_family_weak")


def write_assignments(path, calls):
    with open(path, "w") as fh:
        fh.write("query_id\tbest_family\tbest_bits\trunner_up_family\t"
                 "runner_up_bits\tdelta_bits\tevalue\tconfident\tassignment\n")
        for q in sorted(calls):
            fam, bits, rf, rb, delta, ev, conf, call = calls[q]
            fh.write(f"{q}\t{fam}\t{bits:.1f}\t{rf}\t{rb:.1f}\t{delta:.1f}\t"
                     f"{ev:.2e}\t{int(conf)}\t{call}\n")


# Inputs
def load_fasta_meta(path):
    """transcript -> (product name from the Name:"..." header, protein length)."""
    products, lengths, tid = {}, {}, None
    for line in open(path):
        if line.startswith(">"):
            head = line[1:].rstrip()
            tid = head.split()[0]
            m = re.search(r'Name:"([^"]*)"', head)
            products[tid], lengths[tid] = (m.group(1) if m else ""), 0
        elif tid:
            lengths[tid] += len(line.strip().rstrip("*"))
    return products, lengths


def load_eggnog(path):
    """transcript -> the EGGNOG_COLS values, or {} if emapper has not run."""
    if not path.exists():
        return {}
    out, cols = {}, None
    for line in open(path):
        if line.startswith("#query"):
            cols = line.lstrip("#").rstrip("\n").split("\t")
        elif cols and not line.startswith("#") and line.strip():
            f = line.rstrip("\n").split("\t")
            row = dict(zip(cols, f))
            out[f[0]] = [row.get(c, "-") or "-" for c in EGGNOG_COLS]
    return out


def load_inputs(out, exon_counts):
    """Everything the tables are built from, as one namespace-ish dict."""
    d = {}
    d["loci"] = [(f[0], int(f[1]), int(f[2]), f[3])
                 for f in read_tsv(out / "loci.bed", header=False)]
    d["tracks"] = {f[0]: (int(f[1]), f[2], f[3] if len(f) > 3 else "")
                   for f in read_tsv(out / "locus_tracks.tsv")}

    locus_tx, locus_genes, phmm_hit = defaultdict(list), defaultdict(set), set()
    for locus, tid, gid, phmm in read_tsv(out / "locus_transcripts.tsv"):
        locus_tx[locus].append(tid)
        if gid:
            locus_genes[locus].add(gid)
        if phmm == "1":
            phmm_hit.add(tid)
    d.update(locus_tx=locus_tx, locus_genes=locus_genes, phmm_hit=phmm_hit)

    d["eviann_v2r"] = {f[0] for f in
                       read_tsv(out / "eviann_v2r_labelled.txt", header=False)}
    d["miniprot"] = {f[0]: f[1:] for f in read_tsv(out / "locus_miniprot.tsv")}
    d["locus_of"] = {f[0]: f[1] for f in read_tsv(out / "mrna_locus.tsv")}
    d["exons"] = {f[0].replace("Parent=", ""): int(f[1])
                  for f in read_tsv(exon_counts, header=False) if len(f) > 1}
    d["coords"] = {f[3]: (f[0], int(f[1]), int(f[2]), f[5], f[6])
                   for f in read_tsv(out / "eviann_mrna.bed", header=False)}
    d["attrs"] = {f[0]: f[1:] for f in read_tsv(out / "mrna_attrs.tsv")}

    # Internal stops are not measured from the protein: EviAnn only emits
    # complete ATG..stop ORFs, so the answer is always zero. The real disruption
    # signal is miniprot's, carried per locus as inframe_stops / frameshifts.
    d["products"], d["prot_len"] = load_fasta_meta(out / "candidates.faa")

    # How much of the V2R profile each protein covers -- the fragmentation
    # measure. A half-length model hit is half a gene.
    spans, model_len = defaultdict(list), {}
    for t, mlen, hf, ht in read_tsv(out / "v2r_domains.tsv"):
        spans[t].append((int(hf), int(ht)))
        model_len[t] = int(mlen)
    d["v2r_cov"] = {t: union_len(v) / model_len[t] for t, v in spans.items()}

    d["tm_spans"], _, d["tm_evalue"] = read_domtbl(out / "pfam_7tm.domtbl")

    # tBLASTn 7TM: genomic-space, can hit a locus with no EviAnn model. absent = no hit.
    d["tblastn_7tm"] = {f[0]: f[2] for f in read_tsv(out / "tblastn_7tm_hits.tsv")}
    d["eggnog"] = load_eggnog(out.parent.parent / "eggnog" / "out" / out.name
                              / "v2r.emapper.annotations")
    return d


# Per-locus catalogue
def best_calls(calls, locus_of):
    """locus -> (family, delta, transcript) for its strongest confident call.

    A locus can span a whole tandem array, so it may carry calls from several
    transcripts. Family, margin, transcript and exon count are all read from
    the same one -- mixing a family from one transcript with an exon count from
    another describes a gene that does not exist.
    """
    best = {}
    for q in sorted(calls):
        *_, conf, _call = calls[q]
        delta = calls[q][4]
        locus = locus_of.get(q)
        if conf and locus and (locus not in best or delta > best[locus][1]):
            best[locus] = (calls[q][0], delta, q)
    return best


def catalogue(d, best):
    """Tier every locus.

    T2 = confidently assigned V2R. T1 = a primary engine supports the locus but
    no confident call. T0 = confidently a different Class C family, or nothing.

    There is deliberately no structural gate. Exon count, profile coverage and
    7TM presence are all columns of the master table, so any integrity filter is
    a filter on the output rather than a second definition of T2.
    """
    rows = []
    for chrom, start, end, locus in d["loci"]:
        n_primary, track_list, control = d["tracks"].get(locus, (0, "", ""))
        family, delta, call_tx = best.get(locus, ("", 0.0, ""))
        flags = []

        if family == TARGET and delta >= MIN_DELTA:
            tier = "T2"
        elif family and delta >= MIN_DELTA:
            tier = "T0"
            flags.append(f"assigned_{family}")
        elif n_primary >= 1:
            tier = "T1"
            if not family:
                flags.append("unclassified")
        else:
            tier = "T0"

        # An overlapping olfactory-receptor call is a red flag, not support.
        if control:
            if tier == "T2":
                tier = "T1"
            flags.append(f"control_overlap:{control}")

        features = sorted(set(d["locus_tx"].get(locus, []))
                          | d["locus_genes"].get(locus, set()))
        rows.append([locus, chrom, start, end, ".", end - start, tier,
                     family or "NA", f"{delta:.1f}", d["exons"].get(call_tx, 0),
                     n_primary, track_list, ",".join(features) or "NA", control,
                     ";".join(flags)])
    rows.sort(key=lambda r: (r[1], r[2]))
    return rows


# Master table
def transcript_block(tid, locus, d, calls):
    """The transcript and call columns for one transcript at one locus."""
    fam, bits, _rf, _rb, delta, ev, conf, call = calls.get(
        tid, ("NA", 0.0, "NA", 0.0, 0.0, 1.0, False, "not_scanned"))
    chrom, start, end, strand, gene = d["coords"].get(tid, ("NA", 0, 0, ".", ""))

    cov = d["v2r_cov"].get(tid)
    cov_s = "NA" if cov is None else f"{cov:.3f}"

    if tid in d["tm_spans"]:
        has_7tm, n_tm = "yes", len(d["tm_spans"][tid])
        tm_len, tm_e = union_len(d["tm_spans"][tid]), f"{d['tm_evalue'][tid]:.2e}"
        # pfam is protein-space: comma-joined domain spans (usually one)
        pfam_coords = ",".join(f"{a}-{b}" for a, b in sorted(d["tm_spans"][tid]))
    else:
        has_7tm, n_tm, tm_len, tm_e, pfam_coords = "no", 0, 0, "NA", "NA"

    tblastn_coords = d["tblastn_7tm"].get(tid)
    tblastn_hit = "yes" if tblastn_coords else "no"
    tblastn_coords = tblastn_coords or "NA"

    # Retention. A panel V2R call is kept when the transcript is credited to
    # this locus. A transcript EviAnn labelled V2R that no candidate locus
    # covers is rescued too, because dropping it would delete the only record of
    # that gene -- but a confident call for another family (CASR) contradicts
    # EviAnn rather than merely failing to confirm it, and stays out.
    credited = "yes" if d["locus_of"].get(tid, "UNPLACED") == locus else "no"
    panel_v2r = call == TARGET and credited == "yes"
    if panel_v2r and locus != "UNPLACED":
        is_v2r, evidence, excluded = "yes", "panel", ""
    elif panel_v2r and tid in d["eviann_v2r"]:
        is_v2r, evidence, excluded = "yes", "eviann_rescue", ""
    else:
        is_v2r, evidence, excluded = "no", "none", "not_a_confident_V2R_call"

    # A canonical V2R is Venus-flytrap + cysteine-rich linker + Class C 7TM.
    # eggNOG assigns the architecture independently of our panel, so
    # disagreement is the cleanest flag for a misassigned protein -- ionotropic
    # glutamate receptors especially, which the panel has no profile to compete
    # against.
    egg = d["eggnog"].get(tid, ["NA"] * len(EGGNOG_COLS))
    pfams = egg[EGGNOG_COLS.index("PFAMs")]
    arch = ("NA" if pfams in ("NA", "-") else
            "yes" if "7tm_3" in pfams and "ANF_receptor" in pfams else "no")

    evidence_cls, start_c, stop_c, _n_ex, biotype, label = \
        d["attrs"].get(tid, ("NA",) * 6)
    # EviAnn writes the literal string "function unknown" when it could not
    # transfer a label, which is a different thing from having no label at all.
    label_status = ("none" if label in ("", "NA") else
                    "unknown" if label.lower() == "function unknown"
                    else "annotated")

    product = d["products"].get(tid, "")
    if not product or product == "function unknown":
        agreement = "no_annotation"
    else:
        hints = PRODUCT_HINTS.get(call)
        agreement = ("NA" if not hints else
                     "agree" if any(h in product.lower() for h in hints)
                     else "disagree")

    return [tid, gene or "NA", start, end, strand, end - start,
            d["exons"].get(tid, 0), d["prot_len"].get(tid, 0), cov_s,
            coverage_class(cov), has_7tm, n_tm, tm_len, tm_e,
            tblastn_hit, tblastn_coords, has_7tm, pfam_coords,
            evidence_cls, start_c, stop_c, biotype, label, label_status,
            "yes" if tid in d["phmm_hit"] else "no", credited,
            call, f"{bits:.1f}", f"{delta:.1f}", f"{ev:.2e}", int(conf),
            is_v2r, evidence, excluded, product or "NA", agreement] + egg + [arch]


def seq_type(chrom):
    return "chromosome" if chrom.startswith("ch") else "scaffold"


def master_table(d, calls, tier_of, flags_of, best):
    rows = []
    for chrom, start, end, locus in d["loci"]:
        n_primary, track_list, control = d["tracks"].get(locus, (0, "", ""))
        mp = d["miniprot"].get(locus, ["0", "0", "0", "0"])
        block = [locus, chrom, seq_type(chrom), start, end, end - start,
                 tier_of.get(locus, "NA"), best.get(locus, ("NA",))[0],
                 f"{best.get(locus, ('', 0.0))[1]:.1f}",
                 n_primary, track_list, control,
                 len(d["locus_tx"].get(locus, [])),
                 len(d["locus_genes"].get(locus, set())),
                 *mp, flags_of.get(locus, "")]
        tx = sorted(d["locus_tx"].get(locus, []))
        if tx:
            rows.extend(block + transcript_block(t, locus, d, calls) for t in tx)
        else:
            rows.append(block + ["NA"] * len(TX_COLS))

    # Transcripts EviAnn labelled V2R that fall outside every candidate locus.
    # They were scanned anyway, so they get a row marked UNPLACED rather than
    # being dropped for the accident of sitting between two loci.
    placed = {t for ts in d["locus_tx"].values() for t in ts}
    for tid in sorted(d["eviann_v2r"] - placed):
        chrom, start, end, _s, _g = d["coords"].get(tid, ("NA", 0, 0, ".", ""))
        rows.append(["UNPLACED", chrom, seq_type(chrom), start, end, end - start,
                     "NA", "NA", "0.0", 0, "", "", 1, 1, "0", "0", "0", "0",
                     "no_candidate_locus"]
                    + transcript_block(tid, "UNPLACED", d, calls))

    rows.sort(key=lambda r: (r[COL["chrom"]], r[COL["locus_start"]],
                             r[COL["transcript_id"]]))
    return rows


def chrom_summary(master):
    rows = []
    for c in sorted({r[COL["chrom"]] for r in master},
                    key=lambda c: (not c.startswith("ch"), c)):
        here = [r for r in master if r[COL["chrom"]] == c]
        real = {r[0]: r for r in here if r[0] != "UNPLACED"}
        tiers = [r[COL["locus_tier"]] for r in real.values()]
        v2r = [r for r in here if r[COL["is_v2r"]] == "yes"]
        rows.append([c, seq_type(c), len(real),
                     tiers.count("T2"), tiers.count("T1"), tiers.count("T0"),
                     len(v2r),
                     len({r[COL["gene_id"]] for r in v2r
                          if r[COL["gene_id"]] != "NA"})])
    return rows


# Summary
def summarise(species, d, calls, cat_rows, master, chrom_rows):
    out, v2r = [], [r for r in master if r[COL["is_v2r"]] == "yes"]

    def say(line=""):
        out.append(line)
        print(line)

    def tally(rows, name):
        c = defaultdict(int)
        for r in rows:
            c[r[COL[name]]] += 1
        return c

    def pct(n, total):
        return 100 * n / max(1, total)

    say(f"{species} V2R curation")
    say("=" * 62)
    say("PARAMETERS")
    say(f"  target family              {TARGET}")
    say(f"  confident call             margin >= {MIN_DELTA} bits over runner-up")
    say(f"  E-value ceiling per family {EVALUE_CEIL}")
    say(f"  single-family fallback     >= {MIN_SINGLE_BITS} bits (no runner-up)")
    say(f"  coverage classes           intact >= {COV_INTACT}, "
        f"partial >= {COV_PARTIAL}, else fragmented")
    say("  coverage is a label, not a gate: nothing is dropped for it")
    say()
    say("RESULTS")

    genes = {r[COL["gene_id"]] for r in v2r if r[COL["gene_id"]] != "NA"}
    say(f"candidate loci            {len(d['loci'])}")
    say(f"proteins scanned          {len(calls)}")
    say(f"V2R transcripts           {len(v2r)}  ({len(genes)} genes, "
        f"{len({r[0] for r in v2r})} loci)")

    rejected = defaultdict(int)
    for c in calls.values():
        if c[6] and c[7] not in (TARGET, "single_family_weak"):
            rejected[c[7]] += 1
    say(f"non-V2R Class C rejected  {sum(rejected.values())}  "
        + ", ".join(f"{n} {f}" for f, n in
                    sorted(rejected.items(), key=lambda kv: -kv[1])))
    for label in ("ambiguous", "single_family_weak"):
        say(f"{label:<25} {sum(1 for c in calls.values() if c[7] == label)}")

    tiers = defaultdict(int)
    for r in cat_rows:
        tiers[r[6]] += 1
    for t in ("T2", "T1", "T0"):
        say(f"  {t} loci                {tiers[t]:>5}")

    say()
    say("V2R transcript quality")
    cc = tally(v2r, "coverage_class")
    say(f"  intact / partial / fragmented {cc['intact']} / {cc['partial']} / "
        f"{cc['fragmented']}   (>= {COV_INTACT} / >= {COV_PARTIAL} / below)")
    ve = tally(v2r, "v2r_evidence")
    say(f"  retained by panel / rescued   {ve['panel']} / {ve['eviann_rescue']}")
    tm = tally(v2r, "has_7tm")
    say(f"  7TM domain present            {tm['yes']} yes / {tm['no']} no "
        f"({pct(tm['yes'], len(v2r)):.1f}%)")
    ex = tally(v2r, "n_exons")
    say("  exon counts                   "
        + ", ".join(f"{k}:{v}" for k, v in sorted(ex.items())[:12]))
    arch = tally(v2r, "eggnog_arch_consistent")
    if arch["yes"] or arch["no"]:
        say(f"  eggNOG architecture consistent {arch['yes']} yes / {arch['no']} no"
            "   (7tm_3 + ANF_receptor present)")

    say()
    say("Evidence class EviAnn used for each model")
    ec = tally(v2r, "eviann_evidence")
    for k in ("complete", "protein_only", "transcript_only"):
        say(f"  {k:<16} {ec[k]:>4}  ({pct(ec[k], len(v2r)):.1f}%)")
    say("  complete = transcript AND protein evidence; protein_only = homology")

    say()
    say("Agreement with EviAnn's own homology labels")
    pa = tally(v2r, "product_agreement")
    say(f"  {pa['agree']} agree / {pa['disagree']} disagree / "
        f"{pa['no_annotation']} unannotated   = "
        f"{pct(pa['agree'], pa['agree'] + pa['disagree']):.1f}% of annotated")
    ls = tally(v2r, "eviann_label_status")
    say(f"  EviAnn label: {ls['annotated']} carried a homology label, "
        f"{ls['unknown']} \"function unknown\", {ls['none']} none")

    say()
    say("Which engine found the V2R set")
    for t, n in sorted(tally(v2r, "tracks").items(), key=lambda kv: -kv[1]):
        say(f"  {n:>4}  ({pct(n, len(v2r)):>5.1f}%)  {t}")
    direct = sum(1 for r in v2r if r[COL["found_by_phmm_prot"]] == "yes")
    say(f"  {direct} of {len(v2r)} were hit by the proteome pHMM directly")

    say()
    say("Placement in the assembly")
    on_chrom = sum(1 for r in v2r if r[COL["seq_type"]] == "chromosome")
    say(f"  {on_chrom} V2R transcripts on chromosomes, "
        f"{len(v2r) - on_chrom} on scaffolds")
    for row in chrom_rows:
        if row[1] == "chromosome":
            say(f"  {row[0]:<12} loci {row[2]:>4}  T2 {row[3]:>3} T1 {row[4]:>3} "
                f"T0 {row[5]:>3}  V2R tx {row[6]:>3}  genes {row[7]:>3}")
    scaf = [r for r in chrom_rows if r[1] == "scaffold"]
    if scaf:
        say(f"  {len(scaf)} scaffolds  loci {sum(r[2] for r in scaf):>4}  "
            f"T2 {sum(r[3] for r in scaf):>3} T1 {sum(r[4] for r in scaf):>3} "
            f"T0 {sum(r[5] for r in scaf):>3}  V2R tx {sum(r[6] for r in scaf):>3}"
            f"  genes {sum(r[7] for r in scaf):>3}   (detail in chrom_summary.tsv)")
    return "\n".join(out) + "\n"


def main():
    out, exon_counts = Path(sys.argv[1]), Path(sys.argv[2])
    species = out.name

    calls = {q: assign(h) for q, h in parse_tblout(out / "classC.tblout").items()}
    write_assignments(out / "classC_assign.tsv", calls)

    d = load_inputs(out, exon_counts)
    best = best_calls(calls, d["locus_of"])

    cat_rows = catalogue(d, best)
    write_tsv(out / "v2r_catalogue.tsv", CATALOGUE_HEADER, cat_rows)
    tier_of = {r[0]: r[6] for r in cat_rows}
    flags_of = {r[0]: r[14] for r in cat_rows}

    master = master_table(d, calls, tier_of, flags_of, best)
    write_tsv(out / "v2r_master.tsv", MASTER_HEADER, master)
    write_tsv(out / "v2r_genes.tsv", MASTER_HEADER,
              [r for r in master if r[COL["is_v2r"]] == "yes"])
    write_tsv(out / "classC_all_genes.tsv", MASTER_HEADER,
              [r for r in master if r[COL["confident"]] == 1])

    chrom_rows = chrom_summary(master)
    write_tsv(out / "chrom_summary.tsv",
              ["sequence", "seq_type", "n_loci", "T2", "T1", "T0",
               "n_v2r_transcripts", "n_v2r_genes"], chrom_rows)

    (out / "summary.txt").write_text(
        summarise(species, d, calls, cat_rows, master, chrom_rows))


if __name__ == "__main__":
    main()
