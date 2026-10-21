#!/usr/bin/env python3

"""Parse HyPhy selection-suite JSONs -> compact TSV + summary.txt.
Keeps multi-MB JSONs OUT of the controller context: emit only headline stats,
significant-site tables, and selected-branch lists. Run with -I (untrusted dir).
"""
import sys, json, os, glob

HY, OUT = sys.argv[1], sys.argv[2]
    # python3 -I parse_hyphy <hyphy_dir> <out_dir>
os.makedirs(OUT, exist_ok=True)
summary = []
def S(*a): summary.append(" ".join(str(x) for x in a))


def load(path):
    try:
        with open(path) as f: return json.load(f)
    except Exception as e:
        S("  !! load fail", os.path.basename(path), repr(e)); return None


def find(method_dir, *pats):
    for p in pats:
        g = sorted(glob.glob(os.path.join(HY, method_dir, p)))
        if g: return g[0]
    return None



def mle_rows(j):
    """Return (headers, rows) from a MLE block."""
    mle = j.get("MLE", {})
    hdr = [h[0] for h in mle.get("headers", [])]
    rows = mle.get("content", {}).get("0", [])
    return hdr, rows



def col(hdr, *names):
    for n in names:
        for i, h in enumerate(hdr):
            if h.strip().lower() == n.lower(): return i
    for n in names:                      # substring fallback
        for i, h in enumerate(hdr):
            if n.lower() in h.strip().lower(): return i
    return None


# FEL
f = find("fel", "*.FEL.json")
if f:
    j = load(f)
    if j:
        hdr, rows = mle_rows(j)
        ia, ib, ip = col(hdr,"alpha"), col(hdr,"beta"), col(hdr,"p-value","p value")
        pos=neg=0; out=[]
        for k,r in enumerate(rows,1):
            if ip is None: break
            p=r[ip]; a=r[ia]; b=r[ib]
            if p<=0.1:
                if b>a: pos+=1; out.append((k,a,b,p,"pos"))
                elif b<a: neg+=1; out.append((k,a,b,p,"neg"))
        with open(f"{OUT}/fel_sites.tsv","w") as o:
            o.write("site\talpha\tbeta\tp\tclass\n")
            for s in out: o.write("\t".join(map(str,s))+"\n")
        S(f"FEL (p<=0.1): {pos} positive, {neg} purifying sites  [fel_sites.tsv]")

# MEME
f = find("meme", "*.MEME.json")
if f:
    j = load(f)
    if j:
        hdr, rows = mle_rows(j)
        ip = col(hdr,"p-value","p value")
        ibp= col(hdr,"beta+","&beta;+")  # positive-sel rate
        sig=[]
        for k,r in enumerate(rows,1):
            if ip is None: break
            if r[ip]<=0.05: sig.append((k, r[ip]))
        with open(f"{OUT}/meme_sites.tsv","w") as o:
            o.write("site\tp\n")
            for s in sig: o.write(f"{s[0]}\t{s[1]}\n")
        S(f"MEME (p<=0.05): {len(sig)} episodic-diversifying sites  [meme_sites.tsv]")

# contrast-FEL
f = find("cfel", "*.CFEL.json","*contrast*.json","*.FEL.json")
if f and "cfel" in f:
    j = load(f)
    if j:
        hdr, rows = mle_rows(j)
        ip = col(hdr,"p-value","q-value","p value")
        sig=[(k,r[ip]) for k,r in enumerate(rows,1) if ip is not None and r[ip]<=0.05]
        with open(f"{OUT}/cfel_sites.tsv","w") as o:
            o.write("site\tp\n")
            for s in sig: o.write(f"{s[0]}\t{s[1]}\n")
        S(f"contrast-FEL (p<=0.05): {len(sig)} sites differ Test vs Reference  [cfel_sites.tsv]")

# BUSTED / BUSTED-PH
def busted_report(dirn, label):
    f = find(dirn, "*.BUSTED.json")
    if not f: return
    j = load(f)
    if not j: return
    tr = j.get("test results", {})
    p = tr.get("p-value")
    S(f"{label}: overall p-value = {p}")
    # BUSTED-PH extra tests
    for key in j:
        if key.lower().startswith("test results") and key!="test results":
            S(f"   {key}: p={j[key].get('p-value')}")
busted_report("busted","BUSTED (marine FG episodic selection)")
busted_report("bustedph","BUSTED-PH (trait assoc)")

# RELAX 
f = find("relax","*.RELAX.json")
if f:
    j = load(f)
    if j:
        tr = j.get("test results", {})
        K = tr.get("relaxation or intensification parameter")
        p = tr.get("p-value")
        S(f"RELAX: K={K}  p={p}  ({'intensification K>1' if (K or 0)>1 else 'relaxation K<1'} on Test)")

# aBSREL
f = find("absrel","*.ABSREL.json")
if f:
    j = load(f)
    if j:
        ba = j.get("branch attributes", {}).get("0", {})
        sel=[]
        for br,at in ba.items():
            p = at.get("Corrected P-value", at.get("corrected P-value"))
            if p is not None and p<=0.05: sel.append((br,p))
        sel.sort(key=lambda x:x[1])
        with open(f"{OUT}/absrel_branches.tsv","w") as o:
            o.write("branch\tcorrected_p\n")
            for b,p in sel: o.write(f"{b}\t{p}\n")
        S(f"aBSREL: {len(sel)}/{len(ba)} branches under episodic selection (corr p<=0.05)  [absrel_branches.tsv]")

# GARD 
f = os.path.join(HY,"gard.json")
if os.path.exists(f):
    j = load(f)
    if j:
        bps = j.get("breakpointData") or j.get("breakpoints")
        imp = j.get("improvements")
        S(f"GARD: {len(bps)} segments -> {len(bps)-1} breakpoint(s)")

with open(f"{OUT}/summary.txt","w") as o:
    o.write("HyPhy selection suite — reduced tree (181 tips, 864 codons)\n")
    o.write("Test=marine Hydrophis, Reference=terrestrial snakes\n"+"="*60+"\n")
    o.write("\n".join(summary)+"\n")
print("\n".join(summary))
print(f"\nwrote TSVs + summary.txt to {OUT}")
