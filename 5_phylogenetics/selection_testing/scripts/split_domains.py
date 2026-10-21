#!/usr/bin/env python3
"""
Per-sequence VFT vs 9Cys+7TM split from hmmsearch domtblout (best domain hit per seq per profile).
input from build_domain_hmms.sh    
"""

import sys

vft_f, ncd3g_f, tm7_f, faa_f, out_prefix = sys.argv[1:6] # from build_domain_hmms.sh


def best_hits(domtbl):
    # domtblout columns: target name(0) ... env_from(19) env_to(20) ... score(13)
    best = {}
    for line in open(domtbl):
        if line.startswith("#"):
            continue
        f = line.split()
        name = f[0]
        env_from, env_to = int(f[19]), int(f[20])
        score = float(f[13])
        if name not in best or score > best[name][2]:
            best[name] = (env_from, env_to, score)
    return best

vft = best_hits(vft_f)
ncd3g = best_hits(ncd3g_f)
tm7 = best_hits(tm7_f)

seqs = {}
name = None
for line in open(faa_f):
    line = line.rstrip()
    if line.startswith(">"):
        name = line[1:].split()[0]
        seqs[name] = []
    else:
        seqs[name].append(line)
seqs = {k: "".join(v) for k, v in seqs.items()}

rows = []
n_no_vft = 0
for name, seq in seqs.items():
    L = len(seq)
    v = vft.get(name)
    n = ncd3g.get(name)
    t = tm7.get(name)
    if v is None:
        n_no_vft += 1
        rows.append((name, L, "NA", "NA", "NA", "NA", "dropped_no_VFT"))
        continue
    vft_end = v[1]
    # core start = NCD3G start if present, else 7TM start, else just after VFT
    if n is not None:
        core_start = n[0]
    elif t is not None:
        core_start = t[0]
    else:
        core_start = vft_end + 1
    split = min(vft_end, core_start - 1) if core_start > vft_end else vft_end
    flag = "ok" if core_start > vft_end else "domains_overlap_used_vft_end"
    rows.append((name, L, v[0], vft_end, core_start, split, flag))

with open(f"{out_prefix}.tsv", "w") as fh:
    fh.write("seq\tlen\tvft_start\tvft_end\tcore_start\tsplit\tflag\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")

vft_fa = open(f"{out_prefix}.vft.faa", "w")
core_fa = open(f"{out_prefix}.core.faa", "w")
for name, L, vs, ve, cs, split, flag in rows:
    if flag == "dropped_no_VFT":
        continue
    seq = seqs[name]
    vft_fa.write(f">{name}\n{seq[:split]}\n")
    core_fa.write(f">{name}\n{seq[split:]}\n")
vft_fa.close()
core_fa.close()

print(f"total={len(rows)} dropped_no_VFT={n_no_vft} kept={len(rows)-n_no_vft}")
