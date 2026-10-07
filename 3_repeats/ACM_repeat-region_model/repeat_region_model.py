"""Test whether a genomic region is accumulating and maintaining copies of one repeat family.

For region R and family f, copies (RepeatMasker hits merged by element ID) are tested against
three nulls:

  A  accumulation  copies in R ~ negative binomial fitted to same-size background windows
  C  clustering    inter-copy gaps ~ Exp(lambda) under uniform placement; adjacent close copies
                   share strand with p = 0.5
  M  maintenance   divergence from consensus in R ~ same distribution as f in the background
                   (one-sided Mann-Whitney, R lower)

Maintenance is read through dK/dt = mu - h*K: without homogenisation (h = 0) a copy's divergence
grows with its age, so R and the background share f's age distribution; homogenisation in R
(h > 0) caps copies near K* = mu/h, giving lower and tighter divergence.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

OUT_COLS = ['score', 'div', 'del', 'ins', 'chrom', 'start', 'end', 'left', 'strand',
            'name', 'cls', 'rbeg', 'rend', 'rleft', 'id', 'overlap']


def read_rm_elements(path, families=None, classes=None, max_gap=200, slack=20) -> pd.DataFrame:
    """RepeatMasker .out -> one row per element.

    Neighbouring hits are one element if they share an RM ID, or if they sit within `max_gap` bp
    on the same strand and continue along the consensus (RM often splits one copy under new IDs).
    """
    hits = pd.read_csv(path, sep=r'\s+', skiprows=3, header=None, names=OUT_COLS,
                       usecols=['div', 'chrom', 'start', 'end', 'strand', 'name', 'cls', 'rbeg',
                                'rend', 'rleft', 'id'])
    if families is not None:
        hits = hits[hits.name.isin(families)]
    if classes is not None:
        hits = hits[hits.cls.str.split('/').str[0].isin(classes)]
    num = lambda s: s.astype(str).str.strip('()').astype(int)
    plus = hits.strand == '+'
    hits = hits.assign(cbeg=np.where(plus, num(hits.rbeg), num(hits.rleft)), cend=num(hits.rend),
                       bp=hits.end - hits.start + 1)
    hits = hits.sort_values(['name', 'chrom', 'start'], ignore_index=True)

    prev = hits.shift()
    adjacent = ((hits.name == prev.name) & (hits.chrom == prev.chrom)
                & (hits.strand == prev.strand) & (hits.start - prev.end - 1 <= max_gap))
    continues = np.where(plus, hits.cbeg >= prev.cend - slack, prev.cbeg >= hits.cend - slack)
    same_el = (hits.chrom == prev.chrom) & ((hits.id == prev.id) | (adjacent & continues))
    hits['element'] = (~same_el).cumsum()

    hits['div_bp'] = hits['div'] * hits.bp
    el = hits.groupby('element', sort=False).agg(
        chrom=('chrom', 'first'), name=('name', 'first'), cls=('cls', 'first'), start=('start', 'min'),
        end=('end', 'max'), strand=('strand', 'first'), bp=('bp', 'sum'),
        div_bp=('div_bp', 'sum'), n_frag=('bp', 'size'))
    el['div'] = el.pop('div_bp') / el.bp
    return el.reset_index(drop=True)


def intervals(chroms, starts, ends) -> pd.DataFrame:
    return pd.DataFrame({'chrom': chroms, 'start': starts, 'end': ends})


def in_region(el, region) -> np.ndarray:
    mid = ((el.start + el.end) // 2).to_numpy()
    chrom = el.chrom.to_numpy()
    mask = np.zeros(len(el), bool)
    for r in region.itertuples():
        mask |= (chrom == r.chrom) & (mid >= r.start) & (mid <= r.end)
    return mask


def complement(fai, chroms, region) -> pd.DataFrame:
    """Chromosome segments not covered by region."""
    rows = []
    for c in chroms:
        cuts = region[region.chrom == c].sort_values('start')
        pos = 1
        for r in cuts.itertuples():
            if r.start > pos:
                rows.append((c, pos, r.start - 1))
            pos = max(pos, r.end + 1)
        if pos <= fai[c]:
            rows.append((c, pos, int(fai[c])))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end'])


def accumulation_test(fam, region, background, window=1_000_000):
    """Copies in region vs. a negative binomial fitted to background windows of `window` bp."""
    L_R = int((region.end - region.start + 1).sum())
    n_R = int(in_region(fam, region).sum())
    mids = ((fam.start + fam.end) // 2).to_numpy()
    chrom = fam.chrom.to_numpy()
    counts = np.concatenate([
        np.histogram(mids[chrom == r.chrom], bins=np.arange(r.start, r.end + 1, window))[0]
        for r in background.itertuples() if r.end - r.start + 1 >= 2 * window])
    m, v = counts.mean(), counts.var(ddof=1)
    J = L_R / window
    expected = J * m
    if v > m:
        size = J * m ** 2 / (v - m)
        p = stats.nbinom.sf(n_R - 1, size, size / (size + expected))
    else:
        p = stats.poisson.sf(n_R - 1, expected)
    return dict(L_R=L_R, n_R=n_R, expected=expected,
                enrichment=n_R / expected if expected else np.inf, p_accum=p)


def gaps_to_next(sub):
    """Gap to the next copy on the same sequence, and whether the two share strand."""
    sub = sub.sort_values(['chrom', 'start'])
    same = sub.chrom.to_numpy()[1:] == sub.chrom.to_numpy()[:-1]
    gaps = np.clip(sub.start.to_numpy()[1:] - sub.end.to_numpy()[:-1] - 1, 0, None)[same]
    same_strand = (sub.strand.to_numpy()[1:] == sub.strand.to_numpy()[:-1])[same]
    return gaps, same_strand


def clustering_test(fam, region, background, max_gap=200):
    """Share of copies with a neighbour <= max_gap bp away, region vs. the family elsewhere.

    The exponential expectation under uniform placement is reported for reference only: LTR
    families clump at Mb scale and nest, so even dispersed families beat it.
    """
    sub = fam[in_region(fam, region)]
    if len(sub) < 3:
        return dict(n_gaps=0, frac_close=np.nan, frac_close_bg=np.nan, clustering=np.nan,
                    vs_uniform=np.nan, p_gap=np.nan, frac_same_strand=np.nan)
    gaps, same_strand = gaps_to_next(sub)
    gaps_bg, _ = gaps_to_next(fam[in_region(fam, background)])
    close, close_bg = gaps <= max_gap, gaps_bg <= max_gap
    lam = len(sub) / ((region.end - region.start + 1).sum() - (sub.end - sub.start + 1).sum())
    table = [[close.sum(), (~close).sum()], [close_bg.sum(), (~close_bg).sum()]]
    return dict(
        n_gaps=len(gaps), frac_close=close.mean(), frac_close_bg=close_bg.mean(),
        clustering=close.mean() / close_bg.mean() if close_bg.mean() else np.inf,
        vs_uniform=close.mean() / (1 - np.exp(-lam * max_gap)),
        p_gap=stats.fisher_exact(table, alternative='greater')[1],
        frac_same_strand=same_strand[close].mean() if close.any() else np.nan)


def maintenance_test(fam, region, background):
    """Divergence of copies in region vs. the same family in the background."""
    kR = fam.loc[in_region(fam, region), 'div'].to_numpy()
    kB = fam.loc[in_region(fam, background), 'div'].to_numpy()
    if len(kR) < 3 or len(kB) < 3:
        return dict(K_R=np.nan, K_B=np.nan, dK=np.nan, rank_biserial=np.nan, iqr_ratio=np.nan,
                    p_maint=np.nan)
    u = stats.mannwhitneyu(kR, kB, alternative='less')
    return dict(K_R=np.median(kR), K_B=np.median(kB), dK=np.median(kR) - np.median(kB),
                rank_biserial=1 - 2 * u.statistic / (len(kR) * len(kB)),
                iqr_ratio=stats.iqr(kR) / stats.iqr(kB), p_maint=u.pvalue)


def depth_ratio(bam, region, idxstats, chroms):
    """Reads per bp in region / reads per bp on chromosomes; > 1 suggests collapsed copies."""
    regs = [f'{r.chrom}:{r.start}-{r.end}' for r in region.itertuples()]
    n = int(subprocess.run(['samtools', 'view', '-c', '-F', '0x904', str(bam), *regs],
                           capture_output=True, text=True, check=True).stdout)
    bg = idxstats[idxstats.chrom.isin(chroms)]
    return (n / (region.end - region.start + 1).sum()) / (bg.mapped.sum() / bg.len.sum())


def verdict(r, alpha=0.05):
    acc = r.q_accum < alpha and r.enrichment > 1
    maint = r.q_maint < alpha
    if acc and maint:
        return 'accumulating + maintained'
    if acc:
        return 'accumulating, copies not younger'
    if maint:
        return 'maintained, not enriched'
    return 'no signal'


def test_region(fam, region, background, **kw):
    return {**accumulation_test(fam, region, background, kw.get('window', 1_000_000)),
            **clustering_test(fam, region, background, kw.get('max_gap', 200)),
            **maintenance_test(fam, region, background)}


def add_q_values(res, alpha=0.05):
    for p in ('p_accum', 'p_gap', 'p_maint'):
        ok = res[p].notna()
        res.loc[ok, 'q' + p[1:]] = multipletests(res.loc[ok, p], method='fdr_bh')[1]
    res['verdict'] = res.apply(verdict, axis=1, alpha=alpha)
    res['tandem'] = (res.q_gap < alpha) & (res.clustering > 1)
    return res


def plot_model(res, el, regions, backgrounds, case, outpng):
    """Accumulation-maintenance plane, gap ECDF vs. exponential null, divergence distributions."""
    fam_name, reg_name = case
    colors = {'accumulating + maintained': '#A858A8', 'accumulating, copies not younger': '#FC907C',
              'maintained, not enriched': '#F7C48D', 'no signal': '#BBBBBB'}
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))

    a = ax[0]
    for v, g in res.groupby('verdict'):
        a.scatter(np.log2(g.enrichment.clip(lower=1 / 64)), g.dK, s=28, color=colors[v], label=v,
                  edgecolor='k', linewidth=0.3)
    for r in res[res.region.isin(['unplaced', 'all_qter', 'all_pter'])].itertuples():
        a.annotate(f'{r.family.replace("TE_0000", "")}\n{r.region}', (np.log2(max(r.enrichment, 1 / 64)), r.dK),
                   fontsize=6, xytext=(3, 3), textcoords='offset points')
    a.axvline(0, color='k', lw=0.5); a.axhline(0, color='k', lw=0.5)
    a.set_xlabel('log2 enrichment (accumulation)'); a.set_ylabel('median div. region - background (%)')
    a.set_title('Accumulation vs maintenance', fontsize=11); a.legend(frameon=False, fontsize=7)

    fam = el[el.name == fam_name]
    region, background = regions[reg_name], backgrounds[reg_name]
    sub = fam[in_region(fam, region)]
    lam = len(sub) / ((region.end - region.start + 1).sum() - (sub.end - sub.start + 1).sum())
    x = np.logspace(0, 7, 200)
    a = ax[1]
    for g, color, label in ((gaps_to_next(sub)[0], '#A858A8', reg_name),
                            (gaps_to_next(fam[in_region(fam, background)])[0], '#888888', 'background')):
        a.plot(np.sort(g), np.arange(1, len(g) + 1) / len(g), color=color, label=label)
    a.plot(x, 1 - np.exp(-lam * x), color='k', ls='--', label=f'uniform placement in {reg_name}')
    a.axvline(200, color='k', lw=0.5, ls=':')
    a.set_xscale('log'); a.set_xlabel('gap to next copy (bp)'); a.set_ylabel('cumulative fraction')
    a.set_title(f'Clustering: {fam_name}, {reg_name}', fontsize=11); a.legend(frameon=False, fontsize=8)

    a = ax[2]
    bins = np.arange(0, 41, 1)
    a.hist(fam.loc[in_region(fam, background), 'div'], bins=bins, density=True, alpha=0.5,
           color='#BBBBBB', label='background')
    a.hist(sub['div'], bins=bins, density=True, alpha=0.6, color='#A858A8', label=reg_name)
    a.set_xlabel('divergence from consensus (%)'); a.set_ylabel('density')
    a.set_title(f'Maintenance: {fam_name}', fontsize=11); a.legend(frameon=False, fontsize=8)
    plt.tight_layout(); fig.savefig(outpng, dpi=200); plt.close(fig)


SANDERS = Path('/hpcfs/users/a1864358/sanders_lab')
RESULTS = SANDERS / 'repeatdensity-atv2rs/workspace/res'
SPECIES = ['horn', 'hmaj', 'hcy', 'hcurw', 'hcure']
SP_LAB = {'horn': 'H. ornatus', 'hmaj': 'H. major', 'hcy': 'H. cyanocinctus',
          'hcurw': 'H. curtus (W)', 'hcure': 'H. curtus (E)'}
# class -> (RepeatMasker class prefix, family-name test). LTR: flank families (*_LTR) only, so one
# provirus is not counted again as its *_INT body.
CLASSES = {'LTR': ('LTR', lambda n: n.str.endswith('_LTR')), 'LINE': ('LINE', lambda n: n.notna())}
MIN_COPIES = 300  # merged elements per family in the species
SEXCHECK = SANDERS / 'asm/files/repeat-anno/ra2/results/06_sexcheck'  # horn reads only
COLORS = {'accumulating + maintained': '#A858A8', 'accumulating, copies not younger': '#FC907C',
          'maintained, not enriched': '#F7C48D', 'no signal': '#BBBBBB'}


def build_regions(fai, W=5_000_000):
    """Termini of chromosomes >= 20 Mb, all termini pooled, and unplaced sequence."""
    chroms = list(fai.index[fai >= 20e6])
    unplaced = list(fai.index[fai < 20e6])
    regions = {'unplaced': intervals(unplaced, 1, fai[unplaced].to_numpy())}
    for c in chroms:
        regions[f'{c}_pter'] = intervals([c], 1, W)
        regions[f'{c}_qter'] = intervals([c], int(fai[c]) - W + 1, int(fai[c]))
    regions['all_pter'] = pd.concat([regions[f'{c}_pter'] for c in chroms], ignore_index=True)
    regions['all_qter'] = pd.concat([regions[f'{c}_qter'] for c in chroms], ignore_index=True)
    backgrounds = {k: complement(fai, chroms, v) for k, v in regions.items()}
    return chroms, regions, backgrounds


def run_species(sp, min_copies=MIN_COPIES):
    """Region model for one assembly, every family >= min_copies, LTR and LINE tested separately."""
    fai = pd.read_csv(SANDERS / f'resources/final-asms/10-{sp}-final.renamed.fa.fai', sep='\t',
                      header=None, usecols=[0, 1], names=['chrom', 'len']).set_index('chrom').len
    el_all = read_rm_elements(SANDERS / f'asm/files/repeat-anno/edta_10-{sp}/repeatmasker_out/'
                              f'10-{sp}-final.renamed.fa.out', classes=['LTR', 'LINE'])
    chroms, regions, backgrounds = build_regions(fai)
    if sp == 'horn':
        idx = pd.read_csv(SEXCHECK / 'idxstats.txt', sep='\t', header=None,
                          names=['chrom', 'len', 'mapped', 'unmapped'])
        depth = {k: depth_ratio(SEXCHECK / 'sub.bam', v, idx, chroms) for k, v in regions.items()}
    else:
        depth = dict.fromkeys(regions, np.nan)

    out = []
    for klass, (prefix, keep) in CLASSES.items():
        el = el_all[(el_all.cls.str.split('/').str[0] == prefix) & keep(el_all.name)]
        counts = el.name.value_counts()
        families = list(counts.index[counts >= min_copies])
        print(f'{sp} {klass}: {len(families)} families >= {min_copies} copies', flush=True)
        by_fam = {f: g for f, g in el.groupby('name') if f in families}
        rows = [{'species': sp, 'class': klass, 'family': f, 'region': k, 'depth_ratio': depth[k],
                 **test_region(by_fam[f], region, backgrounds[k])}
                for f in families for k, region in regions.items()]
        out.append(add_q_values(pd.DataFrame(rows)))  # BH over this species x class
    return pd.concat(out, ignore_index=True)


def plot_facets(res, outpng):
    """Accumulation vs maintenance; rows = species, columns = LTR | LINE."""
    fig, axes = plt.subplots(len(SPECIES), 2, figsize=(11, 3.4 * len(SPECIES)), sharex=True,
                             sharey=True, squeeze=False)
    for i, sp in enumerate(SPECIES):
        for j, klass in enumerate(CLASSES):
            a, d = axes[i, j], res[(res.species == sp) & (res['class'] == klass)]
            for v in COLORS:
                g = d[d.verdict == v]
                a.scatter(np.log2(g.enrichment.clip(lower=1 / 64)), g.dK, s=10 if v == 'no signal' else 22,
                          color=COLORS[v], edgecolor='k' if v != 'no signal' else 'none',
                          linewidth=0.3, alpha=0.5 if v == 'no signal' else 0.9, label=v)
            a.axvline(0, color='k', lw=0.5); a.axhline(0, color='k', lw=0.5)
            n_fam = d.family.nunique()
            a.set_title(f'{SP_LAB[sp]} · {klass} ({n_fam} families)', fontsize=10)
            if i == len(SPECIES) - 1:
                a.set_xlabel('log2 enrichment (accumulation)')
            if j == 0:
                a.set_ylabel('median div. region - background (%)')
    axes[0, 1].legend(frameon=False, fontsize=7)
    plt.tight_layout(); fig.savefig(outpng, dpi=200); plt.close(fig)


def main():
    outdir = RESULTS / 'data'
    outdir.mkdir(parents=True, exist_ok=True)
    res = pd.concat([run_species(sp) for sp in SPECIES], ignore_index=True)
    res.to_csv(outdir / 'hydrophis_region_model_results.tsv', sep='\t', index=False)
    plot_facets(res, RESULTS / 'hydrophis_region_model_facets.png')
    print(res.groupby(['species', 'class']).verdict.value_counts().unstack(fill_value=0).to_string())
    return res


if __name__ == '__main__':
    main()
