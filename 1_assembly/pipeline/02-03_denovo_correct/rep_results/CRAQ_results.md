### Representative data from CRAQ output

![craq_genome_summary.png]

Headers example from some of the CRAQ output reports

**_Hydrophis ornatus_**
```bash
Short Report:
#Chr	Covered.Rate	Low-confident.Rate	Avg.CRH	Avg.CSH	Avg.CRE(R-AQI)	Avg.CSE(S-AQI)
Genome	0.999605365774614	1.39007425384429e-05	0	0.00748116689176376	0.631766731707326(93.8777602221448)	0.00374058344588188(99.6266403821509)
```
- R-AQI = 93.88
- S-AQI = 99.63

**_Hydrophis cyanocinctus_**
```bash
Short Report:
#Chr	Covered.Rate	Low-confident.Rate	Avg.CRH	Avg.CSH	Avg.CRE(R-AQI)	Avg.CSE(S-AQI)
Genome	0.999737172607808	9.37072020906479e-05	0	0.0174022221236773	1.88381903066833(82.8298317153782)	0.00380673608955441(99.6200500344962)
```
- R-AQI = 82.83
- S-AQI = 99.62

**_Hydrophis curtus (west)**
```bash
Short Report:
#Chr	Covered.Rate	Low-confident.Rate	Avg.CRH	Avg.CSH	Avg.CRE(R-AQI)	Avg.CSE(S-AQI)
Genome	0.998704698019341	1.88862448978437e-05	0	0.0208434855616719	1.96859302157701(82.1306181191359)	0.0162660926540106(98.3865485842432)
```
- R-AQI = 82.13
- S-AQI = 98.39

---

-> CRAQ assembly QC genome summary

```python
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm

genomes = ["hcy", "hcur-w", "horn"]
raqi    = [82.8, 82.1, 93.9]
saqi    = [99.6, 98.4, 99.6]
covered = [99.974, 99.870, 99.961]
lowconf = [9.37e-3, 1.89e-3, 1.39e-3]

c_metric = cm.lipari(np.linspace(0.25, 0.70, 2))
c_genome = cm.lipari(np.linspace(0.20, 0.80, 3))

plt.rcParams.update({
    "font.family": "serif", "font.size": 10,
    "axes.linewidth": 0.8, "axes.edgecolor": "0.3",
    "axes.grid": True, "grid.color": "0.9", "grid.linewidth": 0.6,
    "axes.axisbelow": True,
})

fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
x = np.arange(len(genomes))

w = 0.38
ax[0].bar(x - w/2, raqi, w, color=c_metric[0], label="R-AQI (structural)")
ax[0].bar(x + w/2, saqi, w, color=c_metric[1], label="S-AQI (small-scale)")
for i, (r, s) in enumerate(zip(raqi, saqi)):
    ax[0].text(i - w/2, r + 1, f"{r:.1f}", ha="center", va="bottom", fontsize=8)
    ax[0].text(i + w/2, s + 1, f"{s:.1f}", ha="center", va="bottom", fontsize=8)
ax[0].set(ylim=(0, 108), ylabel="AQI (0–100, higher = better)")
ax[0].legend(frameon=False, fontsize=8, loc="lower left")

ax[1].bar(x, covered, color=c_genome, width=0.6)
for i, v in enumerate(covered):
    ax[1].text(i, v + 0.004, f"{v:.3f}%", ha="center", va="bottom", fontsize=8)
ax[1].set(ylim=(99.5, 100.02), ylabel="genome covered (%)")

ax[2].bar(x, lowconf, color=c_genome, width=0.6)
for i, v in enumerate(lowconf):
    ax[2].text(i, v, f"{v:.2e}", ha="center", va="bottom", fontsize=8)
ax[2].set(ylabel="low-confidence (%)")

for a in ax:
    a.set_xticks(x); a.set_xticklabels(genomes)
    a.grid(axis="x", visible=False)
    a.spines[["top", "right"]].set_visible(False)
    a.margins(x=0.12)

fig.tight_layout()
fig.savefig("/scratchdata1/users/a1864358/sanders_lab/craq_plots/craq_genome_summary_lipari.png",
            dpi=200, bbox_inches="tight")
print("saved")
```
