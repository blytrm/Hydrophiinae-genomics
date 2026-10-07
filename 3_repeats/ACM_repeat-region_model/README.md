
## ACM Repeat Region Test (Model)

```
INCOMPLETE
```

- Each test compares region R with the same family elsewhere on the chromosomes, so the family’s general insertion habits cancel out.

| Test                | Null hypothesis                                                                                | Statistic                                               |
| ------------------- | ---------------------------------------------------------------------------------------------- | ------------------------------------------------------- |
| A. Accumulation     | Copy count NR∼NB(Jm, Jk)NR​∼NB(Jm,Jk), fitted to background 1 Mb windows (allows for clumping) | Enrichment NR/JmNR​/Jm                                  |
| C. Tandem structure | Share of copies with a neighbour ≤200 bp away is the same as elsewhere                         | Fisher exact test                                       |
| M. Maintenance      | Divergence of copies in R has the same distribution as elsewhere                               | Mann–Whitney test; rank-biserial rr; spread (IQR) ratio |

Maintenance is read through $\frac{dK}{dt} = \mu -hK$ where $K$ is divergence, $\mu$  the substitution rate and $h$ the homogenisation rate:

- With $h =0, K$ just tracks how old each copy is.
- With $h > 0$ , copies are held near $K*=\frac{\mu}{h}$ so divergence is lower and tighter.
- A region is “accumulating + maintained” when it passes both A and M. Test C is reported alongside as the likely mechanism.
- _‘chromosome ends’ = last 5Mb at one end of chromosome_


#### example:

| Family                      | Region           | Enrichment | Close pairs here vs elsewhere | Median divergence here vs elsewhere | IQR ratio | Verdict                               |
| --------------------------- | ---------------- | ---------- | ----------------------------- | ----------------------------------- | --------- | ------------------------------------- |
| `TE_00006440_LTR`           | chromosome ends  | 35×        | 61% vs 25%                    | 3.5% vs 4.5%                        | 0.06      | accumulating + maintained, tandem     |
| `TE_00007750_LTR`           | chromosome ends  | 21×        | 21% vs 30%                    | 14.5% vs 17.6%                      | 0.37      | accumulating + maintained, not tandem |
| `TE_00007750_LTR`           | unplaced contigs | 74×        | 16% vs 27%                    | not lower                           | ~1.1      | accumulating, copies not younger      |
| `TE_00006109_LTR` (control) | either           | ≤1×        | no difference                 | no difference                       | ~0.9      | no signal                             |

>“This script asks: are there unusually many copies of one TE family in a region, and are those copies unusually similar?” 1. reads repeatmasker hits for 3 families 2. merges split pieces into whole elements 3. for each family x each region → **ACM** 4. saves a table + a figure 5. prints summary

`columns`

- enrichment, q_accum → test A
- frac_close, clustering, q_gap, tandem → test C
- K_R, K_B, dK, q_maint → test M
- verdict → the one-line answer
- depth_ratio → <1 means not assembly-collapsed (good for trusting counts)
- q-value: like a p-value corrected for many tests. Rough habit: q < 0.05 = “we call this real.”


`figure`

- **left scatter** → median divergence region - background (%) // enrichment … **right = $\uparrow$ enrichment ; down = younger/more similar ; purple : ==accumulating + maintained [bottom right]**==
- **middle gaps** → purple line rises faster than grey background → copies are closer
- **right histogram** →

