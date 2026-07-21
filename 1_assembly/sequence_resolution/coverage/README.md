# Coverage of new asms

## Z chromosome ID + validation
>_Hydrophis_ is ZW (females heterogametic). In a sequenced female (ZW) the Z is present in one copy vs two for autosomes, so the differentiated Z shows ~0.5× autosomal depth. In a male (ZZ) the Z is diploid → ~1.0× (no coverage signal). Thus, a chromosome that drops to half depth only in the female, and is the same chromosome we named chZ, is positive proof.

**Method**
  * Per-chromosome median depth from mosdepth.global.dist (median = max depth with cumulative fraction ≥ 0.5).
  * Coverage contigs matched to renamed chromosomes by exact length.
  * chZ:autosome ratio = chZ median / median of the seven autosomal macro medians.
  * 100 kb windowed coverage (mosdepth per-base, length-weighted), normalised to the autosomal median, for the per-assembly profile figures.



**Result**
| assembly | species            | platform    | autosome median | chZ median | **chZ : autosome** | inferred sex    |
| -------- | ------------------ | ----------- | --------------- | ---------- | ------------------ | --------------- |
| hcy      | *H. cyanocinctus*  | PacBio HiFi | 20×             | 20×        | **1.00**           | male (ZZ)       |
| hcurw    | *H. curtus* (west) | Nanopore    | 14×             | 13×        | **0.93**           | male (ZZ)       |
| **hmaj** | ***H. major***     | PacBio HiFi | 31×             | 16×        | **0.52**           | **female (ZW)** |
| horn     | *H. ornatus*       | Nanopore    | 58×             | 57×        | **0.98**           | male (ZZ)       |


### Interpretation
- **hmaj (female): proof.** `chZ` is the **only** chromosome at half coverage (0.52×, uniform across its whole length),
- **hcy, hcurw, horn (males): consistent.** `chZ` sits at full diploid depth (~1.0×) like the autosomes — the expected pattern for ZZ males, in which the Z carries **no** coverage signal. Crucially, **no autosome** in any assembly drops to half, so nothing contradicts the chZ assignment. (A minor ~0.85× dip on `ch4` appears in several assemblies — a low-mappability/partially-collapsed autosomal region, not sex-linked.)
