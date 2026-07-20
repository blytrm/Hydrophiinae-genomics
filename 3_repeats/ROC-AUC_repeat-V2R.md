
# ROC-AUC : V2R bin classifier
---
- **A confound check** :
   ## _`Can Assembly/Sequence-Complexity Features Predict V2R Presence?`_
---

### `Plan`
- L1-logistic ; random forest ; gradient boosting <- is predictability = model-specific || property of the data
- ROC-AUC
  - 1Mb bins
  - **Features:**
      - mappability (GenMap)
      - coverage
      - coverage absolute deviation (vs genome median)
      - shannon normalised entropy (RepeatObserver)
      - tandem repeat fraction (TRF)
      - lai (LAI)
      - tBLASTn + pHMM loci / bin
  
- Each row is one **1 Mb genomic bin** on ch2 or chZ.
- **Label** $y \in \{0,1\}$: does the bin contain ≥1 V2R?
- **Features** $X$

---

**`ROC-AUC`** classifier outputs a **score** $s(x) \in [0,1]$ (predicted probability of V2R) -> a threshold $t$ -> call the bin positive if $s(x) \ge t$.

$$
\mathrm{AUC}
= \int_0^1 \mathrm{TPR}\bigl(\mathrm{FPR}^{-1}(u)\bigr)\, du
$$

Equivalently:

$$
\mathrm{AUC}
= P\bigl(s_+ > s_-\bigr)
$$

where $s_+$ is a random **positive** bin’s score and $s_-$ a random **negative**.

-> **_`So AUC is the probability that the model ranks a real V2R bin above an empty one.`_**

---

| AUC | Meaning |
|-----|---------|
| 0.5 | chance ranking (no skill) |
| 1.0 | perfect ranking |
| 0.9 | 90% of random (pos, neg) pairs are ordered correctly |

**PR-AUC** (average precision) is reported as because positives are rare -> its no-skill baseline equals prevalence, not 0.5.

$$
\mathrm{TPR}(t) = P\bigl(s \ge t \mid y=1\bigr)
\qquad
\mathrm{FPR}(t) = P\bigl(s \ge t \mid y=0\bigr)
$$

- **TPR** (true positive rate / sensitivity): fraction of real V2R bins we catch
- **FPR** (false positive rate): fraction of empty bins we falsely call V2R

Sweep $t$ from high → low. Plot TPR against FPR → the **ROC curve**.

---

#### spatial cross-val' -> pseudoreplication i
- neighbouring 1mb = correlated -> inflates AUC
**To address**
  - 10Mb block groupkfold -> contigous blocks never straddle train/test
  - train 1 chr at a time

