# Gapfilling & Final Polish

---

### Gapfilling with [TGS-GapCloser](https://github.com/BGI-Qingdao/TGS-GapCloser)
```
- Using length filtered **long-reads** (ONT)
- The curated assembly's residual scaffold gaps (N-runs introduced during scaffolding) are closed here.
- **Justification**: Scaffolding joins contigs across gaps of estimated size; filling those gaps with real long-read sequence converts ambiguous N-runs into called bases, improving contiguity (scaftig N50) and completeness without altering the validated scaffold order/orientation. A consolidated alternative, `tgs-gap_nextpolish.sh`, performs the equivalent step with TGS-GapCloser (Xu et al., 2020), which is purpose-built for closing gaps in large genomes from low-coverage error-prone long reads, followed immediately by the Step 9 polish.
```

### Final Polish with [NextPolish](https://github.com/Nextomics/NextPolish)
```
- using **short-reads only**
- fixes base errors (SNV/Indel) left in the genome following prior steps in assembly, or residual errors from reads
- A final NextPolish round is run on the curated, gap-filled assembly. This round uses only short reads (sgs_option), with minimap2 -x map-pb/map-ont for the long-read alignment.
- **Justification**: Gap-filling inserts new long-read sequence at gap junctions, and curation creates new contig joins; both introduce fresh, unpolished bases and junction indels. A final polish corrects these, with short reads supplying per-base accuracy and long reads disambiguating repeat-region context, yielding the best achievable consensus QV before release.
```
