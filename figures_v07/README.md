# figures_v07: figure scripts for manuscript v0.7

v0.7 revises figure wording only. Internal plan codes, version numbers and hashes are replaced by plain test names, and planned tests are reported as met or not met. Figure 1 swaps its plan-history and test-outcome panels for donor age (F) and sex composition (G). Data, statistics and printed numbers are unchanged from v0.6; all text is at least 10 pt.

## Run

Inputs are the frozen bundles in `../figures_v05/inputs/` (override with `JBTS_INPUTS`). Run each script from an empty working directory:

```
cd /tmp/out_fig3 && HAHN2021_DIR=/path/to/hahn2021 python /path/to/repo/figures_v07/fig3.py
```

`fig3.py` panel H needs `ClinicalCharacteristicsRNAseqSamples.csv` from Hahn et al. 2021 (not redistributed). `fig7.py` and `ci.py` skip the check against the manuscript text when it is absent. `PUBLIC_PATCH.diff` lists the only changes from the internal scripts: relative paths plus these two release changes.

## Test names in the manuscript and their plan codes

| Manuscript name | Plan code (01_plans/) |
|---|---|
| Original primary endpoint | version 1.0 endpoint |
| Systemic exposure test | G1 |
| GTEx ventilator comparison (supportive) | G2 |
| GC-independent loss test | G3 |
| Correlation / attenuation tests (GC-induced efferocytic genes) | G4a / G4b |
| Aging equivalence test | G5 |
| Circadian test / cAMP test | C2 / K1 |
| Adjustment positive control | A3 |
| GC plan 1 / GC plan 2 | ANALYSIS_PLAN_v1.2 / v1.3 |

## Reproduction check

`reproduction_check.csv`: all 9 scripts ran in a clean checkout. Every number registry is identical to the submitted v0.7 registry. TIFFs are pixel-identical for Figures 1, 2, 6 and S1; the others differ in 0.6-1.3% of pixels, all in text anti-aliasing.
