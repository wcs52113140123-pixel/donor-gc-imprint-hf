# figures_v06: figure scripts for manuscript v0.6

v0.6 re-lays out every figure so that all text is at least 10 pt (journal requirement); data, statistics and printed numbers are unchanged from v0.5.
Panel letters are unchanged. A few in-panel annotations were moved into the figure legends to make room.

## Run

The scripts read the same frozen inputs as v0.5 (`../figures_v05/inputs/`; override with `JBTS_INPUTS`). Run each script from an empty working directory; outputs (TIFF/PNG/PDF, number registry, source-data workbook) are written there.

```
cd /tmp/out_fig3 && HAHN2021_DIR=/path/to/hahn2021 python /path/to/repo/figures_v06/fig3.py
```

- `fig3.py` panel H needs `ClinicalCharacteristicsRNAseqSamples.csv` from the supplementary files of Hahn et al. 2021 (not redistributed); set `HAHN2021_DIR` to its folder.
- `fig7.py` and `ci.py` skip the check of printed numbers against the locked manuscript text when that text is absent (it is not distributed) and record "not checked" in the registry.
- `jbts_style.py` sets the 10-pt floor and provides `min_font_audit(fig, 10)`, which every script runs before saving.

## Changes from the internal scripts

`PUBLIC_PATCH.diff`: absolute input paths replaced with paths relative to this folder, plus the two release changes above (Hahn table location; optional manuscript text). No plotting or statistical code was changed.

## Reproduction check

`reproduction_check.csv`: all 9 scripts were run in a clean checkout. Every number registry was identical to the one submitted with manuscript v0.6. TIFFs were pixel-identical for 5 of 9 figures; Figures 3, 4, 7 and the Central Illustration differed in 0.9-1.3% of pixels, all in text anti-aliasing from the rendering environment.

Fonts: Arial is used when installed; otherwise the metric-compatible Liberation Sans (the submitted files used Liberation Sans).
