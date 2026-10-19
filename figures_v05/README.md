# figures_v05 — regenerate the v0.5 figures

This folder rebuilds every figure of manuscript v0.5 (Figures 1–7, the Central Illustration and Supplemental Figure S1) from the frozen input bundles in `inputs/`. Each script writes its outputs (PDF, PNG, 600-dpi LZW TIFF, a source-data workbook and a number registry) into the **current working directory**.

## Requirements
Python 3.10+ with numpy, pandas, scipy, matplotlib (tested with 3.11.1 / pandas 2.3.3) and openpyxl. Fonts: Arial if installed, otherwise Liberation Sans, otherwise DejaVu Sans (chosen automatically by `jbts_style.py`). With a different font, the layout can shift slightly; numbers do not change.

## Files
| File | Builds | Inputs used (from `inputs/`) |
|---|---|---|
| `fig1.py` | Figure 1 (design, cohorts, analysis plans, test outcomes) | figure_inputs_v05.tar.gz, figure_inputs_v05_perdonor.tar.gz, Supplemental_Tables_S1-S14.xlsx |
| `fig2.py` | Figure 2 (in vitro GC reference) | figure_inputs_v05.tar.gz, Supplemental_Tables_S1-S14.xlsx |
| `fig3.py` | Figure 3 (GC-like imprint by cohort) | figure_inputs_v05.tar.gz, figure_inputs_v05_perdonor.tar.gz, Supplemental_Tables_S1-S14.xlsx; panel H also needs the Hahn et al. 2021 clinical table via `HAHN2021_DIR` |
| `fig4.py` | Figure 4 (circadian and cAMP scores) | same as fig3 |
| `fig5.py` | Figure 5 (identifiability and gene-wise audit) | same as fig3 |
| `fig6.py` | Figure 6 (literature genes) | same as fig3 + perdonor_gene_expression_GC8_clock10_lvad4.csv.gz |
| `fig7.py` | Figure 7 (macrophage modules, GTEx, aging) | same as fig3 (the locked-text check is skipped in the public release) |
| `ci.py` | Central Illustration | — (locked-text check skipped in the public release) |
| `figS1.py` | Supplemental Figure S1 (module class counts) | figure_inputs_v05.tar.gz, Source_data_v04.xlsx |
| `jbts_style.py` | Shared design system (sizes, palette, 183-mm width), loaded by every script | — |
| `PATH_PATCH.diff` | Only difference from the scripts used to build the submitted figures: absolute paths replaced by package-relative paths | — |
| `reproduction_check.csv` | Locked tables vs reruns of the frozen scoring scripts (10 table-level checks, 39 manuscript claims) | — |

### Inputs
- `figure_inputs_v05.tar.gz`: frozen result tables, configs and processed metadata. See `README_figure_inputs_v05.md` and `MANIFEST_figure_inputs_v05.tsv` inside the archive.
- `figure_inputs_v05_perdonor.tar.gz`: per-donor and per-sample scores dumped by re-running the frozen scripts. Only write-out lines were added (see its `README.md` and `code/dump_patches.diff`).
- `perdonor_gene_expression_GC8_clock10_lvad4.csv.gz`: per-donor expression for the Figure 6 signature genes, from the same frozen run.
- `Supplemental_Tables_S1-S14.xlsx`: the submitted supplemental tables.
- `Source_data_v04.xlsx`: v0.4 source data (module class counts reused by Supplemental Figure S1).
- `Manuscript_JACC-BTS_v0.4_text.md`: the locked results text. Scripts check that printed numbers occur verbatim in it.

## Run
```bash
cd figures_v05
for f in fig1 fig2 fig3 fig4 fig5 fig6 fig7 ci figS1; do
  mkdir -p ../out_$f && (cd ../out_$f && python ../figures_v05/$f.py)
done
```
The scripts unpack the input archives into the working directory (`_inputs/` or `data/`; Figure 7 uses a temporary folder). To keep inputs elsewhere, set `JBTS_INPUTS=/path/to/inputs`. Running time is 3–20 s per figure.

## Verification (v0.5 package build, 2026-10-04)
All 9 scripts were run end-to-end in a clean copy of this folder. Every number registry and every source-data workbook (all sheets) was identical to the submitted files. TIFF/PNG images were visually identical. 1–3% of pixels differed only in anti-aliasing of rendered text, because the rendering environment differed. Figure 7 and the Central Illustration report 0 mismatches against the locked text.
