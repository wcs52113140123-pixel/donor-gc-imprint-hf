# Public release notes (v0.5.1)

This public release is the code/data package of manuscript v0.5.1 with three changes, made so that it contains only material the authors can redistribute:

1. **Manuscript text removed.** `figures_v05/inputs/Manuscript_JACC-BTS_v0.4_text.md` (used by `fig7.py` and `ci.py` only to check that printed numbers occur verbatim in the locked text) is not distributed. Both scripts now skip that check when the file is absent and record "not checked (manuscript text not distributed)" in the number registry; figures and numbers are unchanged.
2. **Third-party clinical tables removed.** The two clinical tables from the supplementary files of Hahn et al. 2021 (`ClinicalCharacteristicsRNAseqSamples.csv`, `SamplesTable_ALLpatients_sex_race_age.csv`) were removed from `figures_v05/inputs/figure_inputs_v05.tar.gz` (and from its internal manifest). Only `fig3.py` (panel H, medications) needs the first one: download it from the original publication and set `HAHN2021_DIR` to the folder that contains it.
3. **Internal administrative notes removed** (`00_admin/` result notes from the internal plan workflow); no script reads them. The frozen analysis plans remain in `01_plans/`.

`README.md`, `figures_v05/README.md` and `LICENSE` were updated with the release metadata (DOI, repository, copyright holders) and `CITATION.cff` was added. All other files are byte-identical to the package used for the submitted figures; `MANIFEST.sha256` lists the SHA-256 of every file in this release.
