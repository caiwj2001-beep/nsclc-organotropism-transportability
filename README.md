# nsclc-organotropism-transportability

Analysis code for the manuscript:

> **Organotropism genomics and metastatic burden in non-small cell lung cancer: a 6,516-patient MSK-CHORD analysis**

This repository contains **code only**. No patient-level or derived data are redistributed here.

## Data source

The analysis uses the public MSK-CHORD release, available through cBioPortal:

- Study ID: `msk_chord_2024`
- https://www.cbioportal.org/study/summary?id=msk_chord_2024

Required files (download from the study's *Download* tab):

| File | Used for |
|---|---|
| `data_clinical_patient.csv` | organ annotations, survival, age, sex, smoking, stage |
| `data_clinical_sample.csv` | TMB, FGA, purity, MSI, sample type, panel |
| `data_mutations.csv` | mutation calls (Entrez Gene ID; map to symbols) |
| `data_cna_long.csv` | discrete copy-number calls |

A gene annotation table (`_gene_map.csv`, columns `entrezGeneId,hugoGeneSymbol,fullName`) is needed to map
`entrezGeneId` to HUGO symbols, because `data_mutations.csv` carries no gene symbol column.

## Environment

Python 3.12 with `pandas`, `numpy`, `scipy`, `statsmodels`.
Figures were produced in R 4.6 with `ggplot2`, `gridExtra`, `dplyr`, `scales`.

```
python -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Run order

| Step | Script | Purpose |
|---|---|---|
| 1 | `00_power_calculations.py` | minimum detectable interaction ratios |
| 2 | `01_panel_core_genes.py` | genes observed across MSK-IMPACT 341/410/468/505 (335-gene core) |
| 3 | `02_replicate_lengel_2023.py` | re-derivation of the published LUAD map |
| 4 | `02b_lymphnode_models.py` | nodal models used in the replication step |
| 5 | `03_cohort_landscape.py` | cohort construction, baseline, organ distribution, mutational landscape |
| 6 | `04_cna_matrices.py` | patient-level CNA matrices |
| 7 | `05_single_organ_tropism.py` | single-organ contrast (burden-free tropism) |
| 8 | `06_burden_conditioned_screen.py` | burden-conditioned univariate screening, BH-FDR |
| 9 | `07_multivariable_bootstrap.py` | multivariable models, bootstrap, interaction tests |
| 10 | `08_burden_and_survival.py` | organ-count models, CIN, overall survival |
| 11 | `09_sensitivity_and_power.py` | four sensitivity scenarios, interaction power simulation |
| 12 | `10_figure_data.py` | assembles the plotting tables |
| 13 | `11_figures.R` | builds Figures 1-3 |

Each script reads the raw files from a `DATA_DIR` constant at the top; point it at your local copy before running.
Outputs are written to the working directory and are not included in the repository.

## Notes and pitfalls

- `data_clinical_sample.csv` has no `patientId` column; the sample-to-patient mapping is taken from
  `data_mutations.csv`.
- `SAMPLE_CLASS` is uniformly "Tumor" in this release; primary versus metastatic is `SAMPLE_TYPE`.
- Free-text `METASTATIC_SITE` is available for only about a third of samples and is used for validation,
  not as the primary phenotype.
- `PDL1_POSITIVE` is populated for ~3% of patients and is not used.
- Organ flags are patient-level binaries without lesion counts or dates, so all estimates are cross-sectional.
- Panel generation is a confounder and is included as a covariate in every multivariable model.

## License

Code is released under the MIT License. The underlying data remain subject to the terms of the
MSK-CHORD release and cBioPortal.
