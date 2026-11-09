# Reproducibility and manuscript-concordance audit

## What is executable

Cellular summary statistics; Hill fits and conditional IC50; Bliss; interaction model; concentration holdout; linear benchmarks; exploratory paired biological contrasts with BH and sign-flip; LOBO single-agent Hill; HPrEC MTT summary. Optional GSE83654 rank-based pathway analysis when third-party inputs are obtained separately.

## Not independently verifiable from distributed inputs

1. The original well-level absorbance, blank subtraction and normalization, and whether the biological mean values reproduce the laboratory source files exactly.
2. The original third-party expression-scale interpretation and dye-swap processing for GSE83654; the pipeline follows the previously used descriptive rank workflow.
3. Exact manuscript parameter/metric concordance, including optimizer bounds, bootstrap resampling and confidence intervals: not certified. Numerical results may differ from the manuscript until a side-by-side reconciliation is completed.
4. Source licenses and authorization to distribute the derived biological means must be confirmed by the depositor.
5. Experimental biological n=3; exact two-sided sign-flip tests cannot achieve p<0.05.

## Release gate

Before publishing Zenodo as a *fully reproducing* archive, run `python run_all.py --with-geo`, compare each resulting value with manuscript tables and text, document accepted tolerances and resolve mismatches. Preserve provenance of the external GEO and MSigDB releases (SHA-256 manifest). Until then label this an executable numerical reanalysis package, not a verified exact reproduction.

## Tested numerical concordance (2026-11-08)

The included processed DU-145 input was independently recomputed from the locally available source normalized measurements before packaging: all 360 biological means matched to floating-point precision (maximum absolute difference < 1.5e-14 percentage points). **The source well-level file is not distributed.** This does not validate the upstream laboratory normalization.

Running the current grouped holdout code gives:

| Model | Computed MAE | Computed RMSE | Manuscript MAE | Manuscript RMSE |
|---|---:|---:|---:|---:|
| Hill–Bliss | 5.534594 | 6.693708 | 5.54 | 6.69 |
| Hill–interaction | 1.535542 | 2.320794 | 1.54 | 2.32 |
| Direct linear combination | 2.255747 | 2.885626 | 2.26 | 2.89 |

The Hill–Bliss MAE rounds to **5.53**, while the manuscript reports **5.54**; this small discrepancy must be resolved before asserting exact numerical concordance. The other listed metrics match at two decimals. Single-agent Hill LOBO mean RMSE = **2.3530** and mean MAE = **2.0122**, consistent with manuscript values 2.35 and 2.01.

The optional GSE83654 pathway code was also executed locally against separately available third-party files (not packaged): 21,061 mapped genes and 1,135 eligible sets (50 Hallmark and 1,085 Reactome), matching the manuscript's reported counts. This does not resolve the GEO scale/dye-swap provenance limitation.

No manuscript supplementary CSV tables or figures are included in this archive.
