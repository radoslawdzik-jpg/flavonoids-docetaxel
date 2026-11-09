# Reproducible prostate-cell docetaxel response analyses (Zenodo package)

This archive contains **analysis code and processed biological-experiment-level response data only**. It contains **no original well-level measurements, absorbance readings, GEO source files, manuscript supplementary tables, or figures**. The outputs are computed when the scripts run, rather than being bundled as finished tables.

## Scope and scientific interpretation

DU-145: apigenin/quercetin at 0, 10, 25, 50, 100 µM, with/without 10 nM docetaxel, 24/48 h, XTT/NR/SRB. HPrEC: six measured MTT groups at 24 h. Each derived input is the average of **four technical replicates within one of three independent biological experiments** (biological n=3). Source absorbance and upstream normalization cannot be audited from this archive.

**Do not claim that this archive reproduces every published number without a concordance audit.** It is an executable numerical reanalysis of all supported cellular components and an optional public-GEO workflow. The original phase scripts had different optimization constraints in some versions; exact agreement with manuscript fit parameters is not established. See `docs/REPRODUCIBILITY_AUDIT.md`.

## One-command cellular analysis

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run_all.py
python -m pytest -q
```

On Windows use `.venv\Scripts\activate`. Python 3.10+ is recommended. Outputs are written to `outputs/cellular/`, and a run manifest with input SHA-256 hashes is saved in `outputs/reproduction_manifest.json`.

## Reproduced calculations

- DU-145 and HPrEC mean/SD based on n=3 biological experiments, not n=12 independent observations.
- Hill concentration-response parameters, half-amplitude EC50, and absolute 50%-control IC50 only when reached within 0–100 µM.
- Experiment-matched Bliss expectation and observed-minus-expected deviations.
- Phenomenological Hill–Bliss interaction, residual penalty `0.15*A`, equivalent to `lambda^2*A^2` for `lambda=0.15`.
- Held-out concentration predictions at 10, 25 and 50 µM, with Hill–Bliss, Hill-interaction, linear Bliss and direct combination interpolation benchmarks.
- Exploratory paired biological differences, within-family and global BH adjustments, and exact two-sided sign-flip sensitivity tests (minimum possible p=0.25 for n=3).
- Single-agent Hill leave-one-biological-experiment-out (LOBO) analysis, individual biological-experiment Hill/IC50 fits, and 400-iteration biological-cluster bootstrap diagnostics.

**Loewe additivity is not computed:** only one nonzero docetaxel concentration is available. Internal concentration holdout is not external validation.

## Optional GSE83654 transcriptomics (third-party source not redistributed)

For the transcriptomic component, separately obtain **GSE83654_series_matrix.txt.gz** and **GPL13607-20416.txt** from GEO, and **Hallmark/Reactome v2025.1.Hs gene sets** from MSigDB, under their respective access and licensing terms. Place these four files in `data/external/` with the names below:

```
data/external/GSE83654_series_matrix.txt.gz
data/external/GPL13607-20416.txt
data/external/h.all.v2025.1.Hs.symbols.gmt
data/external/c2.cp.reactome.v2025.1.Hs.symbols.txt
```

Then run:

```bash
python run_all.py --with-geo
```

Alternatively, to fetch the GEO series matrix directly from the public repository:

```bash
python -m src.transcriptomics --download-geo
```

This command still requires the original GPL platform table and licensed gene-set files. No source data are copied into the Zenodo archive. **The GEO analysis is descriptive**, not a replicated differential-expression or inferential enrichment analysis. Platform-specific dye-swap handling and the original expression scale remain a provenance limitation. Avoid calling the values verified log2 fold changes.

## Manuscript concordance

`outputs/holdout_manuscript_concordance.csv` checks three main model metrics against rounded values stated in the manuscript. Hill–interaction and direct linear combination metrics match at two decimals. Hill–Bliss MAE is 5.534594 (rounds to 5.53), whereas the manuscript reports 5.54; see `docs/REPRODUCIBILITY_AUDIT.md`. **This minor mismatch remains open.**

## Output structure

- `outputs/cellular/`: derived numerical CSV files for summaries, fits, interaction and holdout metrics, paired tests, sign-flip sensitivity, LOBO, HPrEC.
- `outputs/transcriptomics/`: gene-level derived time course, pathway rank-biserial results and QC (only if external files supplied).
- `outputs/reproduction_manifest.json`: runtime versions and SHA-256 checksums.

**No images are included or generated.** The code does not export the manuscript's supplementary CSV tables; any runtime numerical CSVs are analysis intermediates/results, not packaged supplementary files.

## Input permissions and restrictions

The processed cell-line data are derived from normalized source measurements. The depositor must confirm redistribution rights for these derived data; processing alone does not automatically remove contractual or copyright restrictions. If such permission is not available, remove the derived input CSV files and revise the reproducibility statement. Do not publish this package without checking this condition.

## Suggested availability statement (only after rights and concordance review)

"Analysis scripts and derived biological-experiment-level normalized responses are deposited in Zenodo. Original well-level measurements and third-party GEO/MSigDB files are not redistributed. The scripts reproduce the downstream numerical analyses from the derived inputs and support independent reconstruction of the exploratory public-transcriptome analysis when the original third-party inputs are obtained separately."
