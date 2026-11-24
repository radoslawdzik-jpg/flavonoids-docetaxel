# Reproducible analyses of apigenin, quercetin and docetaxel responses in prostate cells

**Code and processed-data repository for:** *Apigenin and Quercetin as Potential Modulators of Docetaxel Response in Prostate Cancer Cells: Experimental Evaluation and Concentration–Response Modeling* (manuscript version dated 22 November 2026; publication metadata to be updated).

This repository provides executable Python scripts, experiment-level normalized cellular response data, numerical analysis outputs, manuscript-oriented derived CSV tables and exploratory processed transcriptomic results. **No figures or original well-level measurements are distributed.** It is a downstream numerical reanalysis and transparency package; **exact reproduction of every manuscript value has not been certified**. See [reproducibility audit](docs/REPRODUCIBILITY_AUDIT.md) and [data provenance](docs/DATA_PROVENANCE_AND_LIMITS.md).

## Study scope and experimental units

- **DU-145 prostate cancer cells:** apigenin or quercetin at 0, 10, 25, 50 and 100 µM, with or without docetaxel (10 nM), assessed at 24 and 48 hours using XTT, neutral red (NR) and sulforhodamine B (SRB) assays.
- **HPrEC cells:** selected conditions evaluated using an MTT assay at 24 hours. The HPrEC MTT readout should not be interpreted as directly interchangeable with the DU-145 assays.
- **Replication:** three independent biological experiments per condition, with four technical replicates averaged within each experiment. Analyses use biological **n = 3**, not 12 independent samples.

The distributed normalized biological means do not permit independent verification of original absorbance measurements, blank subtraction or upstream normalization.

## Quick start

Run from the repository root with Python 3.10 or newer:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run_all.py
python -m pytest -q --import-mode=importlib tests/ test_reproduction.py
```

On Windows, activate the environment with `.venv\Scripts\activate`. The main command reads `data/processed/du145_biological_means.csv` and `data/processed/hprec_biological_means.csv` and writes numerical results to `outputs/`. Previously generated CSV/JSON results are also retained in `outputs/`; rerunning the pipeline can overwrite those files. **No figures are generated.** The test suite was reported as **7 passed** on 24 November 2026 in the integration environment; this is not a manuscript-concordance certification.

## Repository contents

| Path | Contents |
| --- | --- |
| `src/`, `run_all.py` | Cellular numerical analysis, validation, optional transcriptomics and generic overlap utilities |
| `tests/`, `test_reproduction.py`, `pytest.ini` | Unit and reproducibility-oriented checks |
| `data/processed/*biological_means.csv` | Normalized response data at biological-experiment level |
| `data/processed/manuscript_derived/` | Previously prepared manuscript-oriented CSV tables: S1, S2, S3 (two files), S5, S6, S7, S8 and S9 |
| `data/processed/transcriptomics/` | Exploratory GSE83654 gene-level time courses, pathway rank effects and QC summaries |
| `data/external/README.md` | Instructions for external inputs; original third-party datasets are excluded |
| `outputs/cellular/` | Previously generated numerical summaries, fitted parameters, holdout predictions/metrics, LOBO and other CSV results |
| `outputs/holdout_manuscript_concordance.csv` | Limited numerical comparison with selected manuscript metrics |
| `outputs/reproduction_manifest.json` | Input hashes and runtime metadata from a prior run |
| `docs/` | Audit and data provenance/redistribution limitations |
| `CITATION.cff`, `LICENSE` | Citation metadata and repository code licence; see data-rights note below |

The manuscript-derived tables are **not a verified complete final supplementary deposit**. Confirmed processed files for S4, S10 and S11 were not available in the package used for integration. Some CSVs may reflect different stages of manuscript preparation.

## Cellular analyses

The pipeline supports:

1. Biological-experiment-level mean and SD summaries for DU-145 and HPrEC.
2. Constrained Hill concentration–response fits, distinguishing the half-amplitude concentration from the absolute 50%-of-control IC50; the latter is reported only when the threshold is attained within the tested 0–100 µM range.
3. Experiment-matched Bliss independence expectations and observed-minus-expected deviations. These deviations are descriptive and **do not establish pharmacological synergy**.
4. Phenomenological Hill–Bliss interaction modelling and concentration-holdout comparisons against alternative numerical benchmarks.
5. Single-agent Hill leave-one-biological-experiment-out (LOBO) analyses, biological-experiment fits and biological-cluster bootstrap diagnostics.
6. Exploratory paired biological contrasts, Benjamini–Hochberg adjustments and exact two-sided sign-flip sensitivity tests. With n = 3 paired biological experiments, the minimum attainable two-sided sign-flip p-value is 0.25.

Only one nonzero docetaxel concentration was tested; **Loewe additivity is not evaluated**, and internal holdout/LOBO assessments are **not external validation**.

## Optional GEO transcriptomics reconstruction

Exploratory processed results for **GSE83654** are supplied under `data/processed/transcriptomics/`. To reconstruct the optional analysis, obtain the following separately, respecting the original source terms:

```text
data/external/GSE83654_series_matrix.txt.gz
data/external/GPL13607-20416.txt
data/external/h.all.v2025.1.Hs.symbols.gmt
data/external/c2.cp.reactome.v2025.1.Hs.symbols.txt
```

The first two are GEO series/platform files; the latter two are Hallmark and Reactome gene-set resources from the specified MSigDB release. Then run:

```bash
python run_all.py --with-geo
```

The earlier exploratory workflow described 21,061 mapped genes, 50 Hallmark sets and 1,085 Reactome sets. The GEO dye-swap observations are technical measurements, not independent biological replicates. The analysis is descriptive and rank-based; do not interpret its values as independently verified log2 fold changes, replicated differential-expression tests or molecular validation of the flavonoid–docetaxel combinations.

### Relationship between manuscript and processed transcriptomic tables

`Table_S9_Gene_level_timecourse.csv` and `gene_timecourse_derived.csv` have the same 21,061 gene identifiers, but use different time-point column names. `Table_S8_Reactome_rank_biserial.csv` and `reactome_timecourse_rank_effects.csv` have the same 1,085 pathway identifiers; the manuscript table additionally contains `set_n` and `overlap_n`. A local integration comparison found maximum absolute numeric differences of approximately **5 × 10⁻⁹** across matched time-point columns, consistent with export rounding. These are alternative representations of overlapping results, **not independent validation cohorts**.

## Optional CTD/drug–gene overlap extension

`src/overlap.py` provides a generic drug–gene overlap workflow. The manuscript discusses a 22-gene CTD intersection, but **the verified CTD association export needed to regenerate that exact intersection is not included**. An authorized association CSV must be supplied separately, with columns `drug,gene_symbol,relation_type,source,evidence_id`, along with gene-set files as applicable.

```bash
python -m src.overlap \
  --associations /path/to/associations.csv \
  --expression data/processed/transcriptomics/gene_timecourse_derived.csv \
  --hallmark /path/to/hallmark.gmt \
  --reactome /path/to/reactome.gmt \
  --out outputs/overlap
```

This example illustrates use of the utility; it does not certify reproduction of the manuscript's 22-gene result.

## Reproducibility and manuscript concordance

A prior audit reported the following grouped holdout metrics:

| Model | Pipeline MAE | Pipeline RMSE | Manuscript MAE | Manuscript RMSE |
| --- | ---: | ---: | ---: | ---: |
| Hill–Bliss | 5.534594 | 6.693708 | 5.54 | 6.69 |
| Hill–interaction | 1.535542 | 2.320794 | 1.54 | 2.32 |
| Direct linear combination | 2.255747 | 2.885626 | 2.26 | 2.89 |

The **Hill–Bliss MAE rounds to 5.53**, whereas the manuscript states **5.54**. This discrepancy remains unresolved. Other numerical results, optimization settings, bootstrap procedures and final manuscript tables still require a systematic concordance audit. The reported tests demonstrate that the tested code paths run successfully, not that every published number is reproduced exactly.

Before a definitive Zenodo release, verify the final manuscript tables, versioned external inputs, output tolerances and release metadata. See `docs/REPRODUCIBILITY_AUDIT.md` for the audit record.

## Data provenance, rights and exclusions

The original well-level measurements, absorbance files, GEO source matrices, GPL platform annotation, MSigDB gene-set files and CTD association export are **not redistributed**. Processed-data redistribution rights must be confirmed by the depositor; the code licence does not automatically license third-party materials or experimental data. Consult `docs/DATA_PROVENANCE_AND_LIMITS.md` before public release.

No PNG, JPEG, SVG, TIFF, PDF or other figure files are intended for this repository. Numerical CSV/JSON outputs are intentionally retained. Python caches, temporary files and local environments should remain excluded via `.gitignore`.

## Citation and release

Use `CITATION.cff` for citation metadata, and update the manuscript bibliographic details, repository release version and Zenodo DOI **after** the corresponding release has been published. Do not cite a DOI that has not yet been assigned. The GitHub repository and any subsequent Zenodo archive should identify the same tagged code/data snapshot.
