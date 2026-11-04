# swot-lake-storage-audit

Code for the paper

> Kaya, Y. Auditing the official SWOT lake storage-change product and its stated uncertainty at 187 gauged reservoirs. Submitted to *Remote Sensing of Environment*.

The analysis evaluates the storage-change attributes (`ds1_l`, `ds1_q`, `ds2_l`, `ds2_q`) and the stated uncertainties (`wse_u`, `area_tot_u`, `ds1_l_u`) of the SWOT LakeSP Prior product against daily reservoir records.

## Data sources (all public)

| Data | Source | Access |
|---|---|---|
| SWOT LakeSP Prior, Version D (PGD0/PID0) and Version C | NASA PO.DAAC, doi:10.5067/SWOT-LAKESP-D, doi:10.5067/SWOT-LAKESP-2.0 | Hydrocron API, no login |
| Gauge-to-PLD mapping; CDEC daily records | Levenson et al. benchmark, doi:10.5281/zenodo.22258932 | GitHub mirror |
| Reservoir elevation, storage, area | USBR RISE (https://data.usbr.gov) | RISE API |
| Reservoir storage, elevation | USGS NWIS daily values, doi:10.5066/F7P55KJN | waterservices.usgs.gov |
| Reservoir elevation, area, storage | Texas Water Development Board (https://www.waterdatafortexas.org) | CSV per reservoir |
| State boundaries (Fig. 1) | U.S. Census Bureau, cb_2024_us_state_20m | place the shapefile in `~/swot_p37_work/map/` |

## Running

```
pip install -r requirements.txt
sh run_all.sh
```

Downloaded and intermediate data (about 200 MB) are written to `~/swot_p37_work` (set in `common.py`). Results are written to `output/` and figures to `figures/`. Random seeds are fixed (`common.SEED`).

## Scripts

| Script | Purpose | Main outputs |
|---|---|---|
| `00_fetch_mapping.py` | gauge-to-PLD mapping table | `raw/gage_swot_mapping.csv` |
| `01_download_lakesp.py` | LakeSP time series for every mapped lake (Version D and C) | `swot/D`, `swot/C` |
| `02_download_insitu.py`, `02b_download_benchmark.py`, `02c_twdb.py` | in-situ records | `insitu/`, `output/twdb_match.csv` |
| `03_build_insitu.py`, `03b_crosscheck.py` | one daily table per lake, despiking, cross-check | `output/insitu_inventory.csv`, `output/insitu_crosscheck.json` |
| `04_match.py` | time matching (Section 2.4) | `matched_D.csv`, `matched_C.csv` |
| `05_product_checks.py` | reproduction of the product equations (Section 4.1) | `output/product_checks.json` |
| `06_storage_eval.py` | pairwise evaluation, estimators, coverage (Sections 4.2–4.5) | `output/storage_eval.json`, `output/lake_metrics.csv` |
| `07_obs_errors.py` | WSE and area errors, Version C vs D (Section 4.4) | `output/obs_errors.json` |
| `08_error_model.py` | empirical error model and detection (Sections 4.6–4.7) | `output/error_model.json`, `output/error_model_cv.csv`, `output/anomaly_series.csv` |
| `09_sensitivity.py` | sensitivity runs (Table S1) | `output/sensitivity.csv` |
| `10_figures.py` | Figures 1–5 | `figures/` |
| `11_numbers.py` | all numbers and tables of the paper | `output/numbers.json`, `output/tables.json` |

`DESIGN.md` lists the methodological decisions and when they were made.

## Notes

- CDEC could not be reached from the processing network; CDEC records are therefore taken from the archived benchmark files, which end in October 2025.
- Hydrocron reports WSE to 1 mm; the reproduction test in `05_product_checks.py` uses the corresponding rounding tolerance.

## License

MIT (see `LICENSE`).
