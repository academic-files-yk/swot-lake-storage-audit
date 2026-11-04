#!/bin/sh
# Reproduce the full analysis: data download -> matching -> evaluation -> figures -> numbers and tables.
# Large intermediate files are written to ~/swot_p37_work (see common.py).
set -e
PY=${PY:-python}
$PY 00_fetch_mapping.py          # gauge-to-PLD mapping (Levenson et al. data release)
$PY 01_download_lakesp.py D C    # SWOT LakeSP Prior time series, Version D and C (Hydrocron)
$PY 02_download_insitu.py usgs usbr
$PY 02b_download_benchmark.py    # archived benchmark files (CDEC records, cross-check)
$PY 02c_twdb.py                  # Texas Water Development Board
$PY 03_build_insitu.py           # one daily in-situ table per lake (despiked)
$PY 03b_crosscheck.py
$PY 04_match.py                  # SWOT observation <-> in-situ value
$PY 05_product_checks.py         # reproduction of ds1_l / ds1_q and field availability
$PY 06_storage_eval.py           # pairwise storage-change evaluation and coverage
$PY 07_obs_errors.py             # WSE and area errors, Version C vs D
$PY 08_error_model.py            # empirical error model, cross-validation, detection
$PY 09_sensitivity.py
$PY 10_figures.py
$PY 11_numbers.py                # every number and table quoted in the paper
