# Design decisions (fixed before results were inspected)

Audit of the official SWOT LakeSP storage-change attributes (`ds1_l`, `ds1_q`, `ds2_*`)
and their stated uncertainties against independent reservoir records.

| ID | Decision |
|---|---|
| D1 | Sample: PLD lakes in the gauge mapping of Levenson et al. (2026 data release) with in-situ storage from USBR RISE, USGS NWIS, TWDB, or CDEC (CDEC from the archived benchmark files, period ends Oct 2025). A lake enters the storage analysis with >= 10 matched valid observations. |
| D2 | One in-situ source per lake, priority TWDB > USBR > USGS > CDEC. Units converted to m and km3 (1 acre-foot = 1233.48 m3). |
| D3 | Daily in-situ values are placed at 12:00 local mean solar time (UTC - lon/15 h) and linearly interpolated to the SWOT time; observations with no in-situ value within 1.5 days on both sides are dropped. |
| D3b | Added after inspecting the example series (Fig. 2) and before any result was recomputed: isolated daily in-situ values that depart from a centred 7-day running median by more than 10 NMAD of the day-to-day differences are removed (storage, stage and area separately). |
| D4 | Valid SWOT observation: `wse`, `area_total`, `ds1_l` not fill; `ice_clim_f` = 0. Filters: F_all (no flag filter), F_q01 (`quality_f` <= 1), F_q0 (`quality_f` = 0, the product's "good"). Main results use F_q0; the others are sensitivity runs. |
| D5 | Storage change is evaluated on consecutive valid observation pairs of each lake (time-ordered after filtering). A second view uses the full anomaly series (per-lake median offset removed). |
| D6 | Errors: e = dS_product - dS_insitu (km3). Per-lake summaries: RMSE, and RMSE normalised by the SD of in-situ dS over the same pairs (nRMSE). Cross-lake medians with 95 % bootstrap intervals (lakes resampled, 2000 draws, fixed seed). |
| D7 | Coverage of stated uncertainty: z = e / sqrt(u_i^2 + u_j^2). Report P(|z| <= 1) (nominal 68.3 %), P(|z| <= 1.96) (nominal 95 %), and the inflation factor k = NMAD(z). Same test for WSE (pair differences, `wse_u`) and area (`area_tot_u`, where in-situ area exists). |
| D8 | Error decomposition with in-situ hypsometry V(h) (monotone PCHIP through 1-cm-binned medians of daily stage-storage pairs): level-only dS = V(h_g,i + dh_SWOT) - V(h_g,i); area-only dS = 0.5 (A_i + A_j) dh_gauge with SWOT areas. |
| D9 | Alternative SWOT-only estimators on the same pairs: (a) pairwise trapezoid 0.5 (A_i + A_j)(h_j - h_i); (b) lake hypsometry fitted to SWOT (wse, area) pairs with a Theil-Sen line, integrated analytically (what the unpopulated `ds2` intends); (c) SWOT WSE with in-situ hypsometry (reference case, not SWOT-only). |
| D10 | Empirical error model from first-order propagation of the product equation with per-observation WSE error sigma_h and relative area error r; sigma_h and r are estimated on one half of the lakes and coverage tested on the other half (two-fold split by lake, fixed seed, repeated 100 times). |
| D11 | Version C vs D: same lake, cycle and pass; paired comparison of WSE anomaly error and area error (Wilcoxon signed-rank). |
| D12 | Software: Python 3.13, NumPy, pandas, SciPy, statsmodels; versions in requirements.txt. |
