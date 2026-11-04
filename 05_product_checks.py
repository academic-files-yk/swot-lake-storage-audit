"""Reproduce the LakeSP ds1 equations from the published inputs and summarise
field availability (all 455 lakes with Version D data, not only gauged ones)."""
import glob, json
import numpy as np, pandas as pd
import importlib
from common import *
mt = importlib.import_module("04_match")

lakes = [p.split("/")[-1][:-4] for p in glob.glob(str(SWOT / "D" / "*.csv"))]
d = pd.concat([x for x in (mt.load_swot("D", l) for l in lakes) if x is not None and len(x)])
res = {"n_lakes": int(d.lake_id.nunique()), "n_obs": int(len(d)),
       "t_first": str(d.time.min()), "t_last": str(d.time.max())}
v = d.dropna(subset=["wse", "area_total"])
res["n_obs_valid"] = int(len(v))
for c in ("ds1_l", "ds1_q", "ds2_l", "ds2_q"):
    res[f"frac_{c}"] = float(v[c].notna().mean())
    res[f"lakes_{c}"] = int(v.loc[v[c].notna(), "lake_id"].nunique())
res["frac_ds1_partial0"] = float(v.loc[v.partial_f == 0, "ds1_l"].notna().mean())
res["frac_ds1_partial1"] = float(v.loc[v.partial_f == 1, "ds1_l"].notna().mean())

# equations (LakeSP ATBD 3.6.5.3): dV relative to the PLD reference, minus p_ds_t0
h, A, hr, Ar = v.wse, v.area_total, v.p_ref_wse, v.p_ref_area
lin = 0.5 * (A + Ar) * (h - hr) / 1e3 - v.p_ds_t0
quad = (h - hr) * (A + Ar + np.sqrt(A * Ar)) / 3 / 1e3 - v.p_ds_t0
ok = v.ds1_l.notna()
el = (v.ds1_l - lin)[ok].abs(); eq = (v.ds1_q - quad)[ok].abs()
res["repro_lin_p99_km3"] = float(el.quantile(.99)); res["repro_lin_frac_1e-6"] = float((el <= 1e-6).mean())
# tolerance from the rounding of the served values (wse and p_ref_wse to 1 mm, areas and ds to 1e-6)
tol = 0.5 * (A + Ar) * 0.001 / 1e3 + 0.5 * np.abs(h - hr) * 1e-6 / 1e3 + 1.5e-6
res["repro_lin_frac_tol"] = float((el <= tol[ok]).mean()); res["repro_quad_frac_tol"] = float((eq <= tol[ok]).mean())
res["repro_n"] = int(ok.sum())
res["repro_quad_p99_km3"] = float(eq.quantile(.99)); res["repro_quad_frac_1e-6"] = float((eq <= 1e-6).mean())
res["ds1l_minus_ds1q_med_rel"] = float(((v.ds1_l - v.ds1_q).abs() / v.ds1_l.abs().clip(1e-9))[ok].median())

# uncertainty: reproduce ds1_l_u by first-order propagation of wse_u and area_tot_u
u_lin = 0.5 * np.sqrt(((A + Ar) * v.wse_u) ** 2 + ((h - hr) * v.area_tot_u) ** 2) / 1e3
ru = (v.ds1_l_u / u_lin)[ok]
res["u_ratio_med"] = float(ru.median()); res["u_ratio_p05"] = float(ru.quantile(.05)); res["u_ratio_p95"] = float(ru.quantile(.95))
res["wse_u_median_m"] = float(v.wse_u.median()); res["wse_u_p95_m"] = float(v.wse_u.quantile(.95))
res["wse_u_frac_le_1cm"] = float((v.wse_u <= 0.0105).mean())
res["area_u_rel_median"] = float((v.area_tot_u / v.area_total).median())
res["area_u_rel_p95"] = float((v.area_tot_u / v.area_total).quantile(.95))
res["ds1_u_rel_median"] = float((v.ds1_l_u / v.p_ref_area / 1e-3)[ok].median())  # km3 per km2 -> m equivalent
res["quality_f_counts"] = v.quality_f.value_counts().sort_index().to_dict()
res["ref_dates_unique_per_lake_max"] = int(v.groupby("lake_id").p_date_t0.nunique().max())
res["p_ref_wse_unique_per_lake_max"] = int(v.groupby("lake_id").p_ref_wse.nunique().max())
json.dump({k: (int(x) if isinstance(x, (np.integer,)) else x) for k, x in res.items()},
          open(OUT / "product_checks.json", "w"), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
