"""Per-observation WSE and area errors (Version D and C) and the Version C vs D
paired comparison (D11)."""
import json
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
from common import *
from lib import *

def obs_errors(ver, F="F_q0"):
    M = pd.read_csv(WORK / f"matched_{ver}.csv", parse_dates=["time"])
    M = M[(M.ice_clim_f == 0) & M.wse.notna() & M.area_total.notna()]
    if F == "F_q0":
        M = M[M.quality_f == 0]
    M = M.sort_values("time").drop_duplicates(["swot_lake_id", "time"])
    M["e_h_raw"] = M.wse - M.g_stage
    M["e_h"] = M.e_h_raw - M.groupby("swot_lake_id").e_h_raw.transform("median")
    M["e_a_rel"] = M.area_total / M.g_area - 1
    M["e_a_anom"] = M.e_a_rel - M.groupby("swot_lake_id").e_a_rel.transform("median")
    M["z_h"] = M.e_h / M.wse_u
    M["z_a"] = (M.area_total - M.g_area) / M.area_tot_u
    return M

if __name__ == "__main__":
    out = {}
    D = obs_errors("D")
    D.to_csv(WORK / "obs_D_F_q0.csv", index=False)
    h = D.dropna(subset=["e_h"]); h = h[h.groupby("swot_lake_id").e_h.transform("size") >= 10]
    a = D.dropna(subset=["e_a_rel"]); a = a[a.groupby("swot_lake_id").e_a_rel.transform("size") >= 10]
    out["wse"] = {"n_obs": len(h), "n_lakes": int(h.swot_lake_id.nunique()),
                  "nmad_m": nmad(h.e_h), "p68_abs_m": float(h.e_h.abs().quantile(.68)),
                  "wse_u_median_m": float(h.wse_u.median()),
                  "p1": float((h.z_h.abs() <= 1).mean()), "p196": float((h.z_h.abs() <= 1.96).mean())}
    size = pd.cut(h.p_ref_area, [0, 10, 100, 1e6], labels=["<10", "10-100", ">100"])
    out["wse_by_size"] = {str(k): {"n_lakes": int(g.swot_lake_id.nunique()), "nmad_m": nmad(g.e_h),
                          "wse_u_med": float(g.wse_u.median())} for k, g in h.groupby(size, observed=True)}
    lake_bias = a.groupby("swot_lake_id").e_a_rel.median()
    out["area"] = {"n_obs": len(a), "n_lakes": int(a.swot_lake_id.nunique()),
                   "sources": a.groupby("source").swot_lake_id.nunique().to_dict(),
                   "bias_rel_median_obs": float(a.e_a_rel.median()),
                   "bias_rel_lake_median": boot_median(lake_bias.values),
                   "frac_lakes_overestimate": float((lake_bias > 0).mean()),
                   "rand_rel_nmad": nmad(a.e_a_anom),
                   "area_u_rel_median": float((a.area_tot_u / a.area_total).median()),
                   "p1": float((a.z_a.abs() <= 1).mean()), "p196": float((a.z_a.abs() <= 1.96).mean())}
    # dark-water fraction and area error
    dfq = pd.qcut(a.dark_frac, 4)
    out["area_by_darkfrac"] = {str(k): {"n": int(len(g)), "bias_rel_med": float(g.e_a_rel.median()),
                               "nmad": nmad(g.e_a_anom)} for k, g in a.groupby(dfq, observed=True)}
    from scipy.stats import spearmanr
    out["area_darkfrac_spearman"] = spearmanr(a.dark_frac, a.e_a_rel, nan_policy="omit")._asdict()
    # Version C vs D on the same lake, cycle and pass (D11)
    C = obs_errors("C")
    key = ["swot_lake_id", "cycle_id", "pass_id"]
    J = D.merge(C, on=key, suffixes=("_D", "_C"))
    jh = J.dropna(subset=["e_h_D", "e_h_C"])
    ja = J.dropna(subset=["e_a_rel_D", "e_a_rel_C"])
    out["C_vs_D"] = {"n_pairs_wse": len(jh), "n_lakes_wse": int(jh.swot_lake_id.nunique()),
        "abs_eh_med_D": float(jh.e_h_D.abs().median()), "abs_eh_med_C": float(jh.e_h_C.abs().median()),
        "p_wse": float(wilcoxon(jh.e_h_D.abs(), jh.e_h_C.abs()).pvalue),
        "n_pairs_area": len(ja), "n_lakes_area": int(ja.swot_lake_id.nunique()),
        "bias_rel_med_D": float(ja.e_a_rel_D.median()), "bias_rel_med_C": float(ja.e_a_rel_C.median()),
        "abs_ea_med_D": float(ja.e_a_rel_D.abs().median()), "abs_ea_med_C": float(ja.e_a_rel_C.abs().median()),
        "p_area": float(wilcoxon(ja.e_a_rel_D.abs(), ja.e_a_rel_C.abs()).pvalue),
        "frac_D_better_area": float((ja.e_a_rel_D.abs() < ja.e_a_rel_C.abs()).mean()),
        "C_ds1_lakes": int(C.dropna(subset=["ds1_l"]).swot_lake_id.nunique()),
        "C_period": [str(C.time.min()), str(C.time.max())]}
    J.to_csv(WORK / "C_vs_D.csv", index=False)
    json.dump(out, open(OUT / "obs_errors.json", "w"), indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))
