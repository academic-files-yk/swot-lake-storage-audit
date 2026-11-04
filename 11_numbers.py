"""Collect every number quoted in the manuscript into output/numbers.json (formatted strings)
and build the table contents (output/tables.json)."""
import json, platform
import numpy as np, pandas as pd, scipy, statsmodels
from scipy.stats import wilcoxon
from common import *
from lib import nmad

J = lambda f: json.load(open(OUT / f))
pc, se, oe, em = J("product_checks.json"), J("storage_eval.json"), J("obs_errors.json"), J("error_model.json")
S = se["F_q0"]
inv = pd.read_csv(OUT / "insitu_inventory.csv")
P = pd.read_csv(WORK / "pairs_F_q0_model.csv")
L = pd.read_csv(OUT / "lake_metrics.csv"); L = L[L["filter"] == "F_q0"]
W = L.pivot(index="swot_lake_id", columns="est", values="nrmse")
N = {}
f2 = lambda x: f"{x:.2f}"; f3 = lambda x: f"{x:.3f}"; pct = lambda x: f"{100 * x:.0f}"; pct1 = lambda x: f"{100 * x:.1f}"
ci = lambda t, f=f2: f"{f(t[0])} ({f(t[1])}–{f(t[2])})"
def pfmt(p):
    return "< 0.001" if p < 0.001 else f"= {p:.3f}" if p < 0.01 else f"= {p:.2f}"

# data
N["n_lakes_all"] = f"{pc['n_lakes']}"; N["n_obs_all"] = f"{pc['n_obs']:,}"
dfmt = lambda x: pd.Timestamp(x).strftime("%-d %B %Y")
N["t_first"] = dfmt(pc["t_first"]); N["t_last"] = dfmt(pc["t_last"])
N["n_map"] = f"{len(pd.read_csv(MAPPING))}"
N["n_insitu_lakes"] = f"{len(inv)}"; N["n_insitu_days"] = f"{inv.n_days.sum():,}"
N["n_despiked"] = f"{inv.n_despiked.sum():,}"; N["pct_despiked"] = f"{100 * inv.n_despiked.sum() / (inv.n_days.sum() + inv.n_despiked.sum()):.1f}"
for s in ("USBR", "USGS", "TWDB", "CDEC"):
    N[f"inv_{s}"] = f"{(inv.source == s).sum()}"; N[f"lakes_{s}"] = f"{S['sources'].get(s, 0)}"
N["n_area_insitu"] = f"{inv.has_area.sum()}"
cc = J("insitu_crosscheck.json")
N["cc_n"], N["cc_frac"] = f"{cc['n_lakes']}", pct(cc["frac_within_0p5pct"])
ar_ = P.groupby("swot_lake_id").pld_area.first()
N["area_min"], N["area_max"] = f"{ar_.min():.2f}", f"{ar_.max():.0f}"
# product checks
N["frac_ds1"] = pct(pc["frac_ds1_l"]); N["lakes_ds1"] = f"{pc['lakes_ds1_l']}"
N["frac_ds1_full"] = pct1(pc["frac_ds1_partial0"]); N["repro_n"] = f"{pc['repro_n']:,}"
N["repro_frac"] = pct1(pc["repro_lin_frac_tol"]); N["u_ratio_med"] = f3(pc["u_ratio_med"])
N["u_ratio_p05"], N["u_ratio_p95"] = f2(pc["u_ratio_p05"]), f2(pc["u_ratio_p95"])
N["wse_u_med_cm"] = f"{100 * pc['wse_u_median_m']:.1f}"; N["area_u_rel_med"] = f"{100 * pc['area_u_rel_median']:.2f}"
N["ds1_lq_rel"] = f"{100 * pc['ds1l_minus_ds1q_med_rel']:.2f}"
qf = pc["quality_f_counts"]; tot = sum(qf.values()); N["frac_q0_all"] = pct(qf["0"] / tot)
# sample
N["n_lakes"] = f"{S['n_lakes']}"; N["n_pairs"] = f"{S['n_pairs']:,}"; N["dt_med"] = f"{S['dt_median_days']:.0f}"
N["n_lakes_lev"] = f"{S['dS_level_only']['n_lakes']}"
for k, e in (("prod", "dS_prod"), ("prodq", "dS_prod_q"), ("trap", "dS_trap"), ("shyp", "dS_shyp"),
             ("lev", "dS_level_only"), ("area", "dS_area_only")):
    N[f"nrmse_{k}"] = ci(S[e]["nrmse"]); N[f"nrmse_{k}_m"] = f2(S[e]["nrmse"][0]); N[f"r_{k}"] = f3(S[e]["r"][0])
for a, b in (("dS_prod", "dS_trap"), ("dS_trap", "dS_shyp"), ("dS_prod", "dS_shyp"), ("dS_level_only", "dS_area_only")):
    x = W[[a, b]].dropna(); t = f"{a[3:]}_{b[3:]}"
    N[f"w_{t}_frac"] = pct(float((x[b] < x[a]).mean())); N[f"w_{t}_p"] = pfmt(wilcoxon(x[a], x[b]).pvalue)
    N[f"w_{t}_n"] = f"{len(x)}"
N["F_q01_prod"], N["F_q01_trap"] = f2(se["F_q01"]["dS_prod"]["nrmse"][0]), f2(se["F_q01"]["dS_trap"]["nrmse"][0])
N["F_all_prod"], N["F_all_trap"] = f2(se["F_all"]["dS_prod"]["nrmse"][0]), f2(se["F_all"]["dS_trap"]["nrmse"][0])
N["F_q01_lakes"], N["F_all_lakes"] = f"{se['F_q01']['n_lakes']}", f"{se['F_all']['n_lakes']}"
# coverage
for k, c in (("dS", S["cov_dS"]), ("dh", S["cov_dh"]), ("dA", S["cov_dA"])):
    N[f"cov_{k}_p1"], N[f"cov_{k}_p2"], N[f"cov_{k}_k"] = pct(c["p1"]), pct(c["p196"]), f"{c['k_nmad']:.0f}"
    N[f"cov_{k}_n"] = f"{c['n']:,}"
N["cov_dA_lakes"] = f"{P.dropna(subset=['dA_g']).swot_lake_id.nunique()}"
N["dh_err_nmad_cm"] = f"{100 * S['dh_err_nmad_m']:.1f}"; N["u_dh_med_cm"] = f"{100 * S['u_dh_median_m']:.1f}"
w = oe["wse"]; a = oe["area"]
N["wse_nmad_cm"] = f"{100 * w['nmad_m']:.1f}"; N["wse_p68_cm"] = f"{100 * w['p68_abs_m']:.1f}"
N["wse_u_med_cal_cm"] = f"{100 * w['wse_u_median_m']:.1f}"; N["wse_n_obs"] = f"{w['n_obs']:,}"; N["wse_n_lakes"] = f"{w['n_lakes']}"
N["wse_p1"], N["wse_p2"] = pct(w["p1"]), pct(w["p196"])
for s_, v in oe["wse_by_size"].items():
    k = {"<10": "small", "10-100": "mid", ">100": "large"}[s_]
    N[f"wse_{k}_cm"], N[f"wse_{k}_u_cm"], N[f"wse_{k}_n"] = f"{100 * v['nmad_m']:.1f}", f"{100 * v['wse_u_med']:.1f}", f"{v['n_lakes']}"
N["area_n_lakes"], N["area_n_obs"] = f"{a['n_lakes']}", f"{a['n_obs']:,}"
N["area_twdb"], N["area_usbr"] = f"{a['sources'].get('TWDB', 0)}", f"{a['sources'].get('USBR', 0)}"
bl = a["bias_rel_lake_median"]; N["area_bias_lake"] = f"{100 * bl[0]:+.1f}% (95% interval {100 * bl[1]:+.1f}% to {100 * bl[2]:+.1f}%)"
N["area_bias_obs"] = f"{100 * a['bias_rel_median_obs']:+.1f}"; N["area_frac_over"] = pct(a["frac_lakes_overestimate"])
N["area_rand_nmad"] = f"{100 * a['rand_rel_nmad']:.1f}"; N["area_u_rel_cal"] = f"{100 * a['area_u_rel_median']:.2f}"
N["area_p1"], N["area_p2"] = pct1(a["p1"]), pct1(a["p196"])
bd = list(oe["area_by_darkfrac"].values())
N["area_nmad_dark_lo"], N["area_nmad_dark_hi"] = f"{100 * bd[0]['nmad']:.1f}", f"{100 * bd[-1]['nmad']:.1f}"
N["area_dark_rho"] = f"{oe['area_darkfrac_spearman']['statistic']:.2f}"
cd = oe["C_vs_D"]
N["cd_n_wse"], N["cd_lakes_wse"] = f"{cd['n_pairs_wse']:,}", f"{cd['n_lakes_wse']}"
N["cd_eh_D"], N["cd_eh_C"], N["cd_p_wse"] = f"{100 * cd['abs_eh_med_D']:.1f}", f"{100 * cd['abs_eh_med_C']:.1f}", pfmt(cd["p_wse"])
N["cd_n_area"], N["cd_lakes_area"] = f"{cd['n_pairs_area']}", f"{cd['n_lakes_area']}"
N["cd_ea_D"], N["cd_ea_C"], N["cd_p_area"] = f"{100 * cd['abs_ea_med_D']:.1f}", f"{100 * cd['abs_ea_med_C']:.1f}", pfmt(cd["p_area"])
N["cd_frac_D_better"] = pct(cd["frac_D_better_area"]); N["cd_end"] = dfmt(cd["C_period"][1])
# error model
N["sig_h_cm"] = f"{100 * em['sigma_h_all']:.1f}"; N["r_pct"] = f"{100 * em['r_all']:.1f}"
N["cv_sig_h"] = f"{100 * em['sigma_h'][1]:.1f}–{100 * em['sigma_h'][2]:.1f}"; N["cv_r"] = f"{100 * em['r'][1]:.1f}–{100 * em['r'][2]:.1f}"
for k in ("prod", "trap", "stated"):
    N[f"cv_{k}_p1"] = f"{pct(em[f'{k}_p1'][0])}% ({pct(em[f'{k}_p1'][1])}–{pct(em[f'{k}_p1'][2])}%)"
    N[f"cv_{k}_p2"] = f"{pct(em[f'{k}_p196'][0])}% ({pct(em[f'{k}_p196'][1])}–{pct(em[f'{k}_p196'][2])}%)"
    N[f"cv_{k}_k"] = f"{em[f'{k}_k'][0]:.2f}" if k != "stated" else f"{em[f'{k}_k'][0]:.0f}"
N["ratio_emp_stated"] = f"{em['ratio_emp_to_stated_median']:.0f}"; N["share_area_var"] = pct(em["share_area_var_prod_median"])
N["z_q95"] = f"{em['model_abs_z_q95']:.1f}"; N["z_q68"] = f"{em['model_abs_z_q68']:.2f}"
zp = (P.dS_prod - P.dS_ins) / P.sig_prod
N["frac_z3"] = pct1(float((zp.abs() > 3).mean()))
d = em["detect_by_dh"]
N["det_5_prod"], N["det_5_trap"], N["det_5_n"] = pct(d["0-0.05"]["sign_prod"]), pct(d["0-0.05"]["sign_trap"]), f"{d['0-0.05']['n']:,}"
N["det_10_prod"], N["det_10_trap"] = pct(d["0.05-0.1"]["sign_prod"]), pct(d["0.05-0.1"]["sign_trap"])
N["det_20_prod"], N["det_20_trap"] = pct(d["0.1-0.2"]["sign_prod"]), pct(d["0.1-0.2"]["sign_trap"])
N["det_50_prod"] = pct(d["0.2-0.5"]["sign_prod"]); N["det_big_prod"] = pct(d["0.5-100"]["sign_prod"])
N["frac_pairs_lt5"] = pct(d["0-0.05"]["n"] / sum(v["n"] for v in d.values()))
N["anom_prod"], N["anom_shyp"] = ci(em["anomaly_nrmse"]), ci(em["anomaly_nrmse_shyp"])
N["anom_r"] = f3(em["anomaly_r"][0]); N["anom_frac_shyp"] = pct(em["anomaly_frac_shyp_better"])
N["x_ref_med"] = f"{em['x_ref_median_m'][0]:.1f}"; N["wse_range_med"] = f"{em['wse_range_median_m'][0]:.1f}"
N["x_over_dh"] = f"{em['pair_x_over_dh_median']:.1f}"
q = list(em["ratio_by_x_quartile"].values())
N["xq1_ratio"], N["xq4_ratio"] = f2(q[0]["ratio_abs_err"]), f2(q[-1]["ratio_abs_err"])
N["xq1_x"], N["xq4_x"] = f2(q[0]["x_med"]), f"{q[-1]['x_med']:.1f}"
N["xq1_worse"], N["xq4_worse"] = pct(q[0]["frac_prod_worse"]), pct(q[-1]["frac_prod_worse"])
# sensitivity
sens = pd.read_csv(OUT / "sensitivity.csv").set_index("case")
for c, k in (("Agency CDEC", "cdec"), ("Agency TWDB", "twdb"), ("Agency USBR", "usbr"), ("Agency USGS", "usgs"),
             ("Without CDEC archive", "nocdec"), ("Pairs <= 11 days apart", "short"), ("Pairs > 11 days apart", "long")):
    N[f"sens_{k}_prod"], N[f"sens_{k}_trap"], N[f"sens_{k}_n"] = f2(sens.loc[c, "dS_prod"]), f2(sens.loc[c, "dS_trap"]), f"{sens.loc[c, 'n_lakes']:.0f}"
N["sens_prod_worse_all"] = "yes" if (sens.dS_prod.iloc[[0] + list(range(3, len(sens)))] > sens.dS_trap.iloc[[0] + list(range(3, len(sens)))]).all() else "no"
N["stated_p1_max"] = pct(sens.stated_p1.max())
N["py"] = platform.python_version(); N["np_v"] = np.__version__; N["pd_v"] = pd.__version__
N["sp_v"] = scipy.__version__; N["sm_v"] = statsmodels.__version__
import matplotlib, geopandas
N["mpl_v"] = matplotlib.__version__; N["gpd_v"] = geopandas.__version__
import re
N = {k: re.sub(r"(?<![\w.])-(?=\d)", "\u2212", v) for k, v in N.items()}   # typographic minus
json.dump(N, open(OUT / "numbers.json", "w"), indent=1, ensure_ascii=False)

# ---------------- tables ----------------
T = {}
src = pd.DataFrame({"Agency": ["USBR (RISE)", "USGS (NWIS)", "TWDB", "CDEC"],
                    "k": ["USBR", "USGS", "TWDB", "CDEC"]})
rows = []
for _, r in src.iterrows():
    ids = P[P.source == r.k].swot_lake_id.unique()
    pl = P[P.source == r.k]
    ar = pl.groupby("swot_lake_id").pld_area.first()
    iv = inv[inv.swot_lake_id.isin(ids)]
    rows.append([r.Agency, f"{len(ids)}", f"{len(pl):,}", f"{ar.median():.1f}", f"{ar.min():.2f}–{ar.max():.0f}",
                 "stage, storage, area" if r.k == "TWDB" else ("stage, storage, area (" + f"{iv.has_area.sum()}" + ")" if r.k == "USBR" else "stage, storage"),
                 "Jun 2023–Oct 2026" if r.k != "CDEC" else "Jul 2023–Oct 2025"])
ar = P.groupby("swot_lake_id").pld_area.first()
rows.append(["Total", f"{P.swot_lake_id.nunique()}", f"{len(P):,}", f"{ar.median():.1f}", f"{ar.min():.2f}–{ar.max():.0f}", "", ""])
T["table1"] = dict(caption="In-situ records and the evaluation sample (filter F_q0: quality_f = 0, no climatological ice, at least 10 consecutive pairs per lake).",
                   columns=["Agency", "Lakes", "Pairs", "Median area (km²)", "Area range (km²)", "Variables", "Period"],
                   rows=rows, note="Areas are the PLD reference areas. Number in parentheses: lakes with daily surface area. CDEC records were taken from the archived benchmark files [@Levenson2026data].")
names = [("dS_prod", "Product ds1_l (Eq. 1)"), ("dS_prod_q", "Product ds1_q (Eq. 2)"), ("dS_trap", "Pairwise trapezoid (Eq. 4)"),
         ("dS_shyp", "SWOT hypsometry (Eq. 6)"), ("dS_level_only", "Level error only"), ("dS_area_only", "Area error only")]
rows = []
for e, lab in names:
    s = S[e]
    rows.append([lab, f"{s['n_lakes']}", f"{s['nrmse'][0]:.2f} ({s['nrmse'][1]:.2f}–{s['nrmse'][2]:.2f})",
                 f"{1e3 * s['rmse'][0]:.2f}", f"{s['r'][0]:.3f}",
                 "–" if e == "dS_prod" else f"{100 * (W[[ 'dS_prod', e]].dropna()[e] < W[['dS_prod', e]].dropna()['dS_prod']).mean():.0f}"])
T["table2"] = dict(caption="Accuracy of storage change between consecutive valid observations (median over lakes, 95% bootstrap interval in parentheses).",
                   columns=["Estimator", "Lakes", "nRMSE", "RMSE (10⁻³ km³)", "r", "Lakes better than ds1_l (%)"],
                   rows=rows, note="nRMSE is the RMSE of the storage-change error divided by the standard deviation of the in-situ storage change of the same pairs. The last two rows are diagnostic: they combine one SWOT input with in-situ data and are not SWOT-only estimators.")
cv = lambda k: f"{pct(em[k][0])}"
rows = [["Level change (wse_u)", f"{S['cov_dh']['n']:,}", pct(S["cov_dh"]["p1"]), pct(S["cov_dh"]["p196"]), f"{S['cov_dh']['k_nmad']:.1f}"],
        ["Area change (area_tot_u)", f"{S['cov_dA']['n']:,}", pct(S["cov_dA"]["p1"]), pct(S["cov_dA"]["p196"]), f"{S['cov_dA']['k_nmad']:.1f}"],
        ["Storage change (ds1_l_u)", f"{S['cov_dS']['n']:,}", pct(S["cov_dS"]["p1"]), pct(S["cov_dS"]["p196"]), f"{S['cov_dS']['k_nmad']:.1f}"],
        ["Storage change, empirical model (Eq. 7), held-out lakes", "–", cv("prod_p1"), cv("prod_p196"), f"{em['prod_k'][0]:.2f}"],
        ["Trapezoid, empirical model, held-out lakes", "–", cv("trap_p1"), cv("trap_p196"), f"{em['trap_k'][0]:.2f}"]]
T["table3"] = dict(caption="Coverage of stated and empirical uncertainties for changes between consecutive observations.",
                   columns=["Quantity and uncertainty", "Pairs", "Within ±1σ (%)", "Within ±1.96σ (%)", "Inflation factor k"],
                   rows=rows, note=f"Nominal coverage is 68.3% and 95%. k is the normalised median absolute deviation of the standardised errors (1 for a calibrated Gaussian error). Empirical-model rows are medians of 100 two-fold splits by lake (parameters fitted on one half of the lakes, coverage computed on the other); area change uses the {P.dropna(subset=['dA_g']).swot_lake_id.nunique()} lakes with in-situ area in the pair sample.")
rows = []
for c in sens.index:
    rows.append([c.replace("Main (F_q0)", "Main sample (F_q0)").replace("Filter F_q01", "Filter: quality_f ≤ 1").replace("Filter F_all", "No quality filter")
                 .replace("Agency ", "Agency: ").replace("PLD area ", "PLD area: ").replace(" km2", " km²").replace("1e+06", "max").replace("<=", "≤"),
                 f"{sens.loc[c, 'n_lakes']:.0f}", f"{sens.loc[c, 'n_pairs']:,.0f}", f"{sens.loc[c, 'dS_prod']:.2f}",
                 f"{sens.loc[c, 'dS_trap']:.2f}", f"{sens.loc[c, 'dS_shyp']:.2f}", pct(sens.loc[c, "stated_p1"])])
T["tableS1"] = dict(caption="Sensitivity of the main results (median per-lake nRMSE).",
                    columns=["Case", "Lakes", "Pairs", "ds1_l", "Trapezoid", "SWOT hypsometry", "ds1_l_u ±1σ coverage (%)"],
                    rows=rows, note="All cases use the F_q0 filter except the two filter rows. The PGD0/PID0 split follows the change of processing baseline on 8 April 2025.")
json.dump(T, open(OUT / "tables.json", "w"), indent=1, ensure_ascii=False)
print(len(N), "numbers;", list(T))
