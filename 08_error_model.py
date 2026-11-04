"""Empirical first-order error model for the product and the pairwise trapezoid,
fitted on half of the lakes and tested on the other half (D10); detectability of
storage change and the anomaly-series view (D5)."""
import json
import numpy as np, pandas as pd
from common import *
from lib import *

P = pd.read_csv(WORK / "pairs_F_q0.csv")
O = pd.read_csv(WORK / "obs_D_F_q0.csv")

def sig_prod(P, sh, r):
    Ai, Aj = P.A_mean - P.dA_swot / 2, P.A_mean + P.dA_swot / 2
    lev = 0.25 * ((Aj + P.Aref) ** 2 + (Ai + P.Aref) ** 2) * sh ** 2
    are = 0.25 * ((P.href_j * Aj) ** 2 + (P.href_i * Ai) ** 2) * r ** 2
    return np.sqrt(lev + are) / 1e3

def sig_trap(P, sh, r):
    Ai, Aj = P.A_mean - P.dA_swot / 2, P.A_mean + P.dA_swot / 2
    lev = 0.25 * (Ai + Aj) ** 2 * 2 * sh ** 2
    are = 0.25 * P.dh_swot ** 2 * (Ai ** 2 + Aj ** 2) * r ** 2
    return np.sqrt(lev + are) / 1e3

def cov(z):
    z = z[np.isfinite(z)]
    return (z.abs() <= 1).mean(), (z.abs() <= 1.96).mean(), nmad(z)

if __name__ == "__main__":
    rng = np.random.default_rng(SEED)
    lakes = P.swot_lake_id.unique()
    alake = O.dropna(subset=["e_a_anom"]).swot_lake_id.unique()
    rows = []
    for k in range(100):
        tr = set(rng.choice(lakes, len(lakes) // 2, replace=False))
        tra = tr & set(alake)          # area lakes inside the same training half
        sh = nmad(O[O.swot_lake_id.isin(tr)].e_h)
        r = nmad(O[O.swot_lake_id.isin(tra)].e_a_anom)
        T = P[~P.swot_lake_id.isin(tr)]
        zp = (T.dS_prod - T.dS_ins) / sig_prod(T, sh, r)
        zt = (T.dS_trap - T.dS_ins) / sig_trap(T, sh, r)
        zs = (T.dS_prod - T.dS_ins) / T.u_prod
        rows.append(dict(rep=k, sigma_h=sh, r=r, *[], **dict(zip(["prod_p1", "prod_p196", "prod_k"], cov(zp))),
                         **dict(zip(["trap_p1", "trap_p196", "trap_k"], cov(zt))),
                         **dict(zip(["stated_p1", "stated_p196", "stated_k"], cov(zs)))))
    R = pd.DataFrame(rows); R.to_csv(OUT / "error_model_cv.csv", index=False)
    out = {c: [float(R[c].median()), float(R[c].quantile(.05)), float(R[c].quantile(.95))] for c in R.columns if c != "rep"}
    # full-sample parameters for reporting and for the figures
    sh, r = nmad(O.e_h), nmad(O.dropna(subset=["e_a_anom"]).e_a_anom)
    out["sigma_h_all"], out["r_all"] = float(sh), float(r)
    P["sig_prod"], P["sig_trap"] = sig_prod(P, sh, r), sig_trap(P, sh, r)
    P["sig_prod_lev"] = sig_prod(P, sh, 0); P["sig_prod_area"] = sig_prod(P, 0, r)
    out["share_area_var_prod_median"] = float((P.sig_prod_area ** 2 / P.sig_prod ** 2).median())
    out["ratio_sigprod_sigtrap_median"] = float((P.sig_prod / P.sig_trap).median())
    out["ratio_emp_to_stated_median"] = float((P.sig_prod / P.u_prod).median())
    # detectability: sign agreement by true level change
    bins = [0, 0.05, 0.10, 0.20, 0.50, 100]
    P["adh"] = P.dh_g.abs()
    det = {}
    for lo, hi in zip(bins[:-1], bins[1:]):
        s = P[(P.adh >= lo) & (P.adh < hi)].dropna(subset=["dh_g"])
        det[f"{lo}-{hi}"] = {"n": int(len(s)),
                             "sign_prod": float((np.sign(s.dS_prod) == np.sign(s.dS_ins)).mean()),
                             "sign_trap": float((np.sign(s.dS_trap) == np.sign(s.dS_ins)).mean())}
    out["detect_by_dh"] = det
    P.to_csv(WORK / "pairs_F_q0_model.csv", index=False)
    # anomaly-series view (D5)
    M = pd.read_csv(WORK / "matched_D.csv", parse_dates=["time"])
    rows = []
    for lid, d in M.groupby("swot_lake_id"):
        d = filt(d, "F_q0")
        if len(d) < 11:
            continue
        e = d.ds1_l - d.g_storage; e = e - e.median()
        sd = d.g_storage.std()
        f = swot_hypsometry(d)
        vs = f(d.wse.median(), d.wse) if f else np.nan * d.wse
        es = vs - d.g_storage; es = es - es.median()
        x = (d.wse - d.p_ref_wse).abs()
        rows.append(dict(swot_lake_id=lid, n=len(d), nrmse=np.sqrt((e ** 2).mean()) / sd,
                         r=np.corrcoef(d.ds1_l, d.g_storage)[0, 1],
                         nrmse_shyp=np.sqrt((es ** 2).mean()) / sd, x_med=x.median(),
                         wse_range=d.wse.quantile(.95) - d.wse.quantile(.05)))
    A = pd.DataFrame(rows); A.to_csv(OUT / "anomaly_series.csv", index=False)
    out["anomaly_nrmse"] = boot_median(A.nrmse); out["anomaly_r"] = boot_median(A.r)
    out["anomaly_n_lakes"] = len(A)
    out["anomaly_nrmse_shyp"] = boot_median(A.nrmse_shyp)
    out["anomaly_frac_shyp_better"] = float((A.nrmse_shyp < A.nrmse).mean())
    out["x_ref_median_m"] = boot_median(A.x_med); out["wse_range_median_m"] = boot_median(A.wse_range)
    xm = np.maximum(P.href_i.abs(), P.href_j.abs())
    out["pair_x_over_dh_median"] = float((xm / P.dh_swot.abs().clip(0.01)).median())
    q = pd.qcut(xm, 4)
    P["e_prod"], P["e_trap"] = P.dS_prod - P.dS_ins, P.dS_trap - P.dS_ins
    out["ratio_by_x_quartile"] = {str(k): {"n": int(len(g)), "x_med": float(xm[g.index].median()),
        "ratio_abs_err": float((g.e_prod.abs() / g.e_trap.abs().clip(1e-9)).median()),
        "frac_prod_worse": float((g.e_prod.abs() > g.e_trap.abs()).mean())}
        for k, g in P.groupby(q, observed=True)}
    # empirical 95 % multiplier of the fitted model
    zp = (P.dS_prod - P.dS_ins) / P.sig_prod
    out["model_abs_z_q68"], out["model_abs_z_q95"] = float(zp.abs().quantile(.683)), float(zp.abs().quantile(.95))
    json.dump(out, open(OUT / "error_model.json", "w"), indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))
