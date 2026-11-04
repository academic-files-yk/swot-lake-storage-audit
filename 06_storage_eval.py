"""Pairwise storage-change evaluation, uncertainty coverage and error decomposition
(D4-D9) for Version D."""
import json, sys
import numpy as np, pandas as pd
from common import *
from lib import *

FILTERS = ["F_q0", "F_q01", "F_all"]
M = pd.read_csv(WORK / "matched_D.csv", parse_dates=["time"])

def pairs_for_lake(d, ins_h, swot_h):
    d = d.reset_index(drop=True)
    i, j = d.index[:-1], d.index[1:]
    P = pd.DataFrame({"swot_lake_id": d.swot_lake_id.iloc[0], "source": d.source.iloc[0],
                      "t_i": d.time.values[i], "t_j": d.time.values[j]})
    g = lambda c: (d[c].values[i], d[c].values[j])
    hi, hj = g("wse"); Ai, Aj = g("area_total"); Si, Sj = g("g_storage")
    gi, gj = g("g_stage"); ui, uj = g("ds1_l_u"); wi, wj = g("wse_u"); ai, aj = g("area_tot_u")
    GAi, GAj = g("g_area")
    P["dt_days"] = (P.t_j - P.t_i).dt.total_seconds() / 86400
    P["dS_ins"] = Sj - Si
    P["dS_prod"] = d.ds1_l.values[j] - d.ds1_l.values[i]
    P["dS_prod_q"] = d.ds1_q.values[j] - d.ds1_q.values[i]
    P["u_prod"] = np.hypot(ui, uj)
    P["dS_trap"] = 0.5 * (Ai + Aj) * (hj - hi) / 1e3
    P["dS_shyp"] = swot_h(hi, hj) if swot_h else np.nan
    P["dh_swot"], P["dh_g"] = hj - hi, gj - gi
    P["u_dh"] = np.hypot(wi, wj)
    P["A_mean"] = 0.5 * (Ai + Aj)
    P["href_i"], P["href_j"] = hi - d.p_ref_wse.values[i], hj - d.p_ref_wse.values[j]
    P["Aref"] = d.p_ref_area.values[i]
    P["dA_swot"], P["dA_g"] = Aj - Ai, GAj - GAi
    P["u_dA"] = np.hypot(ai, aj)
    if ins_h is not None:
        V, (lo, hi_) = ins_h
        P["dS_level_only"] = V(gi + (hj - hi)) - V(gi)
        P["dS_hyp_true"] = V(gj) - V(gi)       # hypsometry closure check
    P["dS_area_only"] = 0.5 * (Ai + Aj) * (gj - gi) / 1e3
    for c in ("dark_frac", "xtrk_dist", "layovr_val", "quality_f", "wse_std"):
        P[c + "_max"] = np.fmax(d[c].values[i], d[c].values[j]) if c != "xtrk_dist" else \
            np.fmin(np.abs(d[c].values[i]), np.abs(d[c].values[j]))
    P["pld_area"] = d.p_ref_area.iloc[0]
    return P

def per_lake(P, est):
    e = P[est] - P.dS_ins
    sd = P.dS_ins.std()
    return dict(rmse=np.sqrt(np.mean(e ** 2)), mae=np.mean(np.abs(e)),
                nrmse=np.sqrt(np.mean(e ** 2)) / sd if sd > 0 else np.nan,
                r=np.corrcoef(P[est], P.dS_ins)[0, 1] if P[est].std() > 0 else np.nan,
                bias=np.median(e))

if __name__ == "__main__":
    summary, lake_rows = {}, []
    for F in FILTERS:
        allP = []
        for lid, d in M.groupby("swot_lake_id"):
            d = filt(d, F)
            if len(d) < 11:
                continue
            ins = pd.read_csv(INSITU / "lake" / f"{lid}.csv")
            allP.append(pairs_for_lake(d, insitu_hypsometry(ins), swot_hypsometry(d)))
        P = pd.concat(allP, ignore_index=True)
        P.to_csv(WORK / f"pairs_{F}.csv", index=False)
        S = {"n_lakes": int(P.swot_lake_id.nunique()), "n_pairs": int(len(P)),
             "dt_median_days": float(P.dt_days.median()),
             "sources": P.groupby("source").swot_lake_id.nunique().to_dict()}
        ests = ["dS_prod", "dS_prod_q", "dS_trap", "dS_shyp", "dS_level_only", "dS_area_only"]
        for est in ests:
            rows = []
            for lid, p in P.groupby("swot_lake_id"):
                p = p.dropna(subset=[est])
                if len(p) >= 10:
                    rows.append(dict(swot_lake_id=lid, est=est, filter=F, n=len(p), **per_lake(p, est)))
            R = pd.DataFrame(rows); lake_rows.append(R)
            S[est] = {k: boot_median(R[k]) for k in ("nrmse", "r", "rmse")} | {"n_lakes": len(R)}
        # coverage of stated uncertainty (D7)
        z = (P.dS_prod - P.dS_ins) / P.u_prod
        zh = ((P.dh_swot - P.dh_g) / P.u_dh).dropna()
        za = ((P.dA_swot - P.dA_g) / P.u_dA).dropna()
        cov = lambda zz: {"n": int(len(zz)), "p1": float((zz.abs() <= 1).mean()),
                          "p196": float((zz.abs() <= 1.96).mean()), "k_nmad": float(nmad(zz)),
                          "med_abs_z": float(zz.abs().median())}
        S["cov_dS"], S["cov_dh"], S["cov_dA"] = cov(z), cov(zh), cov(za)
        # pooled error medians in physical units
        S["dh_err_nmad_m"] = float(nmad(P.dh_swot - P.dh_g))
        S["dh_err_p68_m"] = float((P.dh_swot - P.dh_g).abs().quantile(.68))
        S["u_dh_median_m"] = float(P.u_dh.median())
        summary[F] = S
        print(F, json.dumps({k: S[k] for k in ("n_lakes", "n_pairs")}),
              {e: round(S[e]["nrmse"][0], 3) for e in ests}, S["cov_dS"], S["cov_dh"], S["cov_dA"])
    pd.concat(lake_rows).to_csv(OUT / "lake_metrics.csv", index=False)
    json.dump(summary, open(OUT / "storage_eval.json", "w"), indent=1, default=float)
