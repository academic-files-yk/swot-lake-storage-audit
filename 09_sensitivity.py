"""Sensitivity of the main results to the agency, lake size, pair spacing and filter."""
import json
import numpy as np, pandas as pd
from common import *
from lib import *

def summarise(P, label):
    rows = {}
    for est in ("dS_prod", "dS_trap", "dS_shyp"):
        v = []
        for lid, p in P.groupby("swot_lake_id"):
            p = p.dropna(subset=[est])
            if len(p) >= 10 and p.dS_ins.std() > 0:
                v.append(np.sqrt(np.mean((p[est] - p.dS_ins) ** 2)) / p.dS_ins.std())
        rows[est] = np.median(v) if v else np.nan
        rows["n_lakes"] = len(v)
    z = (P.dS_prod - P.dS_ins) / P.u_prod
    rows["stated_p1"] = float((z.abs() <= 1).mean())
    rows["n_pairs"] = len(P)
    return dict(case=label, **rows)

if __name__ == "__main__":
    res = []
    P = pd.read_csv(WORK / "pairs_F_q0.csv", parse_dates=["t_i", "t_j"])
    res.append(summarise(P, "Main (F_q0)"))
    for F in ("F_q01", "F_all"):
        res.append(summarise(pd.read_csv(WORK / f"pairs_{F}.csv"), f"Filter {F}"))
    for s in ("USBR", "CDEC", "USGS", "TWDB"):
        res.append(summarise(P[P.source == s], f"Agency {s}"))
    res.append(summarise(P[P.source != "CDEC"], "Without CDEC archive"))
    for lo, hi in ((0, 10), (10, 100), (100, 1e6)):
        res.append(summarise(P[(P.pld_area >= lo) & (P.pld_area < hi)], f"PLD area {lo}-{hi:g} km2"))
    res.append(summarise(P[P.dt_days <= 11], "Pairs <= 11 days apart"))
    res.append(summarise(P[P.dt_days > 11], "Pairs > 11 days apart"))
    res.append(summarise(P[P.t_j < "2025-04-08"], "Before 8 Apr 2025 (PGD0 era)"))
    res.append(summarise(P[P.t_i >= "2025-04-08"], "From 8 Apr 2025"))
    R = pd.DataFrame(res); R.to_csv(OUT / "sensitivity.csv", index=False)
    print(R.round(3).to_string())
