"""Build one daily in-situ table per lake: stage (m, local datum), storage (km3),
measured area (km2, where the agency reports it)."""
import json
import numpy as np, pandas as pd
from common import *

ACRE_KM2 = 0.00404685642
m = pd.read_csv(MAPPING)
twdb = pd.read_csv(OUT / "twdb_match.csv").set_index("swot_lake_id").slug.to_dict()
DST = INSITU / "lake"; DST.mkdir(exist_ok=True)

def from_twdb(slug):
    d = pd.read_csv(INSITU / "twdb" / f"{slug}.csv")
    return pd.DataFrame({"date": d.date, "stage": d.water_level * FT_M,
                         "storage": d.reservoir_storage * ACRE_FT_KM3,
                         "area": d.surface_area * ACRE_KM2})

def from_usbr(loc):
    f = INSITU / "usbr" / f"{int(loc)}.csv"
    if not f.exists() or f.stat().st_size == 0:
        return None
    d = pd.read_csv(f)
    out = pd.DataFrame({"date": d.date})
    out["stage"] = d.elev * FT_M if "elev" in d else np.nan
    out["storage"] = d.storage * ACRE_FT_KM3 if "storage" in d else np.nan
    out["area"] = d.area * ACRE_KM2 if "area" in d else np.nan
    return out

def from_usgs(site):
    f = INSITU / "usgs" / f"{int(site):08d}.csv"
    if not f.exists() or f.stat().st_size == 0:
        return None
    d = pd.read_csv(f)
    if "storage_acft" not in d:
        return None
    st = np.nan
    for c in ("elev_ft_navd88", "elev_ft_ngvd29", "elev_ft", "elev_ft_navd88b",
              "elev_ft_ngvd29b", "gage_ft"):
        if c in d and d[c].notna().sum() > 30:
            st = d[c] * FT_M; break
    return pd.DataFrame({"date": d.date, "stage": st, "storage": d.storage_acft * ACRE_FT_KM3,
                         "area": np.nan})

def from_benchmark(lid):
    f = RAW / "benchmark" / f"{lid}_daily.csv"
    if not f.exists():
        return None
    d = pd.read_csv(f, usecols=["date", "stage", "storage"], low_memory=False)
    d["date"] = pd.to_datetime(d.date, format="%m/%d/%y").dt.strftime("%Y-%m-%d")
    d = d.dropna(subset=["storage"]).groupby("date", as_index=False).mean()
    # storage is reported in m3 in the archived files
    return pd.DataFrame({"date": d.date, "stage": d.stage, "storage": d.storage * 1e-9,
                         "area": np.nan})

def despike(d, col, k=10.0):
    """Remove isolated daily values that depart from a centred 7-day running median by more
    than k NMAD of the day-to-day differences (D3b)."""
    x = d[col]
    if x.notna().sum() < 30:
        return d, 0
    med = x.rolling(7, center=True, min_periods=3).median()
    dx = x.diff().dropna()
    s = 1.4826 * np.median(np.abs(dx - dx.median()))
    s = max(s, 1e-9 if col == "storage" else 0.003)
    bad = (x - med).abs() > k * s
    d.loc[bad, col] = np.nan
    return d, int(bad.sum())

rows = []
for _, r in m.iterrows():
    lid = int(r.swot_lake_id)
    cands = []
    if lid in twdb:
        cands.append(("TWDB", lambda: from_twdb(twdb[lid])))
    if pd.notna(r.rise_id):
        cands.append(("USBR", lambda: from_usbr(r.rise_id)))
    if pd.notna(r.nwis_gage_id):
        cands.append(("USGS", lambda: from_usgs(r.nwis_gage_id)))
    if r.agency == "CDEC":
        cands.append(("CDEC", lambda: from_benchmark(lid)))
    for src, fn in cands:
        d = fn()
        if d is None:
            continue
        d = d[d.storage > 0].dropna(subset=["storage"])
        if len(d) < 30:
            continue
        d = d.groupby("date", as_index=False).mean().sort_values("date")
        d = d.set_index(pd.to_datetime(d.date)).asfreq("D")
        d["date"] = d.index.strftime("%Y-%m-%d")
        nsp = 0
        for c in ("storage", "stage", "area"):
            d, k_ = despike(d, c); nsp += k_
        d = d.dropna(subset=["storage"]).reset_index(drop=True)
        d.to_csv(DST / f"{lid}.csv", index=False)
        rows.append(dict(swot_lake_id=lid, source=src, n_days=len(d),
                         has_stage=d.stage.notna().sum() > 30, has_area=d.area.notna().sum() > 30,
                         d0=d.date.iloc[0], d1=d.date.iloc[-1],
                         s_min=d.storage.min(), s_max=d.storage.max(), n_despiked=nsp))
        break
inv = pd.DataFrame(rows)
inv.to_csv(OUT / "insitu_inventory.csv", index=False)
print(inv.source.value_counts(), inv[["has_stage", "has_area"]].sum(), len(inv), inv.n_despiked.sum(), inv.n_days.sum())
