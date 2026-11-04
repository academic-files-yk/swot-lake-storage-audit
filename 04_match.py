"""Match every SWOT LakeSP observation (Version D and C) to the in-situ record (D3)."""
import numpy as np, pandas as pd
from common import *

NUM = ["wse", "wse_u", "wse_r_u", "wse_std", "area_total", "area_tot_u", "area_detct",
       "ds1_l", "ds1_l_u", "ds1_q", "ds1_q_u", "ds2_l", "ds2_l_u", "ds2_q", "ds2_q_u",
       "dark_frac", "xtrk_dist", "layovr_val", "xovr_cal_c", "p_ref_wse", "p_ref_area",
       "p_ds_t0", "p_storage"]

def load_swot(ver, lid):
    f = SWOT / ver / f"{lid}.csv"
    if not f.exists() or f.stat().st_size == 0:
        return None
    d = pd.read_csv(f, dtype={"lake_id": str, "lake_name": str})
    d = d[d.time_str != "no_data"].copy()
    for c in NUM:
        if c in d:
            d[c] = pd.to_numeric(d[c], errors="coerce")
            d.loc[d[c] < FILL, c] = np.nan
    d["time"] = pd.to_datetime(d.time_str, utc=True).dt.tz_localize(None)
    return d

def interp(ins, t, lon, col):
    """Linear interpolation of daily values placed at local solar noon (D3)."""
    s = ins[["date", col]].dropna()
    if len(s) < 2:
        return np.full(len(t), np.nan)
    tt = pd.to_datetime(s.date) + pd.Timedelta(hours=12) - pd.Timedelta(hours=lon / 15)
    x = tt.values.astype("datetime64[s]").astype(float)
    xq = t.values.astype("datetime64[s]").astype(float)
    i = np.searchsorted(x, xq)
    ok = (i > 0) & (i < len(x))
    i = np.clip(i, 1, len(x) - 1)
    gap = (np.minimum(xq - x[i - 1], 1e12) > 1.5 * 86400) | ((x[i] - xq) > 1.5 * 86400)
    y = np.interp(xq, x, s[col].values)
    y[~ok | gap] = np.nan
    return y

if __name__ == "__main__":
    inv = pd.read_csv(OUT / "insitu_inventory.csv")
    m = pd.read_csv(MAPPING).set_index("swot_lake_id")
    for ver in ("D", "C"):
        parts = []
        for _, r in inv.iterrows():
            d = load_swot(ver, r.swot_lake_id)
            if d is None or d.empty:
                continue
            ins = pd.read_csv(INSITU / "lake" / f"{r.swot_lake_id}.csv")
            lon = m.loc[r.swot_lake_id, "lon_swot"]
            for c in ("stage", "storage", "area"):
                d["g_" + c] = interp(ins, d.time, lon, c)
            d["source"] = r.source
            d["swot_lake_id"] = r.swot_lake_id
            parts.append(d)
        allm = pd.concat(parts, ignore_index=True)
        allm.to_csv(WORK / f"matched_{ver}.csv", index=False)
        print(ver, len(allm), allm.swot_lake_id.nunique(), allm.g_storage.notna().sum())
