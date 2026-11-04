"""Shared estimators and statistics."""
import numpy as np, pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.stats import theilslopes
from common import SEED

def filt(d, name):
    d = d[(d.ice_clim_f == 0) & d.wse.notna() & d.area_total.notna() & d.ds1_l.notna()
          & d.g_storage.notna()]
    if name == "F_q0":
        d = d[d.quality_f == 0]
    elif name == "F_q01":
        d = d[d.quality_f <= 1]
    return d.sort_values("time").drop_duplicates("time")

def insitu_hypsometry(ins):
    """Monotone V(h) from daily stage-storage pairs (D8)."""
    s = ins[["stage", "storage"]].dropna()
    if len(s) < 30 or s.stage.nunique() < 10:
        return None
    b = s.groupby((s.stage * 100).round()).storage.median()
    h = b.index.values / 100.0
    v = np.maximum.accumulate(b.values)
    keep = np.r_[True, np.diff(v) > 0]
    h, v = h[keep], v[keep]
    if len(h) < 5:
        return None
    return PchipInterpolator(h, v, extrapolate=True), (h.min(), h.max())

def swot_hypsometry(d):
    """Theil-Sen line A(h) on the lake's valid SWOT pairs; returns V(h1)-V(h0) function (D9b)."""
    h0 = d.wse.median()
    if d.wse.std() < 0.05 or len(d) < 5:
        return None
    b, a, *_ = theilslopes(d.area_total.values, d.wse.values - h0)
    return lambda hi, hj: (a * (hj - hi) + 0.5 * b * ((hj - h0) ** 2 - (hi - h0) ** 2)) / 1e3

def nmad(x):
    x = np.asarray(x); x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) else np.nan

def boot_median(x, n=2000, seed=SEED):
    x = np.asarray(x); x = x[np.isfinite(x)]
    if len(x) < 3:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    bs = np.median(rng.choice(x, (n, len(x))), axis=1)
    return float(np.median(x)), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
