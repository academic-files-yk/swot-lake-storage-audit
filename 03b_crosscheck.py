"""Cross-check of our agency downloads against the archived benchmark files (storage ratio)."""
import json
import numpy as np, pandas as pd
from common import *

inv = pd.read_csv(OUT / "insitu_inventory.csv")
rows = []
for _, r in inv[inv.source != "CDEC"].iterrows():
    f = RAW / "benchmark" / f"{r.swot_lake_id}_daily.csv"
    if not f.exists():
        continue
    b = pd.read_csv(f, usecols=["date", "storage"], low_memory=False)
    b["date"] = pd.to_datetime(b.date, format="%m/%d/%y").dt.strftime("%Y-%m-%d")
    b = b.dropna().groupby("date").storage.mean() * 1e-9
    a = pd.read_csv(INSITU / "lake" / f"{r.swot_lake_id}.csv").set_index("date").storage
    j = pd.concat([a, b], axis=1, join="inner").dropna()
    if len(j) >= 20:
        rows.append(dict(swot_lake_id=r.swot_lake_id, source=r.source, n=len(j),
                         ratio=float((j.iloc[:, 0] / j.iloc[:, 1]).median())))
C = pd.DataFrame(rows); C.to_csv(OUT / "insitu_crosscheck.csv", index=False)
out = {"n_lakes": len(C), "frac_within_0p5pct": float(((C.ratio - 1).abs() <= 0.005).mean()),
       "median_ratio": float(C.ratio.median())}
json.dump(out, open(OUT / "insitu_crosscheck.json", "w"), indent=1)
print(out, C.groupby("source").ratio.median().to_dict())
