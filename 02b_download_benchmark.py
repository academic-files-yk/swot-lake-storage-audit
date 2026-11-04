"""Fetch the archived benchmark time series of Levenson et al. (GitHub mirror of
Zenodo 10.5281/zenodo.22258932). Used only for (i) CDEC gauge records, which could
not be retrieved directly, and (ii) a cross-check of our own gauge downloads."""
from concurrent.futures import ThreadPoolExecutor
import requests, pandas as pd
from common import *

BASE = "https://raw.githubusercontent.com/ericslevenson/swot-s2-lake-storage/main/data/benchmark_timeseries/"
TREE = "https://api.github.com/repos/ericslevenson/swot-s2-lake-storage/git/trees/main?recursive=1"
DST = RAW / "benchmark"; DST.mkdir(exist_ok=True)

def fetch(name):
    out = DST / name
    if out.exists():
        return
    r = requests.get(BASE + name, timeout=120)
    if r.status_code == 200:
        out.write_bytes(r.content)

if __name__ == "__main__":
    names = [x["path"].split("/")[-1] for x in requests.get(TREE, timeout=60).json()["tree"]
             if x["path"].startswith("data/benchmark_timeseries/")]
    with ThreadPoolExecutor(6) as ex:
        list(ex.map(fetch, names))
    print(len(list(DST.glob("*.csv"))), "files")
