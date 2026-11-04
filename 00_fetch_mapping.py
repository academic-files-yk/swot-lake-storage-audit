"""Fetch the gauge-to-PLD mapping table of the Levenson et al. benchmark
(GitHub mirror of Zenodo 10.5281/zenodo.22258932)."""
import requests
from common import *

URL = "https://raw.githubusercontent.com/ericslevenson/swot-s2-lake-storage/main/data/"
for name in ("gage_swot_mapping.csv", "column_guide.csv"):
    r = requests.get(URL + name, timeout=120); r.raise_for_status()
    (RAW / name).write_bytes(r.content)
    print(name, len(r.content), "bytes")
