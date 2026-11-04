"""Shared paths and constants for the LakeSP storage-change audit."""
from pathlib import Path

WORK = Path.home() / "swot_p37_work"          # large data (outside OneDrive)
RAW = WORK / "raw"
SWOT = WORK / "swot"
INSITU = WORK / "insitu"
HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
for p in (RAW, SWOT, INSITU, OUT):
    p.mkdir(parents=True, exist_ok=True)

MAPPING = RAW / "gage_swot_mapping.csv"      # Levenson et al. gauge-PLD table
FILL = -1e11                                  # values below this are fill
START, END = "2023-07-01T00:00:00Z", "2026-10-15T00:00:00Z"
ACRE_FT_KM3 = 1233.48183754752e-9             # 1 acre-foot in km3
FT_M = 0.3048
SEED = 20261018
