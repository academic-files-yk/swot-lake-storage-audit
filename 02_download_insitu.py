"""Download in-situ reservoir records (stage, storage, area) for the mapped lakes.
USGS NWIS daily values, USBR RISE, TWDB (Texas). CDEC is not reachable from the
processing network; CDEC records are taken from the archived benchmark files of
Levenson et al. (see 02b)."""
import io, re, time, json, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd, requests
from common import *

D0, D1 = "2023-06-01", "2026-10-15"
S = requests.Session()

def get(url, **kw):
    for k in range(5):
        try:
            r = S.get(url, timeout=180, **kw)
            if r.status_code == 200:
                return r
            if r.status_code in (400, 404):
                return None
        except Exception:
            pass
        time.sleep(5 * (k + 1))
    return None

# ---------------- USGS NWIS daily values ----------------
NWIS_PARAMS = {"00054": "storage_acft", "00062": "elev_ft", "62614": "elev_ft_ngvd29",
               "62615": "elev_ft_navd88", "00065": "gage_ft", "62616": "elev_ft_ngvd29b",
               "62617": "elev_ft_navd88b"}

def nwis(site):
    site = f"{int(site):08d}" if float(site) < 1e8 else str(int(site))
    out = INSITU / "usgs" / f"{site}.csv"
    if out.exists():
        return site, "cached"
    r = get("https://waterservices.usgs.gov/nwis/dv/", params=dict(
        format="rdb", sites=site, parameterCd=",".join(NWIS_PARAMS), startDT=D0,
        endDT=D1, siteStatus="all"))
    if r is None:
        out.write_text(""); return site, "none"
    lines = [l for l in r.text.splitlines() if not l.startswith("#")]
    if len(lines) < 3:
        out.write_text(""); return site, "empty"
    df = pd.read_csv(io.StringIO("\n".join([lines[0]] + lines[2:])), sep="\t", dtype=str)
    keep = {"datetime": "date"}
    for c in df.columns:
        m = re.match(r"\d+_(\d{5})_00003$", c) or re.match(r"\d+_(\d{5})_32400$", c) \
            or re.match(r"\d+_(\d{5})_00011$", c)
        if m and m.group(1) in NWIS_PARAMS and NWIS_PARAMS[m.group(1)] not in keep.values():
            keep[c] = NWIS_PARAMS[m.group(1)]
    df = df[list(keep)].rename(columns=keep)
    for c in df.columns[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df.to_csv(out, index=False)
    return site, ",".join(df.columns[1:])

# ---------------- USBR RISE ----------------
RISE = "https://data.usbr.gov/rise/api"
H = {"accept": "application/vnd.api+json"}
WANT = {"storage": "Lake/Reservoir Storage", "elev": "Lake/Reservoir Elevation",
        "area": "Lake/Reservoir Area"}

def rise_items(loc):
    r = get(f"{RISE}/location/{loc}", headers=H)
    if r is None:
        return {}
    recs = [x["id"] for x in r.json()["data"]["relationships"]["catalogRecords"]["data"]]
    items = {}
    for rec in recs:
        rr = get("https://data.usbr.gov" + rec, headers=H)
        if rr is None:
            continue
        for it in rr.json()["data"]["relationships"]["catalogItems"]["data"]:
            ri = get("https://data.usbr.gov" + it["id"], headers=H)
            if ri is None:
                continue
            a = ri.json()["data"]["attributes"]
            t = a["itemTitle"]
            if "Daily" not in t or a.get("isModeled"):
                continue
            for k, w in WANT.items():
                if w in t and k not in items:
                    items[k] = (a["_id"], t)
    return items

def rise(loc):
    loc = int(loc)
    out = INSITU / "usbr" / f"{loc}.csv"
    if out.exists():
        return loc, "cached"
    items = rise_items(loc)
    series = []
    for k, (iid, title) in items.items():
        rows, page = [], 1
        while True:
            r = get(f"{RISE}/result", headers=H, params={"itemId": iid, "dateTime[after]": D0,
                    "dateTime[before]": D1, "itemsPerPage": 1000, "page": page})
            if r is None:
                break
            j = r.json()
            rows += [(d["attributes"]["dateTime"][:10], d["attributes"]["result"]) for d in j["data"]]
            if "next" not in j.get("links", {}) or not j["data"]:
                break
            page += 1
        if rows:
            s = pd.DataFrame(rows, columns=["date", k]).groupby("date")[k].mean()
            series.append(s)
    meta = {k: v[1] for k, v in items.items()}
    (INSITU / "usbr" / f"{loc}_items.json").write_text(json.dumps(meta, indent=1))
    if not series:
        out.write_text(""); return loc, "empty"
    pd.concat(series, axis=1).reset_index().to_csv(out, index=False)
    return loc, ",".join(items)

# ---------------- TWDB ----------------
def twdb(slug):
    out = INSITU / "twdb" / f"{slug}.csv"
    if out.exists():
        return slug, "cached"
    r = get(f"https://www.waterdatafortexas.org/reservoirs/individual/{slug}.csv")
    if r is None:
        return slug, "none"
    lines = [l for l in r.text.splitlines() if not l.startswith("#")]
    df = pd.read_csv(io.StringIO("\n".join(lines)))
    df = df[df.date >= D0]
    df.to_csv(out, index=False)
    return slug, len(df)

if __name__ == "__main__":
    m = pd.read_csv(MAPPING)
    for d in ("usgs", "usbr", "twdb"):
        (INSITU / d).mkdir(exist_ok=True)
    which = sys.argv[1:] or ["usgs", "usbr"]
    if "usgs" in which:
        sites = m.nwis_gage_id.dropna().unique()
        with ThreadPoolExecutor(4) as ex:
            res = list(ex.map(nwis, sites))
        pd.DataFrame(res, columns=["site", "params"]).to_csv(OUT / "download_log_usgs.csv", index=False)
        print("usgs", len(res))
    if "usbr" in which:
        locs = m.rise_id.dropna().unique()
        with ThreadPoolExecutor(4) as ex:
            res = list(ex.map(rise, locs))
        pd.DataFrame(res, columns=["location", "items"]).to_csv(OUT / "download_log_usbr.csv", index=False)
        print("usbr", len(res))
