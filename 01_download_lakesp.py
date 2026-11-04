"""Download SWOT LakeSP Prior time series (Version D and Version C) from Hydrocron
for every PLD lake in the gauge mapping table."""
import io, json, time, sys
from concurrent.futures import ThreadPoolExecutor
import pandas as pd, requests
from common import *

URL = "https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries"
FIELDS = ("lake_id,time_str,time_tai,cycle_id,pass_id,crid,wse,wse_u,wse_r_u,wse_std,"
          "area_total,area_tot_u,area_detct,area_det_u,layovr_val,xtrk_dist,"
          "ds1_l,ds1_l_u,ds1_q,ds1_q_u,ds2_l,ds2_l_u,ds2_q,ds2_q_u,"
          "quality_f,dark_frac,ice_clim_f,ice_dyn_f,partial_f,xovr_cal_q,xovr_cal_c,"
          "geoid_hght,p_ref_wse,p_ref_area,p_date_t0,p_ds_t0,p_storage,p_lon,p_lat,lake_name")
COLL = {"D": "SWOT_L2_HR_LakeSP_D", "C": "SWOT_L2_HR_LakeSP_2.0"}

def fetch(lake, ver):
    out = SWOT / ver / f"{lake}.csv"
    if out.exists():
        return lake, "cached"
    params = dict(feature="PriorLake", feature_id=str(lake), start_time=START,
                  end_time=END, output="csv", fields=FIELDS, collection_name=COLL[ver])
    for k in range(5):
        try:
            r = requests.get(URL, params=params, timeout=300)
            if r.status_code == 200:
                df = pd.read_csv(io.StringIO(r.json()["results"]["csv"]), dtype={"lake_id": str})
                df.to_csv(out, index=False)
                return lake, len(df)
            if r.status_code in (400, 404):
                out.write_text("")
                return lake, f"http{r.status_code}"
        except Exception as e:
            err = e
        time.sleep(5 * (k + 1))
    return lake, "failed"

if __name__ == "__main__":
    ids = pd.read_csv(MAPPING).swot_lake_id.astype(str).unique()
    for ver in sys.argv[1:] or ["D", "C"]:
        (SWOT / ver).mkdir(exist_ok=True)
        with ThreadPoolExecutor(4) as ex:
            res = list(ex.map(lambda i: fetch(i, ver), ids))
        log = pd.DataFrame(res, columns=["lake_id", "status"])
        log.to_csv(OUT / f"download_log_lakesp_{ver}.csv", index=False)
        print(ver, log.status.astype(str).value_counts().head())
