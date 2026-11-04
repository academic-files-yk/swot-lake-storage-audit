"""TWDB daily level/area/storage for mapped Texas (and border) reservoirs.
The PLD-to-TWDB correspondence was set by reservoir name and checked against the
PLD centroid (table written to output/twdb_match.csv)."""
import pandas as pd
from common import *
import importlib
dl = importlib.import_module("02_download_insitu")

SLUG = {7510028163: "choke-canyon", 7510029063: "medina", 7510071263: "arlington",
        7510074813: "worth", 7510091793: "bridgeport", 7510097803: "o-h-ivie",
        7510099483: "j-b-thomas", 7510100243: "e-v-spence", 7510100453: "nasworthy",
        7510103363: "tawakoni", 7510103403: "cedar-creek", 7510103423: "richland-chambers",
        7510104533: "eagle-mountain", 7510104563: "corpus-christi", 7510105803: "houston",
        7510124273: "conroe", 7510138853: "livingston", 7510140043: "fork",
        7510150343: "nacogdoches", 7510154643: "toledo-bend", 7420196883: "greenbelt",
        7520024523: "elephant-butte", 7420024963: "caddo"}

if __name__ == "__main__":
    m = pd.read_csv(MAPPING).set_index("swot_lake_id")
    rows = []
    for lid, slug in SLUG.items():
        s, n = dl.twdb(slug)
        rows.append(dict(swot_lake_id=lid, slug=slug, n=n, gage_name=m.loc[lid, "gage_name"],
                         names_swot=m.loc[lid, "names_swot"], lat=m.loc[lid, "lat_swot"],
                         lon=m.loc[lid, "lon_swot"], pld_area=m.loc[lid, "poly_area_swot"]))
    pd.DataFrame(rows).to_csv(OUT / "twdb_match.csv", index=False)
    print(pd.DataFrame(rows)[["slug", "n", "names_swot", "pld_area"]].to_string())
