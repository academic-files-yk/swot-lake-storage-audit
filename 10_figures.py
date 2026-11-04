"""Figures 1-5 (function names follow the analysis order; file names follow the citation order in the text)."""
import json, sys
import numpy as np, pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from common import *
from lib import *
from figstyle import *

FIG = HERE / "figures"
style()
P = pd.read_csv(WORK / "pairs_F_q0_model.csv")
O = pd.read_csv(WORK / "obs_D_F_q0.csv")
OH = O.dropna(subset=["e_h"]); OH = OH[OH.groupby("swot_lake_id").e_h.transform("size") >= 10]   # as in 07_obs_errors
L = pd.read_csv(OUT / "lake_metrics.csv"); L = L[L["filter"] == "F_q0"]
lakes = P.groupby("swot_lake_id").agg(source=("source", "first"), area=("pld_area", "first"),
                                      n=("dS_ins", "size")).reset_index()
m = pd.read_csv(MAPPING).set_index("swot_lake_id")
lakes["lon"], lakes["lat"] = m.loc[lakes.swot_lake_id, "lon_swot"].values, m.loc[lakes.swot_lake_id, "lat_swot"].values

def fig1():
    fig = plt.figure(figsize=(FULL, 72 * MM))
    ax = fig.add_axes([0.02, 0.02, 0.64, 0.86]); bx = fig.add_axes([0.74, 0.17, 0.24, 0.72])
    st = gpd.read_file(WORK / "map" / "cb_2024_us_state_20m.shp")
    st = st[~st.STUSPS.isin(["AK", "HI", "PR"])].to_crs(5070)
    st.plot(ax=ax, color="#f4f4f2", edgecolor="#b5b5b0", linewidth=0.3)
    g = gpd.GeoDataFrame(lakes, geometry=gpd.points_from_xy(lakes.lon, lakes.lat), crs=4326).to_crs(5070)
    for s, c in AGENCY.items():
        q = g[g.source == s]
        ax.scatter(q.geometry.x, q.geometry.y, s=6 + 7 * np.log10(q.area.clip(0.1) * 10), c=c,
                   edgecolor="white", linewidth=0.3, zorder=3, label=f"{s} ({len(q)})")
    for lon in (-120, -110, -100, -90, -80):
        ln = gpd.GeoSeries(gpd.points_from_xy([lon] * 50, np.linspace(24, 50, 50)), crs=4326).to_crs(5070)
        ax.plot(ln.x, ln.y, color=GRID, lw=0.4, zorder=1)
        p = gpd.GeoSeries(gpd.points_from_xy([lon], [24.5]), crs=4326).to_crs(5070)
        if lon > -120:
            ax.text(p.x[0], p.y[0] - 9e4, f"{-lon}°W", ha="center", va="top", fontsize=6, color=INK2)
    for lat in (30, 40, 50):
        ln = gpd.GeoSeries(gpd.points_from_xy(np.linspace(-125, -66, 80), [lat] * 80), crs=4326).to_crs(5070)
        ax.plot(ln.x, ln.y, color=GRID, lw=0.4, zorder=1)
        p = gpd.GeoSeries(gpd.points_from_xy([-124.0], [lat]), crs=4326).to_crs(5070)
        ax.text(p.x[0] - 5e4, p.y[0], f"{lat}°N", ha="right", va="center", fontsize=6, color=INK2)
    ax.plot([1.2e6, 1.7e6], [0.35e6] * 2, color=INK, lw=1.2)
    ax.text(1.45e6, 0.4e6, "500 km", ha="center", va="bottom", fontsize=6)
    ax.set_xlim(-2.65e6, 2.35e6); ax.set_ylim(0.05e6, 3.25e6); ax.set_axis_off()
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.0), fontsize=6.5, handletextpad=0.2,
              markerscale=1.3, title="Agency (lakes)", title_fontsize=6.5)
    title(ax, "a", f"Gauged reservoirs evaluated (n = {len(lakes)})")
    bins = np.logspace(-1, 3, 17)
    bottom = np.zeros(len(bins) - 1)
    for s, c in AGENCY.items():
        h, _ = np.histogram(lakes[lakes.source == s].area, bins)
        bx.bar(bins[:-1], h, width=np.diff(bins), align="edge", bottom=bottom, color=c,
               edgecolor="white", linewidth=0.4)
        bottom += h
    bx.set_xscale("log"); bx.set_xlabel("PLD reference area (km$^2$)"); bx.set_ylabel("Lakes")
    bx.set_xticks([0.1, 1, 10, 100, 1000]); bx.set_xticklabels(["0.1", "1", "10", "100", "1000"])
    title(bx, "b", "Lake size")
    save(fig, "Fig1", FIG)

def fig2():
    M = pd.read_csv(WORK / "matched_D.csv", parse_dates=["time"])
    A = pd.read_csv(OUT / "anomaly_series.csv")
    A = A[A.n >= 40].sort_values("nrmse")
    picks = [int(A.iloc[int(round(q * (len(A) - 1)))].swot_lake_id) for q in (0.25, 0.5, 0.75)]
    fig, axs = plt.subplots(1, 3, figsize=(FULL, 52 * MM), constrained_layout=True)
    for k, (ax, lid) in enumerate(zip(axs, picks)):
        d = filt(M[M.swot_lake_id == lid], "F_q0")
        ins = pd.read_csv(INSITU / "lake" / f"{lid}.csv", parse_dates=["date"])
        ins = ins[(ins.date >= d.time.min() - pd.Timedelta(days=10)) & (ins.date <= d.time.max() + pd.Timedelta(days=10))]
        ref = d.g_storage.median()
        f = swot_hypsometry(d)
        ax.plot(ins.date, (ins.storage - ref) * 1e3, color=REF, lw=1.0, label="In situ")
        off = (d.ds1_l - d.g_storage).median()
        ax.plot(d.time, (d.ds1_l - off - ref) * 1e3, "o", ms=2.6, mfc=PAL[0], mec="white", mew=0.3, label="LakeSP ds1_l")
        vs = f(d.wse.median(), d.wse); off2 = (vs - d.g_storage).median()
        ax.plot(d.time, (vs - off2 - ref) * 1e3, "^", ms=2.6, mfc="none", mec=PAL[1], mew=0.6, label="SWOT hypsometry")
        a = A.set_index("swot_lake_id").loc[lid]
        name = str(m.loc[lid, "names_swot"]).split(";")[0].title() if pd.notna(m.loc[lid, "names_swot"]) else str(lid)
        title(ax, "abc"[k], f"{name} ({m.loc[lid, 'poly_area_swot']:.0f} km$^2$)")
        ax.set_xlabel(f"nRMSE {a.nrmse:.2f} (product), {a.nrmse_shyp:.2f} (hypsometry)", fontsize=6.5, color=INK2)
        ax.set_ylabel("Storage anomaly (10$^6$ m$^3$)" if k == 0 else "")
        ax.xaxis.set_major_locator(plt.matplotlib.dates.YearLocator())
        ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter("%Y"))
    axs[1].legend(loc="lower left", fontsize=6.3, handlelength=1.2)
    save(fig, "Fig2", FIG)
    return picks

def fig3():
    fig, axs = plt.subplots(1, 2, figsize=(FULL, 68 * MM), constrained_layout=True)
    ax = axs[0]
    zh = ((P.dh_swot - P.dh_g) / P.u_dh).abs().dropna()
    za = ((P.dA_swot - P.dA_g) / P.u_dA).abs().dropna()
    zs = ((P.dS_prod - P.dS_ins) / P.u_prod).abs()
    zm = ((P.dS_prod - P.dS_ins) / P.sig_prod).abs()
    from scipy.stats import norm
    p = np.linspace(0.01, 0.99, 99); q = norm.ppf(0.5 + p / 2)
    for z, lab, c, ls in ((zh, "Level change, stated wse_u", PAL[2], "-"),
                          (za, "Area change, stated area_tot_u", PAL[3], "-"),
                          (zs, "Storage change, stated ds1_l_u", PAL[0], "-"),
                          (zm, "Storage change, empirical model", PAL[0], "--")):
        ax.plot(p, [(z <= qq).mean() for qq in q], color=c, ls=ls, lw=1.3, label=lab)
    ax.plot([0, 1], [0, 1], color=REF, lw=0.8, ls=":")
    ax.set_xlabel("Nominal coverage"); ax.set_ylabel("Observed coverage")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=6.3, handlelength=1.8)
    title(ax, "a", "Coverage of stated intervals (pairs)")
    ax = axs[1]
    h = OH
    hb = ax.hexbin(h.wse_u.clip(1e-3), h.e_h.abs().clip(1e-3), xscale="log", yscale="log", gridsize=36,
                   bins="log", cmap="Blues", mincnt=1, linewidths=0.1, extent=(-3, 0.5, -3, 1))
    ax.plot([1e-3, 3], [1e-3, 3], color=INK2, lw=0.8, ls=":")
    ax.axhline(nmad(h.e_h), color=PAL[1], lw=1.0, ls="--")
    ax.text(1.2e-3, nmad(h.e_h) * 1.25, f"NMAD of error {nmad(h.e_h) * 100:.1f} cm", color=INK, fontsize=6.5,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.85))
    ax.set_xlabel("Stated uncertainty wse_u (m)"); ax.set_ylabel("|WSE anomaly error| (m)")
    cb = fig.colorbar(hb, ax=ax, shrink=0.8, pad=0.02); cb.set_label("Observations", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    title(ax, "b", f"WSE error vs stated uncertainty (n = {len(h):,})")
    save(fig, "Fig5", FIG)

def fig4():
    fig, axs = plt.subplots(1, 2, figsize=(FULL, 66 * MM), constrained_layout=True,
                            gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axs[0]
    names = [("dS_prod", "Product\nds1_l"), ("dS_trap", "Pairwise\ntrapezoid"),
             ("dS_shyp", "SWOT\nhypsometry"), ("dS_level_only", "Level error\nonly"),
             ("dS_area_only", "Area error\nonly")]
    data = [L[L.est == e].nrmse.dropna().values for e, _ in names]
    cols = [PAL[0], PAL[1], PAL[4], REF, REF]
    bp = ax.boxplot(data, widths=0.55, showfliers=False, patch_artist=True, medianprops=dict(color=INK, lw=1))
    for b, c in zip(bp["boxes"], cols):
        b.set(facecolor=c, alpha=0.35, edgecolor=c, lw=0.8)
    for w in bp["whiskers"] + bp["caps"]:
        w.set(color=INK2, lw=0.6)
    for i, x in enumerate(data):
        ax.text(i + 1, np.median(x), f"{np.median(x):.2f}", ha="left", va="bottom", fontsize=6.3,
                color=INK, transform=ax.transData)
    ax.set_xticks(range(1, 6)); ax.set_xticklabels([n for _, n in names], fontsize=6.5)
    ax.axvline(3.5, color=GRID, lw=0.8)
    ax.set_ylabel("nRMSE of storage change"); ax.set_ylim(0, 1.0)
    title(ax, "a", f"Per-lake error (n = {len(data[0])} lakes)")
    ax = axs[1]
    xm = np.maximum(P.href_i.abs(), P.href_j.abs())
    P["xb"] = pd.qcut(xm, 8)
    rng = np.random.default_rng(SEED); xs, ys, lo, hi = [], [], [], []
    for k, g in P.groupby("xb", observed=True):
        r = (g.e_prod.abs() / g.e_trap.abs().clip(1e-9)).values if "e_prod" in g else \
            ((g.dS_prod - g.dS_ins).abs() / (g.dS_trap - g.dS_ins).abs().clip(1e-9)).values
        bs = [np.median(rng.choice(r, len(r))) for _ in range(500)]
        xs.append(xm[g.index].median()); ys.append(np.median(r)); lo.append(np.percentile(bs, 2.5)); hi.append(np.percentile(bs, 97.5))
    ax.errorbar(xs, ys, yerr=[np.array(ys) - lo, np.array(hi) - ys], fmt="o", ms=4, color=PAL[0],
                ecolor=PAL[0], elinewidth=0.8, capsize=1.5)
    ax.axhline(1, color=REF, lw=0.8, ls=":")
    ax.set_xscale("log"); ax.set_xlabel("Distance of the pair from the PLD reference level (m)")
    ax.set_ylabel("Median |error| ratio, product / trapezoid")
    ax.set_xticks([0.1, 0.3, 1, 3, 10]); ax.set_xticklabels(["0.1", "0.3", "1", "3", "10"])
    title(ax, "b", "Error growth away from the reference")
    save(fig, "Fig3", FIG)

def fig5():
    fig, axs = plt.subplots(1, 3, figsize=(FULL, 56 * MM), constrained_layout=True)
    ax = axs[0]
    h = OH
    b = pd.cut(h.p_ref_area, [0, 1, 3, 10, 30, 100, 1e4])
    g = h.groupby(b, observed=True)
    x = np.arange(len(g)); ax.bar(x - 0.18, [nmad(v.e_h) * 100 for _, v in g], 0.36, color=PAL[0], label="Observed (NMAD)")
    ax.bar(x + 0.18, [v.wse_u.median() * 100 for _, v in g], 0.36, color=PAL[3], label="Stated (median wse_u)")
    ax.set_xticks(x); ax.set_xticklabels(["<1", "1–3", "3–\n10", "10–\n30", "30–\n100", ">100"], fontsize=6.5)
    ax.set_xlabel("PLD reference area (km$^2$)"); ax.set_ylabel("WSE error (cm)")
    ax.set_ylim(0, 10.5); ax.legend(fontsize=6.3, loc="upper right", ncol=1)
    title(ax, "a", "Level error by lake size")
    ax = axs[1]
    a = O.dropna(subset=["e_a_rel"])
    a = a[a.groupby("swot_lake_id").e_a_rel.transform("size") >= 10]
    b = pd.qcut(a.dark_frac, 6)
    g = a.groupby(b, observed=True)
    xm_ = [v.dark_frac.median() * 100 for _, v in g]
    ax.plot(xm_, [v.e_a_rel.median() * 100 for _, v in g], "o-", color=PAL[1], ms=3.5, lw=1, label="Median (bias)")
    ax.fill_between(xm_, [v.e_a_rel.quantile(.25) * 100 for _, v in g], [v.e_a_rel.quantile(.75) * 100 for _, v in g],
                    color=PAL[1], alpha=0.18, lw=0, label="Interquartile range")
    ax.axhline(0, color=REF, lw=0.8, ls=":")
    ax.set_xscale("log"); ax.set_xticks([0.3, 1, 3, 10]); ax.set_xticklabels(["0.3", "1", "3", "10"])
    ax.set_xlabel("Dark-water fraction (%)"); ax.set_ylabel("Area error (%)")
    ax.legend(fontsize=6.3, loc="upper left")
    title(ax, "b", f"Area error, {a.swot_lake_id.nunique()} lakes")
    ax = axs[2]
    bins = [0, 0.05, 0.10, 0.20, 0.50, 100]; labs = ["<5", "5–10", "10–20", "20–50", ">50"]
    adh = P.dh_g.abs(); ok = P.dh_g.notna()
    for est, c, mk, lab in (("dS_prod", PAL[0], "o", "Product ds1_l"), ("dS_trap", PAL[1], "^", "Pairwise trapezoid")):
        v = [(np.sign(P.loc[ok & (adh >= lo) & (adh < hi), est]) == np.sign(P.loc[ok & (adh >= lo) & (adh < hi), "dS_ins"])).mean() * 100
             for lo, hi in zip(bins[:-1], bins[1:])]
        ax.plot(range(5), v, marker=mk, color=c, ms=3.5, lw=1, label=lab)
    ax.axhline(50, color=REF, lw=0.8, ls=":")
    ax.set_xticks(range(5)); ax.set_xticklabels(labs, fontsize=6.5); ax.set_ylim(45, 102)
    ax.set_xlabel("In-situ level change (cm)"); ax.set_ylabel("Correct sign of storage change (%)")
    ax.legend(fontsize=6.3, loc="center right", bbox_to_anchor=(1.0, 0.32))
    title(ax, "c", "Detection of fill or drawdown")
    save(fig, "Fig4", FIG)

if __name__ == "__main__":
    which = sys.argv[1:] or ["1", "2", "3", "4", "5"]
    for w in which:
        r = globals()[f"fig{w}"]()
        if w == "2":
            json.dump([int(x) for x in r], open(OUT / "fig2_lakes.json", "w"))
