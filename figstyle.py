"""Journal figure style, text-overlap check and PDF + 600 dpi PNG/TIFF export."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MM = 1 / 25.4
FULL, COL = 165 * MM, 90 * MM
INK, INK2, GRID, REF = "#1a1a1a", "#555555", "#e3e3e3", "#7a7974"
PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#6b4fbb"]
AGENCY = {"USBR": PAL[0], "CDEC": PAL[1], "USGS": PAL[2], "TWDB": PAL[3]}

def style():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "legend.fontsize": 7, "axes.edgecolor": INK2, "axes.labelcolor": INK,
        "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False, "grid.color": GRID,
        "grid.linewidth": 0.5, "legend.frameon": False, "savefig.facecolor": "white",
        "pdf.fonttype": 42, "ps.fonttype": 42, "mathtext.default": "regular"})

def title(ax, letter, text):
    ax.set_title(f"({letter}) {text}", loc="left", fontsize=8, color=INK, pad=4)

def check(fig, name):
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    hidden = set()
    for ax in fig.axes:
        if not ax.axison:
            hidden.update(id(t) for t in ax.get_xticklabels() + ax.get_yticklabels())
            hidden.update({id(ax.xaxis.label), id(ax.yaxis.label), id(ax.xaxis.offsetText), id(ax.yaxis.offsetText)})
        for axis, lim in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            lo, hi = min(lim), max(lim)
            for t in axis.get_major_ticks():
                if not lo - 1e-9 <= t.get_loc() <= hi + 1e-9:
                    hidden.update({id(t.label1), id(t.label2)})
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()
             and id(t) not in hidden and t.get_window_extent(r).width > 0]
    boxes = [(t, t.get_window_extent(r).padded(-0.5)) for t in texts]
    issues = [f"outside canvas: {t.get_text()[:30]!r}" for t, b in boxes
              if b.x0 < -1 or b.y0 < -1 or b.x1 > W + 1 or b.y1 > H + 1]
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            if boxes[i][1].overlaps(boxes[j][1]):
                issues.append(f"overlap: {boxes[i][0].get_text()[:25]!r} <-> {boxes[j][0].get_text()[:25]!r}")
    print(f"{name}: {len(texts)} text items, {'no issues' if not issues else str(len(issues)) + ' issues'}")
    for s in issues:
        print("   ", s)
    return issues

def save(fig, name, folder):
    from PIL import Image
    check(fig, name)
    os.makedirs(folder, exist_ok=True)
    fig.savefig(os.path.join(folder, f"{name}.pdf"))
    png = os.path.join(folder, f"{name}.png")
    fig.savefig(png, dpi=600)
    im = Image.open(png); im.load()
    im.convert("RGB").save(os.path.join(folder, f"{name}.tif"), dpi=(600, 600), compression="tiff_lzw")
    im.save(png, dpi=(600, 600))
    plt.close(fig)
