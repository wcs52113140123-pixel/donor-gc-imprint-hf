#!/usr/bin/env python
"""Supplemental Figure S1 (JACC: BTS v0.5). In vitro GC class of the genes in the 5 frozen macrophage modules
(content of old Figure 1B). 183 x 78 mm, panels A-B. Standalone; reads only the input artifacts below.
Writes Supplemental_Figure_S1.{pdf,png,tif}, _source_data.xlsx, _number_registry.csv.
"""
import os, tarfile
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties

import os  # v0.5 package: inputs resolved relative to this script (override with JBTS_INPUTS)
HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.environ.get("JBTS_INPUTS", os.path.join(HERE, "inputs"))
STYLE = os.path.join(HERE, "jbts_style.py")
TAR_FI = os.path.join(INP, "figure_inputs_v05.tar.gz")
SD04 = os.path.join(INP, "Source_data_v04.xlsx")

exec(open(STYLE).read())
apply_jbts_style()
mpl.rcParams["savefig.pad_inches"] = 0
mpl.rcParams["hatch.linewidth"] = 0.5
IN = "_inputs"
if not os.path.isdir(f"{IN}/figure_inputs_v05"):
    os.makedirs(IN, exist_ok=True)
    with tarfile.open(TAR_FI) as t:
        t.extractall(IN)
SR = f"{IN}/figure_inputs_v05/04_results/validation_v1.2/steroid_reference/"
FS_T, FS_A, FS_S = SIZES["panel_title"], SIZES["tick"], SIZES["min_gene_label"]
REG, SRC = [], {}


def reg(panel, text, value, source, computation, match):
    REG.append(dict(figure="Supplemental Figure S1", panel=panel, printed_text=text, value=value,
                    source_file=source, computation=computation, matches_locked_text=match))


def tw(s, size=FS_S, weight="normal", style="normal"):
    fp = FontProperties(family=mpl.rcParams["font.sans-serif"][0], size=size, weight=weight, style=style)
    return TextPath((0, 0), s, prop=fp).get_extents().width / 72 * 25.4


# ------------------------------------------------------------------ data
gcc = pd.read_csv(SR + "gc_gene_classification.csv")
MODS = ["RES_EFF_core", "RES_EFF_lit", "MHCII_core", "ISG_IFN", "MAC_PAN_ID"]
NAMES = {"RES_EFF_core": "Resident efferocytic (core)", "RES_EFF_lit": "Efferocytosis (literature)",
         "MHCII_core": "MHC class II", "ISG_IFN": "Interferon-stimulated", "MAC_PAN_ID": "Pan-macrophage identity"}
CLS = ["GC_induced", "GC_repressed", "GC_insensitive", "ambiguous", "not_measured"]
LAB = {"GC_induced": "GC-induced", "GC_repressed": "GC-repressed", "GC_insensitive": "GC-insensitive",
       "ambiguous": "Ambiguous", "not_measured": "Not measured in ≥2 datasets"}
FC = {"GC_induced": GC_CLASS["induced"], "GC_repressed": GC_CLASS["repressed"],
      "GC_insensitive": GC_CLASS["insensitive"], "ambiguous": GC_CLASS["ambiguous"],
      "not_measured": GC_CLASS["not_measured"]}
g = gcc[gcc.module.isin(MODS)].copy()
g["cls"] = g["class"].fillna("not_measured").replace({"not_measured_in_>=2": "not_measured"})
ct = g.groupby(["module", "cls"]).size().unstack(fill_value=0).reindex(index=MODS, columns=CLS, fill_value=0)
# cross-check against the v0.4 source data for old Fig 1B
old = pd.read_excel(SD04, "Fig1B").set_index("module").rename(columns={"not_measured_in_>=2": "not_measured"})
assert (old.reindex(index=MODS, columns=CLS, fill_value=0).values == ct.values).all()
assert ct.loc["RES_EFF_lit"].sum() == 32 and ct.loc["RES_EFF_lit", "not_measured"] == 2

# ------------------------------------------------------------------ canvas
W, H = WIDTH_MM, 78.0
fig = plt.figure(figsize=(W * MM, H * MM))
fig.patches.append(Rectangle((0, 0), 1, 1, transform=fig.transFigure, fc="white", ec="none", zorder=-10))


def axmm(x, y, w, h):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


def header(x, y, w, h, letter, title):
    fr = axmm(x, y, w, h); fr.axis("off")
    panel_label(fr, letter, x=0, y=1 - 4.6 / h)
    fr.text(5.6 / w, 1 - 0.75 / h, title, transform=fr.transAxes, fontsize=FS_T, fontweight="bold", va="top",
            ha="left")


# ======================================================================= A  class counts per module
AX, AW = 0.0, 70.0
header(AX, 0, AW, H, "A", "In vitro GC class by module")
ax = axmm(AX + 36.0, 10.0, AW - 36.0 - 1.0, H - 10.0 - 9.6)
yy = np.arange(len(MODS))
left = np.zeros(len(MODS))
for c in CLS:
    v = ct[c].values
    ax.barh(yy, v, left=left, height=0.66, color=FC[c], lw=0.6,
            edgecolor=AMBIGUOUS_EDGE if c in ("ambiguous", "not_measured") else FC[c],
            hatch="////" if c == "not_measured" else None, zorder=2)
    for i, (l_, n_) in enumerate(zip(left, v)):
        if n_ >= 2:
            ax.text(l_ + n_ / 2, i, str(n_), fontsize=FS_S, ha="center", va="center",
                    color="white" if c in ("GC_induced", "GC_repressed") else TEXT, zorder=3)
    left += v
ax.set_yticks(yy)
ylabs = []
for m in MODS:
    n_tot = int(ct.loc[m].sum()); n_meas = n_tot - int(ct.loc[m, "not_measured"])
    s = f"{NAMES[m]}\n{n_meas} of {n_tot} measured" if n_meas < n_tot else f"{NAMES[m]}\n{n_tot} genes"
    ylabs.append(s)
    reg("A", s.replace("\n", "; "), f"{n_meas}/{n_tot}", "gc_gene_classification.csv",
        "genes per module; measured = classified in >=2 in vitro datasets",
        "yes" if m == "RES_EFF_lit" else "new")
ax.set_yticklabels(ylabs, fontsize=FS_S, linespacing=1.0)
ax.set_ylim(len(MODS) - 0.45, -0.75)
ax.set_xlim(0, 32.5); ax.set_xticks([0, 10, 20, 30])
ax.set_xlabel("Genes (n)", labelpad=1.0)
ax.tick_params(axis="y", length=0, pad=1.5)
for m in MODS:
    for c in CLS:
        if ct.loc[m, c] > 0:
            reg("A", f"{NAMES[m]}: {ct.loc[m, c]} {LAB[c]}", int(ct.loc[m, c]), "gc_gene_classification.csv",
                "count of genes per class", "yes" if (m, c) in {("RES_EFF_core", "GC_induced"),
                                                                 ("MHCII_core", "GC_repressed"),
                                                                 ("MHCII_core", "ambiguous")} else "new")
SRC["S1A_class_counts"] = ct.reset_index().rename(columns={"module": "module_id"}).assign(
    module=lambda d: d.module_id.map(NAMES))

# ======================================================================= B  gene tiles
BX = AW + 3.0; BW = W - BX
header(BX, 0, BW, H, "B", "Genes in each module, colored by in vitro class")
b = axmm(BX, 0, BW, H); b.set_xlim(0, BW); b.set_ylim(H, 0); b.axis("off")
ORDER = {c: i for i, c in enumerate(CLS)}
y0, row_h, xpad = 9.6, 4.15, 1.2
rec = []
y = y0
for m in MODS:
    t = g[g.module == m].assign(o=lambda d: d.cls.map(ORDER)).sort_values(["o", "meta_log2FC"],
                                                                         ascending=[True, False])
    b.text(xpad, y + 0.6, NAMES[m], fontsize=FS_S, fontweight="bold", va="center", ha="left")
    y += row_h * 0.9
    x = xpad
    for _, r_ in t.iterrows():
        wbox = tw(r_.gene, FS_S, style="italic") + 1.2
        if x + wbox > BW - 0.4:
            x = xpad; y += row_h
        c = r_.cls
        b.add_patch(Rectangle((x, y - 1.75), wbox, 3.5, fc=FC[c], lw=0.6,
                              ec=AMBIGUOUS_EDGE if c in ("ambiguous", "not_measured") else "none",
                              hatch="////" if c == "not_measured" else None, zorder=2))
        b.text(x + wbox / 2, y + 0.05, r_.gene, fontsize=FS_S, fontstyle="italic", ha="center", va="center",
               color="white" if c in ("GC_induced", "GC_repressed") else TEXT, zorder=3)
        rec.append(dict(module=NAMES[m], gene=r_.gene, gc_class=c, meta_log2FC=r_.meta_log2FC, k=r_.k))
        x += wbox + 0.45
    y += row_h * 1.05
assert y < H - 7.6, y
# shared key (bottom right)
xk = xpad; yk = H - 5.6
for c in CLS:
    if c == "ambiguous":
        xk = xpad; yk += 3.4
    b.add_patch(Rectangle((xk, yk - 1.1), 2.2, 2.2, fc=FC[c], lw=0.6,
                          ec=AMBIGUOUS_EDGE if c in ("ambiguous", "not_measured") else "none",
                          hatch="////" if c == "not_measured" else None))
    b.text(xk + 2.8, yk, LAB[c], fontsize=FS_S, va="center", ha="left")
    xk += 2.8 + tw(LAB[c]) + 2.2
assert xk - 2.2 < BW, xk
SRC["S1B_genes"] = pd.DataFrame(rec)

# ------------------------------------------------------------------ verify
fig.canvas.draw(); r = fig.canvas.get_renderer(); k = 25.4 / fig.dpi
items = [(t.get_text()[:30], t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
         if t.get_text().strip() and t.get_visible()]
out = [s for s, bb in items if bb.x0 * k < -0.15 or bb.x1 * k > W + 0.15 or bb.y0 * k < -0.15 or bb.y1 * k > H + 0.15]
ov = [(items[i][0], items[j][0]) for i in range(len(items)) for j in range(i + 1, len(items))
      if (min(items[i][1].x1, items[j][1].x1) - max(items[i][1].x0, items[j][1].x0)) * k > 0.15
      and (min(items[i][1].y1, items[j][1].y1) - max(items[i][1].y0, items[j][1].y0)) * k > 0.15]
print("outside:", out, "overlaps:", ov)
sizes = sorted({round(t.get_fontsize(), 1) for t in fig.findobj(mpl.text.Text) if t.get_text().strip()})
print("font sizes:", sizes); assert min(sizes) >= 8

# ------------------------------------------------------------------ save
stem = "Supplemental_Figure_S1"
save_figure(fig, stem)
from PIL import Image as _Im
_Im.open(f"{stem}.tif").convert("RGB").save(f"{stem}.tif", compression="tiff_lzw", dpi=(600, 600))
with pd.ExcelWriter(f"{stem}_source_data.xlsx") as xw:
    for kk, v in SRC.items():
        v.to_excel(xw, sheet_name=kk[:31], index=False)
pd.DataFrame(REG).to_csv(f"{stem}_number_registry.csv", index=False)
print("saved", fig.get_size_inches() / MM)
