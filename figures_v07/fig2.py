#!/usr/bin/env python
"""Figure 2 (JACC: BTS v0.7; reader-facing wording: dataset check, met / not met). Resident efferocytic macrophage genes are glucocorticoid-inducible in human
monocytes and macrophages in vitro.  183 x 240 mm, 8 panels (A-H).

Standalone: reads only the input artifacts listed below; writes Figure_2.{pdf,png,tif},
Figure_2_source_data.xlsx and Figure_2_number_registry.csv. All plotted values come from the files.
"""
import os, tarfile, itertools
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Polygon, Circle
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
import matplotlib.patheffects as pe

import os  # public package: inputs resolved relative to this script (override with JBTS_INPUTS)
HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.environ.get("JBTS_INPUTS", os.path.join(HERE, "..", "figures_v05", "inputs"))
STYLE = os.path.join(HERE, "jbts_style.py")
TAR_FI = os.path.join(INP, "figure_inputs_v05.tar.gz")
XL = os.path.join(INP, "Supplemental_Tables_S1-S14.xlsx")

exec(open(STYLE).read())
FAM = apply_jbts_style()
mpl.rcParams["savefig.pad_inches"] = 0
mpl.rcParams["hatch.linewidth"] = 0.5

IN = "_inputs"
if not os.path.isdir(f"{IN}/figure_inputs_v05"):
    os.makedirs(IN, exist_ok=True)
    with tarfile.open(TAR_FI) as t:
        t.extractall(IN)
SR = f"{IN}/figure_inputs_v05/04_results/validation_v1.2/steroid_reference/"

FS_T, FS_A, FS_S = SIZES["panel_title"], SIZES["tick"], SIZES["min_any"]   # 11 / 10 / 10 pt (v0.6: >=10 pt)
REG, SRC = [], {}


def reg(panel, text, value, source, computation, match):
    REG.append(dict(figure="Figure 2", panel=panel, printed_text=text, value=value,
                    source_file=source, computation=computation, matches_locked_text=match))


def tw(s, size=FS_S, weight="normal", style="normal"):
    fp = FontProperties(family=mpl.rcParams["font.sans-serif"][0], size=size, weight=weight, style=style)
    return TextPath((0, 0), s, prop=fp).get_extents().width / 72 * 25.4


# ------------------------------------------------------------------ data
DS = ["GSE85333", "GSE112101", "E-MEXP-2598", "E-MTAB-7572"]
S2 = pd.read_excel(XL, "S2_datasets").set_index("accession")
S4b = pd.read_excel(XL, "S4b_invitro_QC").set_index("dataset")
S5 = pd.read_excel(XL, "S5_GC_classification")
S6 = pd.read_excel(XL, "S6_leave_one_out").set_index("gene")
S3 = pd.read_excel(XL, "S3_gene_sets")
pdl = pd.read_csv(SR + "per_dataset_log2fc.csv.gz")
meta = pd.read_csv(SR + "meta_all_genes.csv.gz")
gcc = pd.read_csv(SR + "gc_gene_classification.csv")
assert (S5.set_index("gene")["class"].sort_index() == meta.set_index("gene")["class"].sort_index()).all()
cls = meta.set_index("gene")["class"]
KEY = {"GC_induced": "induced", "GC_repressed": "repressed", "GC_insensitive": "insensitive", "ambiguous": "ambiguous"}
counts = meta["class"].value_counts()
assert counts.to_dict() == {"GC_insensitive": 8750, "ambiguous": 3716, "GC_repressed": 508, "GC_induced": 394}

# ------------------------------------------------------------------ canvas / grid
W, H = WIDTH_MM, 240.0
GUT = 3.0
ROWS = [52.0, 52.0, 63.0, 64.0]
assert abs(sum(ROWS) + 3 * GUT - H) < 1e-9
fig = plt.figure(figsize=(W * MM, H * MM))
colw = (W - 11 * GUT) / 12
gs = fig.add_gridspec(4, 12, left=0, right=1, bottom=0, top=1, height_ratios=ROWS,
                      wspace=GUT / colw, hspace=GUT / np.mean(ROWS))
fig.patches.append(Rectangle((0, 0), 1, 1, transform=fig.transFigure, fc="white", ec="none", zorder=-10))


def cell(r0, r1, c0, c1):
    bb = gs[r0:r1, c0:c1].get_position(fig)
    return dict(x=bb.x0 * W, y=(1 - bb.y1) * H, w=bb.width * W, h=bb.height * H)


def axmm(x, y, w, h, **kw):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **kw)


def cellx(r0, r1, x, w):
    """Cell spanning grid rows r0..r1 at a custom x/width (mm); used where the 12-column grid is too coarse."""
    bb = gs[r0:r1, 0:12].get_position(fig)
    return dict(x=x, y=(1 - bb.y1) * H, w=w, h=bb.height * H)


# bottom block (rows 2-3): D/E 56 mm | F 63 mm | G 25 mm | H 30 mm, 3-mm gutters -> 183 mm
BX = dict(DE=(0.0, 56.0), F=(59.0, 63.0), G=(125.0, 25.0), H=(153.0, 30.0))
assert abs(BX["H"][0] + BX["H"][1] - W) < 1e-9


def canvas(c):
    ax = axmm(c["x"], c["y"], c["w"], c["h"])
    ax.set_xlim(0, c["w"]); ax.set_ylim(c["h"], 0); ax.axis("off")
    return ax


def header(c, letter, title, title_dx=5.6):
    fr = axmm(c["x"], c["y"], c["w"], c["h"]); fr.axis("off")
    panel_label(fr, letter, x=0, y=1 - 4.6 / c["h"])
    fr.text(title_dx / c["w"], 1 - 0.75 / c["h"], title, transform=fr.transAxes, fontsize=FS_T,
            fontweight="bold", va="top", ha="left", linespacing=1.05)
    return fr


def t8(ax, x, y, s, **kw):
    kw.setdefault("fontsize", FS_S); kw.setdefault("va", "center"); kw.setdefault("ha", "left")
    kw.setdefault("linespacing", 1.05)
    return ax.text(x, y, s, **kw)


def rbox(ax, x, y, w, h, fc="white", ec="#7F7F7F", lw=0.6, ls="-", r=1.0, z=1):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw, ls=ls,
                       zorder=z)
    ax.add_patch(p); return p


def arrow(ax, x0, y0, x1, y1, color="#555555", lw=0.9, ms=6):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                mutation_scale=ms, shrinkA=0, shrinkB=0), zorder=3)


def cls_style(k):
    """face, edge for an in vitro class key."""
    if k == "ambiguous":
        return GC_CLASS["ambiguous"], AMBIGUOUS_EDGE
    return GC_CLASS[k], "none"


def swatch(ax, x, y, k, s=2.2):
    fc, ec = cls_style(k) if k != "not_measured" else (GC_CLASS["not_measured"], "#9A9A9A")
    ax.add_patch(Rectangle((x, y - s / 2), s, s, fc=fc, ec=ec, lw=0.6,
                           hatch="////" if k == "not_measured" else None, zorder=4))


# piecewise-linear -log10 P axis: linear to BRK, compressed by FACTOR above (marked on the axis)
def make_scale(brk, factor):
    def fwd(y):
        y = np.asarray(y, float); return np.where(y <= brk, y, brk + (y - brk) * factor)
    def inv(y):
        y = np.asarray(y, float); return np.where(y <= brk, y, brk + (y - brk) / factor)
    return fwd, inv


def mark_break(ax, brk, fwd):
    """Double slash on the y spine where the -log10 P scale changes slope."""
    y0, y1 = ax.get_ylim()
    yf = float((fwd(brk) - fwd(y0)) / (fwd(y1) - fwd(y0)))
    ax.plot([0, 0], [yf - 0.018, yf + 0.018], transform=ax.transAxes, color="white", lw=1.6, clip_on=False, zorder=6)
    for off in (-0.018, 0.018):
        ax.plot([-0.035, 0.035], [yf + off - 0.012, yf + off + 0.012], transform=ax.transAxes, color=TEXT, lw=0.6,
                clip_on=False, zorder=7)


# ======================================================================= A  GC reference schematic
cA = cell(0, 1, 0, 7)
header(cA, "A", "Four datasets meeting QC define the GC reference")
a = canvas(cA)
info = {"GSE85333": ("MDM", "prednisolone", "24 h"),
        "GSE112101": ("Monocytes", "methylprednisolone", "6 h"),
        "E-MEXP-2598": ("Monocytes", "dexamethasone", "36 h"),
        "E-MTAB-7572": ("MDM", "dexamethasone", "24 h")}
donors = {d: int(str(S2.loc[d, "samples"]).split()[0]) for d in DS}
assert [donors[d] for d in DS] == [6, 4, 4, 3]
NUM = {d: str(i + 1) for i, d in enumerate(DS)}          # dataset keys 1-4, reused in C, G and H


def badge(ax, x, y, s, r=1.75):
    ax.add_patch(Circle((x, y), r, fc="#555555", ec="none", zorder=5))
    ax.text(x, y + 0.05, s, fontsize=FS_S, color="white", fontweight="bold", ha="center", va="center", zorder=6)


cols = [("Dataset", [d for d in DS], "bold"), ("Cells", [info[d][0] for d in DS], "normal"),
        ("GC (vs vehicle)", [info[d][1] for d in DS], "normal"), ("Time", [info[d][2] for d in DS], "normal"),
        ("Donors", [str(donors[d]) for d in DS], "normal")]
cxs = [5.0]
for (lab, vals, wt) in cols[:-1]:
    cxs.append(cxs[-1] + max([tw(lab, FS_S, "bold")] + [tw(v, FS_S, wt) for v in vals]) + 2.2)
x_tab_end = cxs[-1] + tw("Donors", FS_S, "bold")
assert x_tab_end < cA["w"] - 0.5, x_tab_end
ty0, LP = 8.4, 4.3
for k, (lab, vals, wt) in enumerate(cols):
    t8(a, cxs[k], ty0, lab, fontweight="bold", color="#444444")
a.plot([0.4, x_tab_end + 0.4], [ty0 + 2.3, ty0 + 2.3], color=TEXT, lw=0.6)
recA = []
for i, d in enumerate(DS):
    yy = ty0 + 4.7 + i * LP
    badge(a, 2.2, yy, NUM[d])
    for k, (lab, vals, wt) in enumerate(cols):
        t8(a, cxs[k], yy, vals[i], fontweight=wt)
    cell_, agent, time = info[d]
    recA.append(dict(key=NUM[d], dataset=d, cells=cell_, agent=agent, time=time, donors=donors[d],
                     design=S4b.loc[d, "design"], dataset_check_met=bool(S4b.loc[d, "pass"])))
    reg("A", f"{NUM[d]} {d}: {agent}, {time}, {donors[d]} donors", donors[d], "S2_datasets; S4b_invitro_QC",
        "as tabulated", "yes (agent/cell, T); donors new (S)")
# meta-analysis line: funnel icon + text
ym = ty0 + 4.7 + 3 * LP + 4.7
fx = 3.0
a.add_patch(Polygon([[fx - 2.4, ym - 2.0], [fx + 2.4, ym - 2.0], [fx + 0.7, ym + 0.5], [fx + 0.7, ym + 2.0],
                     [fx - 0.7, ym + 2.0], [fx - 0.7, ym + 0.5]], closed=True, fc="#DADADA", ec=TEXT, lw=0.6,
                    zorder=3))
s1 = "Random-effects meta-analysis"
t8(a, 6.6, ym, s1)
arrow(a, 6.6 + tw(s1) + 1.6, ym, 6.6 + tw(s1) + 5.2, ym)
t8(a, 6.6 + tw(s1) + 6.0, ym, "13,368 genes classified", fontweight="bold")
assert 6.6 + tw(s1) + 6.0 + tw("13,368 genes classified", FS_S, "bold") < cA["w"] - 0.5
reg("A", "13,368 genes classified", 13368, "meta_all_genes.csv.gz", "row count", "yes")
# classification rule strip (schematic number line, maps -1.2..1.2)
ry0, RP = ym + 4.2, 4.2
rules = [("induced", "GC-induced", "≥0.5, P<0.01, ≥2/3 same sign", (0.78, 0.56, 1.0)),
         ("repressed", "GC-repressed", "≤−0.5, P<0.01, ≥2/3 same sign", (-0.78, -1.0, -0.56)),
         ("insensitive", "GC-insensitive", "95% CI within ±0.5", (0.05, -0.3, 0.4)),
         ("ambiguous", "Ambiguous", "none of the above", (0.5, -0.15, 1.1))]
xl0, xl1 = 1.0, 27.0
def X(v): return xl0 + (v + 1.2) / 2.4 * (xl1 - xl0)
ytk = ry0 + 3 * RP + 3.8
for v in (-0.5, 0.5):
    a.plot([X(v), X(v)], [ry0 - 1.6, ytk - 2.0], ls=(0, (1, 1.2)), color="#7F7F7F", lw=0.6, zorder=1)
a.plot([X(0), X(0)], [ry0 - 1.6, ytk - 2.0], color="#BDBDBD", lw=0.5, zorder=1)
xr = xl1 + 2.0 + max(tw(r_[1], FS_S, "bold") for r_ in rules) + 2.0
for i, (k, name, rule, (m_, lo, hi)) in enumerate(rules):
    yy = ry0 + i * RP
    a.plot([X(lo), X(hi)], [yy, yy], color=TEXT if k != "insensitive" else "#7F7F7F", lw=0.8, zorder=3)
    fc, ec = cls_style(k)
    a.add_patch(Circle((X(m_), yy), 1.0, fc=fc, ec=ec if ec != "none" else TEXT, lw=0.6, zorder=4))
    t8(a, xl1 + 2.0, yy, name, fontweight="bold")
    t8(a, xr, yy, rule)
    assert xr + tw(rule) < cA["w"] - 0.5
for v, s in [(-0.5, "−0.5"), (0, "0"), (0.5, "0.5")]:
    t8(a, X(v), ytk, s, ha="center", color="#444444")
t8(a, xl1 + 2.0, ytk, "Meta log2FC (schematic)", color="#555555", fontstyle="italic")
assert ytk + 2.1 < cA["h"] + 2.2, ytk
for s_, v_ in [("log2FC ≥0.5", 0.5), ("P<0.01", 0.01), ("≥2/3 same sign", "2/3"), ("95% CI within ±0.5", 0.5)]:
    reg("A", s_, v_, "Manuscript v0.4 Methods (In vitro GC reference)", "rule as stated", "yes")
SRC["2A_datasets"] = pd.DataFrame(recA)

# ======================================================================= B  positive-control check
cB = cell(0, 1, 7, 12)
header(cB, "B", "Dataset check: induced in 4 of 6")
QC_ORDER = DS + ["GSE150049", "GSE35683"]
B_LW = max(tw(d) for d in DS + ["GSE150049", "GSE35683"]) + 2.4
B_RW = max(tw("FDR 0.21–0.24"), tw("4.03 / 2.97")) + 2.0
axB = axmm(cB["x"] + B_LW, cB["y"] + 16.0, cB["w"] - B_LW - B_RW, cB["h"] - 16.0 - 10.4)
recB = []
for i, d in enumerate(QC_ORDER):
    r = S4b.loc[d]; ok = bool(r["pass"])
    for gname, mk_, dy in [("FKBP5", "o", -0.24), ("TSC22D3", "s", 0.24)]:
        lf, fdr = r[f"{gname}_logFC"], r[f"{gname}_FDR"]
        filled = fdr < 0.05
        recB.append(dict(dataset=d, gene=gname, log2FC=lf, FDR=fdr, filled_FDR_lt_0_05=filled, dataset_check_met=ok))
        axB.plot(lf, i + dy, mk_, ms=4.2, mfc=TEXT if filled else "white", mec=TEXT if ok else "#8C8C8C",
                 mew=0.8, zorder=3)
    if ok:
        s = f"{r['FKBP5_logFC']:.2f} / {r['TSC22D3_logFC']:.2f}"
        axB.text(1.02, i, s, transform=axB.get_yaxis_transform(), fontsize=FS_S, va="center", ha="left")
        reg("B", f"{d} FKBP5 / TSC22D3 {s}", s, "S4b_invitro_QC", "logFC rounded to 2 dp", "new (S: S-R2)")
axB.text(1.0 + B_RW / axB.get_position().width / W, -0.65, "FKBP5 / TSC22D3", transform=axB.get_yaxis_transform(),
         fontsize=FS_S, va="bottom", ha="right", color="#555555", fontstyle="italic")
r = S4b.loc["GSE150049"]
fdr_rng = f"FDR {min(r.FKBP5_FDR, r.TSC22D3_FDR):.2f}–{max(r.FKBP5_FDR, r.TSC22D3_FDR):.2f}"
axB.text(1.02, 4, fdr_rng, transform=axB.get_yaxis_transform(), fontsize=FS_S, va="center",
         ha="left", color="#555555", linespacing=1.0)
reg("B", f"GSE150049 not met: {fdr_rng}", fdr_rng, "S4b_invitro_QC", "min-max of the two FDRs", "new (S: S-R2)")
axB.text(1.02, 5, "stimulated", transform=axB.get_yaxis_transform(), fontsize=FS_S, va="center",
         ha="left", color="#555555")
axB.axvline(1.0, ls=(0, (3, 2)), color="#7F7F7F", lw=0.6, zorder=1)
axB.axhline(3.5, color=GRID, lw=0.6, zorder=0)
axB.set_yticks(range(len(QC_ORDER))); axB.set_yticklabels(QC_ORDER, fontsize=FS_S)
for tl, d in zip(axB.get_yticklabels(), QC_ORDER):
    if d not in DS:
        tl.set_color("#7F7F7F")
axB.set_ylim(len(QC_ORDER) - 0.5, -0.6)
axB.set_xlim(0, 4.6); axB.set_xticks([0, 1, 2, 3, 4])
axB.set_xlabel("log2 fold change", labelpad=1.0)
axB.tick_params(axis="y", length=0, pad=1.5)
kb = canvas(cB)
xx = 1.4
for mk_, lab, fc_, st in [("o", "FKBP5", TEXT, "italic"), ("s", "TSC22D3", TEXT, "italic"), ("o", "FDR ≥0.05", "white", "normal")]:
    kb.plot(xx, 8.4, mk_, ms=4.6, mfc=fc_, mec=TEXT, mew=0.8)
    t8(kb, xx + 2.0, 8.4, lab, fontstyle=st)
    xx += 2.0 + tw(lab, FS_S, style=st) + 3.2
assert xx < cB["w"]
SRC["2B_dataset_check"] = pd.DataFrame(recB)

# ======================================================================= C  per-dataset volcanoes
cC = cell(1, 2, 0, 12)
header(cC, "C", "The GC response replicates in each dataset")
kc = canvas(cC)
# class key (first panel using the class palette; shared by C, D, F, H)
xx = cC["w"] - 1.0
items = [("induced", "Induced"), ("repressed", "Repressed"), ("insensitive", "Insensitive"),
         ("ambiguous", "Ambiguous"), ("unclassified", "Not classified")]
widths = [3.3 + tw(l) + 2.6 for _, l in items]
xx = 5.6                                   # second line, under the panel title
KY = 7.0
for (k, lab), wdt in zip(items, widths):
    if k == "unclassified":
        kc.add_patch(Circle((xx + 1.1, KY), 0.85, fc="#E4E4E4", ec="none", zorder=4))
    else:
        fc, ec = cls_style(k)
        kc.add_patch(Circle((xx + 1.1, KY), 1.0, fc=fc, ec=ec if ec != "none" else fc, lw=0.6, zorder=4))
    t8(kc, xx + 2.9, KY, lab)
    xx += wdt
LAB = ["FKBP5", "TSC22D3", "ZBTB16", "CD163", "MERTK", "HLA-DQA1", "CIITA"]
BRK, FAC = 15.0, 0.25
fwd, inv = make_scale(BRK, FAC)
XL_, XR_ = -10.0, 11.0
LEADER_CROSSINGS = {}
pad_l, gap = 13.0, 2.0
fw = (cC["w"] - pad_l - 3 * gap - 0.5) / 4
yv0, hv = 15.0, cC["h"] - 15.0 - 10.0
order = ["GC_insensitive", "unclassified", "ambiguous", "GC_repressed", "GC_induced"]
recC = []
axes_C = []
for j, d in enumerate(DS):
    ax = axmm(cC["x"] + pad_l + j * (fw + gap), cC["y"] + yv0, fw, hv)
    axes_C.append(ax)
    t = pdl[pdl.dataset == d].copy()
    t["class"] = t.gene.map(cls).fillna("unclassified")
    t["nlp"] = -np.log10(t.P)
    t["x_plot"] = t.logFC.clip(XL_ + 0.15, XR_ - 0.15)
    t["clamped"] = (t.logFC < XL_) | (t.logFC > XR_)
    ax.set_yscale("function", functions=(fwd, inv))
    for c_ in order:
        s = t[t["class"] == c_]
        if c_ == "unclassified":
            ax.scatter(s.x_plot, s.nlp, s=0.35, c="#ECECEC", lw=0, rasterized=True, zorder=0)
        elif c_ == "ambiguous":
            ax.scatter(s.x_plot, s.nlp, s=1.8, facecolors="white", edgecolors="#6A6A6A", lw=0.35,
                       rasterized=True, zorder=2)
        else:
            k = KEY[c_]
            ax.scatter(s.x_plot, s.nlp, s=1.4 if k == "insensitive" else 2.4, c=GC_CLASS[k], lw=0,
                       rasterized=True, zorder=1 if k == "insensitive" else 3)
    cl = t[t.clamped]
    if len(cl):
        ax.scatter(cl.x_plot, cl.nlp, marker=">", s=6, c=[GC_CLASS[KEY.get(c_, "insensitive")] if c_ in KEY else
                   "#BDBDBD" for c_ in cl["class"]], lw=0.3, edgecolors=TEXT, zorder=4)
    for v in (-0.5, 0.5):
        ax.axvline(v, ls=(0, (1, 1.2)), color="#9A9A9A", lw=0.5, zorder=0)
    ax.set_xlim(XL_, XR_); ax.set_ylim(0, 35)
    ax.set_yticks([0, 5, 10, 15, 35]); ax.set_xticks([-5, 0, 5, 10])
    mark_break(ax, BRK, fwd)
    if j == 0:
        ax.set_ylabel("−log10 P", labelpad=1.0)
    else:
        ax.set_yticklabels([])
    ax.set_xlabel("log2 fold change" if j in (1, 2) else "", labelpad=1.0)
    # facet header and n genes
    cell_, agent, time = info[d]
    badge(kc, pad_l + j * (fw + gap) + 1.75, yv0 - 2.8, NUM[d])     # agent per dataset: panel A table
    t8(kc, pad_l + j * (fw + gap) + 4.4, yv0 - 2.75, d, fontweight="bold")
    ng = len(t)
    assert ng == S4b.loc[d, "n_genes"]
    ax.text(0.03, 0.975, f"{ng:,} genes", transform=ax.transAxes, fontsize=FS_S, ha="left", va="top",
            color="#444444")
    reg("C", f"{d}: {ng:,} genes", ng, "per_dataset_log2fc.csv.gz (= S4b n_genes)", "rows per dataset",
        "new (S: S4b n_genes)")
    # labels: induced genes in a right-hand column, repressed genes in a left-hand column below the axis break.
    # Slot order is chosen by exhaustive search to minimise leader crossings, then total leader length.
    lab = t[t.gene.isin(LAB)].copy()
    lab["yy"] = fwd(lab.nlp)
    ux = (XR_ - XL_) / fw                                   # data units per mm (x)
    uy = float(fwd(35)) / hv                                # forward-y units per mm
    for grp, xcol, ha, (ylo, yhi) in [(lab[lab.logFC > 0], XR_ - 0.25, "right", (float(fwd(1.6)), float(fwd(35)) - 1.3 * 3.1 * uy)),
                                      (lab[lab.logFC < 0], XL_ + 0.25, "left", (1.4, BRK - 1.6))]:
        k = len(grp)
        if k == 0:
            continue
        step = 4.3 * uy                                     # 3.25 mm between label baselines
        c = float(np.clip(grp.yy.mean(), ylo + (k - 1) * step / 2, yhi - (k - 1) * step / 2))
        slots = [c + (j - (k - 1) / 2) * step for j in range(k)]
        pts = [(float(r_.x_plot), float(r_.yy), r_.gene) for _, r_ in grp.iterrows()]
        def ends(g):
            w_ = (tw(g, FS_S, style="italic") + 0.6) * ux
            return xcol - w_ if ha == "right" else xcol + w_
        def seg_cross(p1, p2, p3, p4):
            def o(a_, b_, c_): return (b_[0] - a_[0]) * (c_[1] - a_[1]) - (b_[1] - a_[1]) * (c_[0] - a_[0])
            return (o(p1, p2, p3) * o(p1, p2, p4) < 0) and (o(p3, p4, p1) * o(p3, p4, p2) < 0)
        best = None
        for perm in itertools.permutations(range(k)):
            segs = [((ends(pts[i_][2]) / ux, slots[perm[i_]] / uy), (pts[i_][0] / ux, pts[i_][1] / uy)) for i_ in range(k)]
            ncross = sum(seg_cross(*segs[a_], *segs[b_]) for a_ in range(k) for b_ in range(a_ + 1, k))
            length = sum(np.hypot(p[0][0] - p[1][0], p[0][1] - p[1][1]) for p in segs)
            key = (ncross, length)
            if best is None or key < best[0]:
                best = (key, perm)
        for i_, (px, py, g) in enumerate(pts):
            ys_ = slots[best[1][i_]]
            ax.text(xcol, float(inv(ys_)), g, fontsize=FS_S, fontstyle="italic", ha=ha, va="center", zorder=8,
                    path_effects=[pe.withStroke(linewidth=2.2, foreground="white")],
                    bbox=dict(boxstyle="square,pad=0.04", fc="white", ec="none", alpha=0.85))
            ax.plot([ends(g), px], [float(inv(ys_)), float(inv(py))], color="#444444", lw=0.45, zorder=7)
            ax.plot(px, float(inv(py)), "o", ms=2.6, mfc="none", mec=TEXT, mew=0.6, zorder=7)
        LEADER_CROSSINGS[d] = LEADER_CROSSINGS.get(d, 0) + best[0][0]
    recC.append(t[["gene", "logFC", "P", "nlp", "class", "clamped"]].assign(dataset=d))
print("leader crossings per facet:", LEADER_CROSSINGS)
SRC["2C_volcano_per_dataset"] = pd.concat(recC).rename(columns={"nlp": "neg_log10_P", "class": "meta_class"})

# ======================================================================= D  genome-wide meta volcano
cD = cellx(2, 3, *BX["DE"])
header(cD, "D", "13,368 genes classified")
axD = axmm(cD["x"] + 13.0, cD["y"] + 21.6, cD["w"] - 13.0 - 1.2, cD["h"] - 21.6 - 10.4)
BRK_D, FAC_D = 20.0, 0.1
fwdD, invD = make_scale(BRK_D, FAC_D)
axD.set_yscale("function", functions=(fwdD, invD))
m = meta.copy(); m["nlp"] = -np.log10(m.P)
for c_ in ["GC_insensitive", "ambiguous", "GC_repressed", "GC_induced"]:
    s = m[m["class"] == c_]; k = KEY[c_]
    if k == "ambiguous":
        axD.scatter(s.meta_log2FC, s.nlp, s=2.0, facecolors="white", edgecolors="#6A6A6A", lw=0.35,
                    rasterized=True, zorder=2)
    else:
        axD.scatter(s.meta_log2FC, s.nlp, s=1.4 if k == "insensitive" else 2.6, c=GC_CLASS[k], lw=0,
                    rasterized=True, zorder=1 if k == "insensitive" else 3)
for v in (-0.5, 0.5):
    axD.axvline(v, ls=(0, (1, 1.2)), color="#9A9A9A", lw=0.5, zorder=0)
axD.set_xlim(-3.4, 6.6); axD.set_xticks([-2, 0, 2, 4, 6])
ymaxD = float(np.ceil(m.nlp.max() / 10) * 10)
axD.set_ylim(0, ymaxD + 5)
axD.set_yticks([0, 10, 20, 100, 190])
mark_break(axD, BRK_D, fwdD)
axD.set_xlabel("Meta log2 fold change", labelpad=1.0)
axD.set_ylabel("−log10 P", labelpad=1.0)
for gname, xt_, yt_, ha in [("CD163", 4.0, 30.0, "left"), ("TSC22D3", 0.9, 75.0, "left"),
                            ("FKBP5", 4.0, 14.0, "left"), ("ZBTB16", 6.4, 7.0, "right")]:
    rr = m[m.gene == gname].iloc[0]
    axD.annotate(gname, xy=(rr.meta_log2FC, rr.nlp), xytext=(xt_, yt_), fontsize=FS_S, fontstyle="italic", ha=ha,
                 va="center", zorder=6, bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none", alpha=0.8),
                 arrowprops=dict(arrowstyle="-", lw=0.45, color="#555555", shrinkA=0, shrinkB=1.2))
    axD.plot(rr.meta_log2FC, rr.nlp, "o", ms=2.6, mfc="none", mec=TEXT, mew=0.6, zorder=6)
# class counts (replaces the inset bar: two of four segments would be <1 mm wide)
kd = canvas(cD)
cls_order = ["GC_induced", "GC_repressed", "GC_insensitive", "ambiguous"]
vals = [int(counts[c_]) for c_ in cls_order]
names_ = {"GC_induced": "Induced", "GC_repressed": "Repressed", "GC_insensitive": "Insensitive", "ambiguous": "Ambiguous"}
xx, yy = 1.0, 8.6
for c_, v in zip(cls_order, vals):
    lab_ = f"{names_[c_]} {v:,}"
    wd = 3.3 + tw(lab_) + 2.6
    if xx + wd - 2.6 > cD["w"] - 0.3:
        xx, yy = 1.0, yy + 4.3
    swatch(kd, xx, yy, KEY[c_], s=2.6)
    t8(kd, xx + 3.3, yy, lab_)
    reg("D", lab_, v, "meta_all_genes.csv.gz / S5", "class counts", "yes")
    xx += wd
D_KEY_BOTTOM = yy + 2.2
assert D_KEY_BOTTOM <= 20.0, D_KEY_BOTTOM
SRC["2D_meta_volcano"] = m[["gene", "k", "meta_log2FC", "ci_lo", "ci_hi", "P", "nlp", "class"]].rename(
    columns={"nlp": "neg_log10_P"})
SRC["2D_inset_counts"] = pd.DataFrame({"class": cls_order, "n_genes": vals})

# ======================================================================= E  between-dataset concordance
cE = cellx(3, 4, *BX["DE"])
header(cE, "E", "Genome-wide agreement")
wide = pdl[pdl.dataset.isin(DS)].pivot_table(index="gene", columns="dataset", values="logFC")[DS]
labs2 = dict(NUM)              # diagonal: dataset key 1-4 (as in panel A)
mx0, my0 = 0.6, 7.6
mw = cE["w"] - mx0 - 12.6; mh = cE["h"] - my0 - 0.6
cw_, ch_ = mw / 4, mh / 4
rcmap = LinearSegmentedColormap.from_list("r_grey", ["#F4F4F4", "#9A9A9A", "#3A3A3A"])
rnorm = Normalize(0.4, 0.7)
recE = []
for i in range(4):
    for j in range(4):
        ax = axmm(cE["x"] + mx0 + j * cw_ + 0.3, cE["y"] + my0 + i * ch_ + 0.3, cw_ - 0.6, ch_ - 0.6)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_linewidth(0.4); sp.set_color("#9A9A9A")
        if i == j:
            ax.text(0.5, 0.5, labs2[DS[i]], transform=ax.transAxes, fontsize=FS_T, ha="center", va="center",
                    fontweight="bold")
            assert max(tw(z) for z in labs2[DS[i]].split("\n")) < cw_ - 0.8
            for sp in ax.spines.values():
                sp.set_visible(False)
        else:
            s = wide[[DS[j], DS[i]]].dropna()
            rr = float(s.corr().iloc[0, 1])
            if i > j:
                ax.set_facecolor(rcmap(rnorm(rr)))
            else:
                ax.hexbin(s[DS[j]].clip(-4, 6), s[DS[i]].clip(-4, 6), gridsize=14, extent=(-4, 6, -4, 6),
                          bins="log", cmap="Greys",
                          mincnt=1, lw=0.1, rasterized=True)
                ax.plot([-4, 6], [-4, 6], color="#9A9A9A", lw=0.4, ls=(0, (1, 1)), zorder=3)
                ax.set_xlim(-4, 6); ax.set_ylim(-4, 6)
            if i > j:
                recE.append(dict(dataset_row=DS[i], dataset_col=DS[j], n_shared_genes=len(s), pearson_r=rr))
cax = axmm(cE["x"] + cE["w"] - 10.4, cE["y"] + my0 + 4.6, 1.8, mh - 9.0)
cb = mpl.colorbar.ColorbarBase(cax, cmap=rcmap, norm=rnorm, orientation="vertical")
cb.set_ticks([0.4, 0.55, 0.7]); cb.ax.tick_params(labelsize=FS_S, length=1.5, pad=1)
cb.outline.set_linewidth(0.4)
cax.text(0.5, 1.03, "r", transform=cax.transAxes, ha="center", va="bottom", fontsize=FS_S, fontstyle="italic")
SRC["2E_concordance"] = pd.DataFrame(recE)

# ======================================================================= F/G/H  module genes
cF = cellx(2, 4, *BX["F"]); cG = cellx(2, 4, *BX["G"]); cH = cellx(2, 4, *BX["H"])
F_WORDS = "Most resident efferocytic genes are GC-induced (7 of 9); most MHC-II genes are GC-repressed (6 of 10)".split()
F_LINES, cur = [], ""
for w_ in F_WORDS:                                   # greedy wrap at the measured title width
    trial = (cur + " " + w_).strip()
    if tw(trial, FS_T, "bold") <= cF["w"] - 5.6 - 0.4:
        cur = trial
    else:
        F_LINES.append(cur); cur = w_
F_LINES.append(cur)
header(cF, "F", "\n".join(F_LINES))
header(cG, "G", "Per dataset")
header(cH, "H", "Stable class")
gs3 = S3.groupby("gene_set").gene.apply(list)
groups = [("GC\nscore", gs3["GC_specific3"]),
          ("Resident efferocytic", list(gcc[gcc.module == "RES_EFF_core"].sort_values("meta_log2FC",
                                                                                    ascending=False).gene)),
          ("MHC-II", list(gcc[gcc.module == "MHCII_core"].sort_values("meta_log2FC").gene))]
rows = []
y = 0.0
for gname_, genes in groups:
    for g_ in genes:
        rows.append((gname_, g_, y)); y += 1.0
    y += 0.55
nrow_span = y - 0.55
TOPZ, BOTZ = len(F_LINES) * 4.9 + 2.4, 11.4
yt0 = cF["y"] + TOPZ; hrow = cF["h"] - TOPZ - BOTZ
def ylim_(ax): ax.set_ylim(nrow_span - 0.5, -0.5)
ms5 = S5.set_index("gene")
PITCH = hrow / nrow_span
assert PITCH >= 4.15, PITCH                             # 10-pt row labels need >=4.15 mm per row
BRX = 10.6                                              # bracket x (mm); group labels left of it
F_X0 = BRX + 0.8 + max(tw(g_, FS_S, style="italic") for g_ in sum([gl for _, gl in groups], [])) + 1.6
# F forest
axF = axmm(cF["x"] + F_X0, yt0, cF["w"] - F_X0 - 0.6, hrow)
recF = []
for grp, g_, yy in rows:
    r_ = ms5.loc[g_]
    k = KEY[r_["class"]]
    fc, ec = cls_style(k)
    axF.plot([r_.ci_lo, r_.ci_hi], [yy, yy], color=TEXT, lw=0.8, zorder=2)
    axF.plot(r_.meta_log2FC, yy, "o", ms=4.4, mfc=fc, mec=TEXT if k in ("ambiguous",) else TEXT, mew=0.5, zorder=3)
    recF.append(dict(group=grp.replace("\n", " "), gene=g_, meta_log2FC=r_.meta_log2FC, ci_lo=r_.ci_lo,
                     ci_hi=r_.ci_hi, P=r_.P, k=r_.k, gc_class=r_["class"]))
for v in (-0.5, 0.5):
    axF.axvline(v, ls=(0, (1, 1.2)), color="#7F7F7F", lw=0.6, zorder=0)
axF.axvline(0, color="#BDBDBD", lw=0.5, zorder=0)
axF.set_xlim(-2.2, 7.6); axF.set_xticks([-2, 0, 2, 4, 6])
ylim_(axF); axF.set_yticks([yy for _, _, yy in rows])
axF.set_yticklabels([g_ for _, g_, _ in rows], fontsize=FS_S, fontstyle="italic")
axF.tick_params(axis="y", length=0, pad=1.2)
axF.spines["left"].set_visible(False)
axF.set_xlabel("Meta log2FC", labelpad=1.0)
# value labels for the two headline genes
for g_ in ["CD163", "LYVE1"]:
    yy = [r[2] for r in rows if r[1] == g_][0]; r_ = ms5.loc[g_]
    axF.text(r_.ci_hi + 0.25, yy, f"{r_.meta_log2FC:.2f}", fontsize=FS_S, va="center", ha="left")
    reg("F", f"{g_} {r_.meta_log2FC:.2f}", round(r_.meta_log2FC, 2), "S5_GC_classification", "meta_log2FC, 2 dp", "yes")
res = gcc[gcc.module == "RES_EFF_core"]; mhc = gcc[gcc.module == "MHCII_core"]
n_ind = int((res["class"] == "GC_induced").sum()); n_rep = int((mhc["class"] == "GC_repressed").sum())
n_amb = int((mhc["class"] == "ambiguous").sum())
assert (n_ind, len(res), n_rep, len(mhc), n_amb) == (7, 9, 6, 10, 4)
yr = [r[2] for r in rows if r[0].startswith("Resident")]
ym = [r[2] for r in rows if r[0] == "MHC-II"]
reg("F", "title: 7 of 9 (GC-induced)", "7/9", "gc_gene_classification.csv (RES_EFF_core)", "count class", "yes")
reg("F", "title: 6 of 10 (GC-repressed)", "6/10", "gc_gene_classification.csv (MHCII_core)", "count class", "yes")
# group brackets (left of gene labels)
fb = canvas(cF)
for grp, genes in groups:
    ys_ = [r[2] for r in rows if r[0] == grp]
    def yy_mm(v): return TOPZ + (v + 0.5) / nrow_span * hrow
    y_a, y_b = yy_mm(ys_[0] - 0.4), yy_mm(ys_[-1] + 0.4)
    fb.plot([BRX + 0.6, BRX, BRX, BRX + 0.6], [y_a, y_a, y_b, y_b], color=TEXT, lw=0.7)
    fb.text(BRX / 2, (y_a + y_b) / 2, grp, rotation=90, fontsize=FS_S, ha="center", va="center", fontweight="bold",
            linespacing=0.95, multialignment="center")
    assert max(tw(z, FS_S, "bold") for z in grp.split("\n")) < (y_b - y_a) + 0.5, grp
SRC["2F_forest"] = pd.DataFrame(recF)
# G heatmap: per-dataset log2FC
axG = axmm(cG["x"] + 0.6, yt0, cG["w"] - 1.2, hrow)
divc = diverging_cmap()
gnorm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
wide_all = wide
recG = []
for grp, g_, yy in rows:
    for j, d in enumerate(DS):
        v = wide_all.loc[g_, d] if g_ in wide_all.index else np.nan
        if np.isnan(v):
            axG.add_patch(Rectangle((j - 0.5, yy - 0.5), 1, 1, fc=GC_CLASS["not_measured"], ec="#9A9A9A", lw=0.3,
                                    hatch="////"))
        else:
            axG.add_patch(Rectangle((j - 0.5, yy - 0.5), 1, 1, fc=divc(gnorm(np.clip(v, -2, 2))), ec="white", lw=0.3))
        recG.append(dict(gene=g_, dataset=d, log2FC=v, displayed_capped=np.clip(v, -2, 2) if not np.isnan(v) else np.nan))
axG.set_xlim(-0.5, 3.5); ylim_(axG)
axG.set_xticks([]); axG.set_yticks([])
for sp in axG.spines.values():
    sp.set_visible(False)
for j, d in enumerate(DS):
    axG.text(j, -0.62, NUM[d], fontsize=FS_S, ha="center", va="bottom", fontweight="bold")
caxG = axmm(cG["x"] + 3.4, cG["y"] + cG["h"] - 10.0, cG["w"] - 6.8, 1.6)
cbG = mpl.colorbar.ColorbarBase(caxG, cmap=divc, norm=gnorm, orientation="horizontal")
cbG.set_ticks([-2, 0, 2]); cbG.set_ticklabels(["≤−2", "0", "≥2"]); cbG.ax.tick_params(labelsize=FS_S, length=1.5, pad=1)
cbG.outline.set_linewidth(0.4)
caxG.text(0.5, -2.6, "log2FC", transform=caxG.transAxes, ha="center", va="top", fontsize=FS_S)
SRC["2G_per_dataset"] = pd.DataFrame(recG)
# H leave-one-dataset-out
axH = axmm(cH["x"] + 0.6, yt0, cH["w"] - 1.2, hrow)
hcols = ["full"] + [f"drop_{d}" for d in DS]
hlabs = ["All"] + [f"−{NUM[d]}" for d in DS]
recH = []
STABLE = ["CD163", "C1QA", "C1QB", "MERTK"]
for grp, g_, yy in rows:
    if g_ not in S6.index:
        continue
    for j, c_ in enumerate(hcols):
        k = KEY[S6.loc[g_, c_]]
        fc, ec = cls_style(k)
        axH.add_patch(Rectangle((j - 0.5 + 0.06, yy - 0.5 + 0.08), 0.88, 0.84, fc=fc,
                                ec=AMBIGUOUS_EDGE if k == "ambiguous" else "white", lw=0.5))
        recH.append(dict(gene=g_, run=c_, gc_class=S6.loc[g_, c_]))
    if g_ in STABLE:
        axH.add_patch(Rectangle((-0.5, yy - 0.5), 5, 1, fc="none", ec=TEXT, lw=1.1, zorder=5))
gc_rows = [r[2] for r in rows if r[0] == "GC\nscore"]
axH.text(2.0, np.mean(gc_rows), "excluded\n(define the\nscore)", fontsize=FS_S, ha="center", va="center", color="#555555", linespacing=1.0,
         fontstyle="italic")
assert all((S6.loc[g_, hcols] == "GC_induced").all() for g_ in STABLE)
assert all(S6.loc[g_, "n_changed"] == 0 for g_ in STABLE)
reg("H", "Bold outline: CD163, C1QA, C1QB, MERTK GC-induced in every run", "4 genes", "S6_leave_one_out",
    "rows with all 5 runs == GC_induced", "yes")
axH.set_xlim(-0.6, 4.6); ylim_(axH); axH.set_xticks([]); axH.set_yticks([])
for sp in axH.spines.values():
    sp.set_visible(False)
for j, s in enumerate(hlabs):
    axH.text(j, -0.62, s, fontsize=FS_S, ha="center", va="bottom", fontweight="bold")
fh = canvas(cH)
fh.add_patch(Rectangle((1.0, cH["h"] - 8.2), 4.0, 2.6, fc="none", ec=TEXT, lw=1.1))
t8(fh, 6.0, cH["h"] - 6.9, "GC-induced\nin all 5 runs", linespacing=1.0)
assert 6.0 + tw("GC-induced") < cH["w"]
SRC["2H_leave_one_out"] = pd.DataFrame(recH)

# ------------------------------------------------------------------ render-then-verify (geometry)
def layout_check(fig, W, H, tol=0.15):
    """sec 9.1: text-text overlaps and text outside the canvas (mm)."""
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    dpi = fig.dpi; k = 25.4 / dpi
    items = []
    for t in fig.findobj(mpl.text.Text):
        if not t.get_visible() or not t.get_text().strip():
            continue
        try:
            bb = t.get_window_extent(r)
        except Exception:
            continue
        items.append((t.get_text().replace("\n", "|")[:40], bb))
    out = [(s, round(bb.x0 * k, 1), round(bb.x1 * k, 1)) for s, bb in items
           if bb.x0 * k < -tol or bb.x1 * k > W + tol or bb.y0 * k < -tol or bb.y1 * k > H + tol]
    ov = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a_, b_ = items[i][1], items[j][1]
            ox = min(a_.x1, b_.x1) - max(a_.x0, b_.x0); oy = min(a_.y1, b_.y1) - max(a_.y0, b_.y0)
            if ox * k > tol and oy * k > tol:
                ov.append((items[i][0], items[j][0], round(ox * k, 1), round(oy * k, 1)))
    return out, ov


_out, _ov = layout_check(fig, W, H)
print("text outside canvas:", _out)
print("text-text overlaps (annotation boxes include leader lines):", _ov)

_sizes = sorted({round(t.get_fontsize(), 1) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()})
print("font sizes used (pt):", _sizes); assert min(_sizes) >= 8

_bad = min_font_audit(fig, 10); print("min_font_audit:", _bad); assert _bad == []

# ------------------------------------------------------------------ save
save_figure(fig, "Figure_2", tight=False)   # canvas drawn at exactly 183 mm
# JACC: RGB (no alpha) TIFF, LZW, 600 dpi
from PIL import Image as _Im
_t = _Im.open("Figure_2.tif"); _t.convert("RGB").save("Figure_2.tif", compression="tiff_lzw", dpi=(600, 600))

fig.savefig("Figure_2_preview_full.png", dpi=300)
with pd.ExcelWriter("Figure_2_source_data.xlsx") as xw:
    for k, v in SRC.items():
        v.to_excel(xw, sheet_name=k[:31], index=False)
_reg = pd.DataFrame(REG)
_reg["status_note"] = _reg.matches_locked_text
_reg["matches_locked_text"] = _reg.matches_locked_text.str.split(" ").str[0].str.strip()
assert set(_reg.matches_locked_text) <= {"yes", "no", "new"}, set(_reg.matches_locked_text)
_reg.to_csv("Figure_2_number_registry.csv", index=False)
print("saved Figure_2", fig.get_size_inches() / MM)
