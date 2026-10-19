#!/usr/bin/env python
"""Figure 1 (JACC: BTS v0.5). Study design: an in vitro glucocorticoid reference applied to
donor-controlled heart failure cohorts.  183 x 150 mm, 7 panels (A-G).

Standalone: reads only the input artifacts listed below (paths are the project artifact store);
writes Figure_1.{pdf,png,tif}, Figure_1_source_data.xlsx, Figure_1_number_registry.csv.
Every plotted / printed value is read from the files; schematic text follows the v0.4 manuscript.
"""
import os, json, tarfile
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Polygon, Circle, PathPatch
from matplotlib.path import Path
from matplotlib.lines import Line2D

import os  # v0.5 package: inputs resolved relative to this script (override with JBTS_INPUTS)
HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.environ.get("JBTS_INPUTS", os.path.join(HERE, "inputs"))
STYLE = os.path.join(HERE, "jbts_style.py")
TAR_FI = os.path.join(INP, "figure_inputs_v05.tar.gz")
TAR_PD = os.path.join(INP, "figure_inputs_v05_perdonor.tar.gz")
XL = os.path.join(INP, "Supplemental_Tables_S1-S14.xlsx")

exec(open(STYLE).read())
FAM = apply_jbts_style()
mpl.rcParams["savefig.pad_inches"] = 0          # canvas stays exactly 183 mm wide
mpl.rcParams["hatch.linewidth"] = 0.5

IN = "_inputs"
for tar, root in [(TAR_FI, "figure_inputs_v05"), (TAR_PD, "figure_inputs_v05_perdonor")]:
    if not os.path.isdir(f"{IN}/{root}"):
        os.makedirs(IN, exist_ok=True)
        with tarfile.open(tar) as t:
            t.extractall(IN)
FI = f"{IN}/figure_inputs_v05/"
PD = f"{IN}/figure_inputs_v05_perdonor/"

FS_T, FS_A, FS_S = SIZES["panel_title"], SIZES["tick"], SIZES["min_gene_label"]   # 10 / 9 / 8 pt
REG, SRC = [], {}


def reg(panel, text, value, source, computation, match):
    REG.append(dict(figure="Figure 1", panel=panel, printed_text=text, value=value,
                    source_file=source, computation=computation, matches_locked_text=match))


# ------------------------------------------------------------------ data
S2 = pd.read_excel(XL, "S2_datasets")
S3 = pd.read_excel(XL, "S3_gene_sets")
S5 = pd.read_excel(XL, "S5_GC_classification")
S1 = pd.read_excel(XL, "S1_deviations")

meta = pd.read_csv(FI + "02_data/processed/GSE290367_pseudobulk/pseudobulk_L1_meta.csv")
dup = pd.read_csv(FI + "04_results/validation_v1.1/donor_overlap/gse290367_duplicate_individuals.csv")
g = meta[(meta.region == "LV") & meta.disease.isin(["ND", "DCM"]) & (meta.age >= 20) & (meta.age_status == "postnatal")]
g = g[~g.donor_id.isin(dup.record_drop)].drop_duplicates("donor_id")
gse_ct = pd.crosstab(g.study, g.disease).reindex(columns=["ND", "DCM"], fill_value=0)
assert gse_ct.ND.sum() == 139 and gse_ct.DCM.sum() == 67, gse_ct.sum()

scp = pd.read_csv(PD + "v13_rerun/perdonor_scores_SCP3342_HFpEF.csv").drop_duplicates("unit_id")
assert len(scp) == 43
scp["pool_n"] = scp.pool.str.replace("pool", "").astype(int)

dec_g1 = json.load(open(FI + "04_results/validation_v1.2/gc_tests_G1/G1_decision.json"))
dec_v13 = json.load(open(FI + "04_results/validation_v1.3/v13_decisions.json"))
dec_g25 = json.load(open(PD + "v12_rerun/G2_G5_decisions.json"))

# ------------------------------------------------------------------ canvas / grid
W, H = WIDTH_MM, 150.0
GUT = 3.0
ROWS = [44.0, 45.0, 36.0, 16.0]                 # mm; sum + 3 gutters = 150
assert abs(sum(ROWS) + 3 * GUT - H) < 1e-9
fig = plt.figure(figsize=(W * MM, H * MM))
colw = (W - 11 * GUT) / 12
gs = fig.add_gridspec(4, 12, left=0, right=1, bottom=0, top=1, height_ratios=ROWS,
                      wspace=GUT / colw, hspace=GUT / np.mean(ROWS))
fig.patches.append(Rectangle((0, 0), 1, 1, transform=fig.transFigure, fc="white", ec="none", zorder=-10))


def cell(r, c0, c1):
    bb = gs[r, c0:c1].get_position(fig)
    return dict(x=bb.x0 * W, y=(1 - bb.y1) * H, w=bb.width * W, h=bb.height * H)


def axmm(x, y, w, h, **kw):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **kw)


def canvas(c, dy=0.0):
    """mm canvas over a cell; y increases downward."""
    ax = axmm(c["x"], c["y"] + dy, c["w"], c["h"] - dy)
    ax.set_xlim(0, c["w"]); ax.set_ylim(c["h"] - dy, 0); ax.axis("off")
    return ax


def header(c, letter, title, title_dx=5.6):
    fr = axmm(c["x"], c["y"], c["w"], c["h"]); fr.axis("off")
    panel_label(fr, letter, x=0, y=1 - 4.6 / c["h"])
    fr.text(title_dx / c["w"], 1 - 0.75 / c["h"], title, transform=fr.transAxes, fontsize=FS_T,
            fontweight="bold", va="top", ha="left", linespacing=1.05)
    return fr


def rbox(ax, x, y, w, h, fc="white", ec="#7F7F7F", lw=0.6, ls="-", r=1.0, z=1):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw,
                       ls=ls, zorder=z)
    ax.add_patch(p); return p


def arrow(ax, x0, y0, x1, y1, color="#555555", lw=0.9, ms=6, ls="-"):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=ms, ls=ls,
                                shrinkA=0, shrinkB=0), zorder=3)


def num(ax, x, y, k, r=1.55):
    ax.add_patch(Circle((x, y), r, fc=TEXT, ec="none", zorder=4))
    ax.text(x, y + 0.05, str(k), color="white", fontsize=FS_S, fontweight="bold", ha="center",
            va="center", zorder=5)


def heart_path(cx, cy, s):
    t = np.linspace(0, 2 * np.pi, 160)
    x = 16 * np.sin(t) ** 3
    y = 13 * np.cos(t) - 5 * np.cos(2 * t) - 2 * np.cos(3 * t) - np.cos(4 * t)
    x = cx + x / 34 * s; y = cy - (y + 2.5) / 34 * s      # y down
    return Path(np.c_[x, y], closed=True)


def heart(ax, cx, cy, s, fc, ec=TEXT, stripes=None):
    p = heart_path(cx, cy, s)
    if stripes:
        pp = PathPatch(p, fc="none", ec=ec, lw=0.7, zorder=3); ax.add_patch(pp)
        x0 = cx - s / 2; dw = s / len(stripes)
        for i, col in enumerate(stripes):
            rr = Rectangle((x0 + i * dw, cy - s / 2 - 0.5), dw, s + 1, fc=col, ec="none", zorder=2)
            ax.add_patch(rr); rr.set_clip_path(pp)
    else:
        ax.add_patch(PathPatch(p, fc=fc, ec=ec, lw=0.7, zorder=3))


def mk(ax, x, y, cohort, s=3.4, z=4):
    ax.plot([x], [y], marker=COHORT_MARKER[cohort], ms=s, mfc=COHORT[cohort], mec=TEXT, mew=0.4,
            ls="none", zorder=z)


def t8(ax, x, y, s, **kw):
    kw.setdefault("fontsize", FS_S); kw.setdefault("va", "center"); kw.setdefault("ha", "left")
    kw.setdefault("linespacing", 1.1)
    return ax.text(x, y, s, **kw)


from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties


def tw(s, size=FS_S, weight="normal", style="normal"):
    """Rendered width of a string in mm."""
    fp = FontProperties(family=mpl.rcParams["font.sans-serif"][0], size=size, weight=weight, style=style)
    return TextPath((0, 0), s, prop=fp).get_extents().width / 72 * 25.4


# ======================================================================= A  design schematic
cA = cell(0, 0, 12)
header(cA, "A", "Study design and analysis logic")
a = canvas(cA)
TOP0, TOP1, BOT0, BOT1 = 6.4, 25.6, 27.2, 43.4
LN = 3.1                                       # 8-pt line pitch, mm
a.add_patch(FancyBboxPatch((0.0, TOP0 - 0.9), cA["w"], BOT1 - TOP0 + 1.5,
                           boxstyle="round,pad=0,rounding_size=1.6", fc="#F2F2F2", ec="none", zorder=0))


def chip(x, y, s, ha="left"):
    """Figure-reference chip (figure and panels that show this step)."""
    w = tw(s, FS_S, "bold") + 1.6
    x0 = x if ha == "left" else x - w
    a.add_patch(FancyBboxPatch((x0, y - 1.65), w, 3.3, boxstyle="round,pad=0,rounding_size=0.8",
                               fc="#555555", ec="none", zorder=6))
    a.text(x0 + w / 2, y + 0.05, s, fontsize=FS_S, color="white", fontweight="bold", ha="center",
           va="center", zorder=7)
    return x0, x0 + w


def fits(x0, s, right, **kw):
    assert x0 + tw(s, FS_S, **kw) <= right + 1e-6, (s, round(x0 + tw(s, FS_S, **kw) - right, 2))


# --- context (not a numbered step): donor controls
X1, W1 = 0.8, 35.0
a.add_patch(FancyBboxPatch((X1, TOP0), W1, TOP1 - TOP0, boxstyle="round,pad=0,rounding_size=1.0",
                           fc="#E3E3E3", ec="none", zorder=1))
t8(a, X1 + 1.4, TOP0 + 2.35, "Donor controls", fontweight="bold")
t8(a, X1 + 1.4, TOP0 + 5.35, "(background)", fontstyle="italic", color="#444444")
bg_lines = ["ICU management often", "includes corticosteroids", "(refs 10–12); exposure", "not recorded in any cohort"]
for i, s in enumerate(bg_lines):
    t8(a, X1 + 1.4, TOP0 + 8.6 + i * 3.0, s)
    fits(X1 + 1.4, s, X1 + W1 - 0.6)
assert TOP0 + 8.6 + 3 * 3.0 + 1.5 < TOP1
# --- step 1: in vitro reference (bottom lane)
a.add_patch(FancyBboxPatch((X1, BOT0), W1, BOT1 - BOT0, boxstyle="round,pad=0,rounding_size=1.0",
                           fc="white", ec="#7F7F7F", lw=0.6, zorder=1))
num(a, X1 + 2.3, BOT0 + 2.3, 1)
t8(a, X1 + 4.4, BOT0 + 2.35, "4 in vitro datasets", fontweight="bold")
for i in range(4):
    cx = X1 + 3.3 + i * 4.5; cy = BOT0 + 7.4
    a.add_patch(Circle((cx, cy), 1.8, fc="white", ec=TEXT, lw=0.6, zorder=3))
    for dx, dy in [(-0.6, -0.4), (0.55, -0.55), (0.0, 0.6)]:
        a.add_patch(Circle((cx + dx, cy + dy), 0.42, fc=GC_CLASS["induced"], ec="none", zorder=4))
iv_lab = "GC vs vehicle"
fits(X1 + 1.4, iv_lab, X1 + W1 - 0.6)
t8(a, X1 + 1.4, BOT1 - 3.2, iv_lab)
# --- hearts: donor control vs failing heart
hx_c, hx_f, hy = 40.6, 50.4, TOP0 + 6.2
heart(a, hx_c, hy, 6.6, CONTROL)
heart(a, hx_f, hy, 6.6, None, stripes=[COHORT[k] for k in ["SCP3342", "Hahn", "GSE290367", "MAGNet"]])
a.plot([X1 + W1, hx_c - 3.6], [hy, hy], ls=(0, (1, 1.2)), color="#7F7F7F", lw=0.8, zorder=2)  # context, no arrow
t8(a, (hx_c + hx_f) / 2, hy + 0.2, "vs", ha="center", color="#555555")
t8(a, hx_c, TOP1 - 4.4, "Donor\ncontrol", ha="center", linespacing=1.0)
t8(a, hx_f, TOP1 - 4.4, "Failing\nheart", ha="center", linespacing=1.0)
assert hx_c - tw("control") / 2 > X1 + W1 + 0.4
X2, W2 = 56.6, 38.6
arrow(a, hx_f + 3.6, hy, X2, hy)
fx = (X1 + W1 + X2) / 2; fy = BOT0 + 1.6
a.add_patch(Polygon([[fx - 3.0, fy], [fx + 3.0, fy], [fx + 0.8, fy + 3.2], [fx + 0.8, fy + 5.0], [fx - 0.8, fy + 5.0],
                     [fx - 0.8, fy + 3.2]], closed=True, fc="#DADADA", ec=TEXT, lw=0.6, zorder=3))
t8(a, fx, BOT1 - 4.6, "Random-effects\nmeta-analysis", ha="center", linespacing=1.0)
assert fx - tw("Random-effects") / 2 > X1 + W1 + 0.3 and fx + tw("Random-effects") / 2 < X2 - 0.3
arrow(a, X1 + W1, BOT0 + 4.1, fx - 2.7, BOT0 + 4.1)
arrow(a, fx + 2.7, BOT0 + 4.1, X2, BOT0 + 4.1)
# --- cohorts (top) and classification (bottom)
rbox(a, X2, TOP0, W2, TOP1 - TOP0)
t8(a, X2 + 1.4, TOP0 + 2.35, "Transcriptomes", fontweight="bold")
ca, cb = X2 + 1.4, X2 + 20.4
t8(a, ca, TOP0 + 5.8, "HFpEF", color="#555555"); t8(a, cb, TOP0 + 5.8, "DCM", color="#555555")
for j, (coh, name, xx) in enumerate([("SCP3342", "SCP3342", ca), ("Hahn", "Hahn biopsy", ca),
                                     ("GSE290367", "GSE290367", cb), ("MAGNet", "MAGNet", cb)]):
    yy = TOP0 + 9.0 + (j % 2) * LN
    mk(a, xx + 1.0, yy, coh); t8(a, xx + 2.5, yy, name)
    fits(xx + 2.5, name, (cb - 0.6) if xx == ca else (X2 + W2 - 0.3))
yy = TOP0 + 9.0 + 2 * LN + 1.1
a.plot([X2 + 1.2, X2 + W2 - 1.2], [yy - 1.95, yy - 1.95], color="#BDBDBD", lw=0.5, zorder=2)
mk(a, ca + 1.0, yy, "GTEx"); t8(a, ca + 2.5, yy, "GTEx v11 (ventilator proxy)")
fits(ca + 2.5, "GTEx v11 (ventilator proxy)", X2 + W2 - 0.4)
rbox(a, X2, BOT0, W2, BOT1 - BOT0)
t8(a, X2 + 1.4, BOT0 + 2.35, "13,368 genes", fontweight="bold")
c6 = chip(X2 + W2 - 1.0, BOT0 + 2.5, "Fig 6", ha="right")
c2 = chip(c6[0] - 0.8, BOT0 + 2.5, "Fig 2", ha="right")
assert X2 + 1.4 + tw("13,368 genes", FS_S, "bold") < c2[0] - 0.5
reg("A", "13,368 genes", 13368, "S5_GC_classification", "row count", "yes")
for i, (k, lab) in enumerate([("induced", "Induced"), ("repressed", "Repressed"),
                              ("insensitive", "Insensitive"), ("ambiguous", "Ambiguous")]):
    xx = X2 + 1.4 + (i % 2) * 17.6; yy = BOT0 + 7.2 + (i // 2) * 3.8
    a.add_patch(Rectangle((xx, yy - 1.1), 2.2, 2.2, fc=GC_CLASS[k],
                          ec=AMBIGUOUS_EDGE if k == "ambiguous" else "none", lw=0.6, zorder=3))
    t8(a, xx + 2.9, yy, lab)
    fits(xx + 2.9, lab, X2 + W2 - 0.6)
reg("A", "4 in vitro datasets", 4, "S2_datasets", "rows with role == 'in vitro GC reference'", "yes")
# --- step 2: scores (spans both lanes)
X3, W3 = 97.0, 25.4
arrow(a, X2 + W2, hy, X3, hy + 2.0)
arrow(a, X2 + W2, BOT0 + 4.1, X3, BOT0 + 2.1)
rbox(a, X3, TOP0, W3, BOT1 - TOP0)
num(a, X3 + 2.3, TOP0 + 2.3, 2)
t8(a, X3 + 4.4, TOP0 + 2.35, "Scores", fontweight="bold")
sc_rows = [("GC score", dict(fontweight="bold"), 6.4), ("FKBP5, TSC22D3,", dict(fontstyle="italic"), 9.5),
           ("ZBTB16", dict(fontstyle="italic"), 12.6), ("8-gene GC", {}, 17.0), ("response", {}, 20.1),
           ("Circadian, cAMP", {}, 24.5), ("Macrophage", {}, 28.9), ("modules", {}, 32.0)]
for s, kw, dy in sc_rows:
    t8(a, X3 + 1.4, TOP0 + dy, s, **kw)
    fits(X3 + 1.4, s, X3 + W3 - 0.6, **{k.replace("font", ""): v for k, v in kw.items()})
c7 = chip(X3 + 1.4 + tw("modules") + 1.2, TOP0 + 32.0, "Fig 7")
assert c7[1] < X3 + W3 - 0.4
# --- step 3: contrast
X4, W4 = 124.2, 17.4
arrow(a, X3 + W3, (TOP0 + BOT1) / 2, X4, (TOP0 + BOT1) / 2)
rbox(a, X4, TOP0, W4, BOT1 - TOP0)
num(a, X4 + 2.3, TOP0 + 2.3, 3)
t8(a, X4 + 4.4, TOP0 + 2.35, "Contrast", fontweight="bold")
fits(X4 + 4.4, "Contrast", X4 + W4 - 0.4, weight="bold")
for i, s in enumerate(["Control", "minus case,", "per cohort", "and cell type"]):
    t8(a, X4 + 1.2, TOP0 + 6.2 + i * LN, s); fits(X4 + 1.2, s, X4 + W4 - 0.4)
chip(X4 + 1.2, TOP0 + 6.2 + 4 * LN + 0.5, "Fig 3")
for i, s in enumerate(["Adjusted for", "circadian or", "cAMP score"]):
    t8(a, X4 + 1.2, TOP0 + 23.2 + i * LN, s); fits(X4 + 1.2, s, X4 + W4 - 0.4)
chip(X4 + 1.2, TOP0 + 23.2 + 3 * LN + 0.5, "Fig 4")
assert TOP0 + 23.2 + 3 * LN + 0.5 + 1.65 < BOT1
# --- step 4: identifiability gate (post hoc, dashed) and outcomes
XG0 = 143.4; gw, gh = 9.4, 8.5
gx, gy = XG0 + gw, (TOP0 + TOP1) / 2
arrow(a, X4 + W4, gy, XG0, gy)
a.add_patch(Polygon([[gx - gw, gy], [gx, gy - gh], [gx + gw, gy], [gx, gy + gh]], closed=True, fc="white",
                    ec=TEXT, lw=0.8, ls=(0, (2.2, 1.2)), zorder=3))      # dashed: gate added post hoc
num(a, gx, gy - 4.3, 4, r=1.3)
t8(a, gx, gy + 0.6, "Identifiable?", ha="center", fontweight="bold", zorder=5)
XY = gx + gw + 1.6; WY = cA["w"] - 0.2 - XY
arrow(a, gx + gw, gy, XY, gy)
rbox(a, XY, TOP0, WY, TOP1 - TOP0, ec=TEXT)
t8(a, XY + 1.1, TOP0 + 3.0, "Yes:", fontweight="bold")
t8(a, XY + 1.3 + tw("Yes:", FS_S, "bold") + 0.6, TOP0 + 3.0, "adjust")
for i, s in enumerate(["for GC score;", "audit genes"]):
    t8(a, XY + 1.1, TOP0 + 3.0 + (i + 1) * LN, s); fits(XY + 1.1, s, XY + WY - 0.4)
cy5 = chip(XY + 1.1, TOP0 + 3.0 + 3 * LN + 1.0, "Fig 5E–I")
assert cy5[1] < XY + WY - 0.3
arrow(a, gx, gy + gh, gx, BOT0)
rbox(a, XG0, BOT0, cA["w"] - 0.2 - XG0, BOT1 - BOT0, ec=TEXT, ls=(0, (2.2, 1.2)))
t8(a, XG0 + 1.2, BOT0 + 3.0, "No:", fontweight="bold")
t8(a, XG0 + 1.4 + tw("No:", FS_S, "bold") + 0.6, BOT0 + 3.0, "report as non-identifiable")
fits(XG0 + 2.0 + tw("No:", FS_S, "bold"), "report as non-identifiable", cA["w"] - 0.6)
t8(a, XG0 + 1.2, BOT0 + 3.0 + LN, "AUC ≥0.9 or ≤0.1, or VIF ≥5")
t8(a, XG0 + 1.2, BOT0 + 3.0 + 2 * LN, "Illustrative; post hoc", color="#555555", fontstyle="italic")
cn5 = chip(cA["w"] - 1.2, BOT0 + 3.0 + 3 * LN + 0.6, "Fig 5A–D", ha="right")
assert BOT0 + 3.0 + 3 * LN + 0.6 + 1.65 < BOT1
for s_, v_ in [("AUC ≥0.9 or ≤0.1", "0.9; 0.1"), ("VIF ≥5", 5)]:
    reg("A", s_, v_, "Manuscript v0.4 Methods (Identifiability)", "threshold as stated", "yes")
reg("A", "refs 10–12", "10-12", "Manuscript v0.4 Introduction", "reference numbers as cited", "yes")
SRC["1A_schematic"] = pd.DataFrame({
    "element": ["context: donor controls", "in vitro datasets", "genes classified", "identifiability AUC threshold",
                "identifiability VIF threshold", "GC score genes", "step-to-figure map"],
    "value": ["ICU management often includes corticosteroids (refs 10-12); exposure not recorded in any cohort",
              4, 13368, "AUC >=0.9 or <=0.1", "VIF >=5", "FKBP5, TSC22D3, ZBTB16",
              "genes classified: Fig 2 (reference), Fig 6 (literature genes); macrophage modules: Fig 7; contrast: "
              "Fig 3; circadian/cAMP adjustment: Fig 4; non-identifiable: Fig 5A-D; adjust and audit: Fig 5E-I"],
    "source": ["Manuscript v0.4 Introduction", "S2_datasets", "S5_GC_classification (n rows)", "Manuscript v0.4 Methods",
               "Manuscript v0.4 Methods", "S3_gene_sets GC_specific3", "storyboard v0.5 + lead instructions"]})

# ======================================================================= B  cohort table
cB = cell(1, 0, 7)
header(cB, "B", "Cohorts and their design constraints")
b = canvas(cB)
CX = dict(name=4.0, tissue=22.6, ctrl=54.0, case=69.4, constr=71.0)
hy0 = 7.4
for k, lab, ha in [("name", "Cohort", "left"), ("tissue", "Tissue / assay", "left"), ("ctrl", "Controls", "right"),
                   ("case", "Cases", "right"), ("constr", "Design constraint", "left")]:
    t8(b, CX[k], hy0, lab, fontweight="bold", ha=ha, color="#444444")
b.plot([0.4, cB["w"]], [hy0 + 1.9, hy0 + 1.9], color=TEXT, lw=0.6)
gtex_n = {d["dataset"]: d["n"] for d in dec_g25["G5"]["rows"]}     # other (non-ventilator) deaths
trows = [
    ("SCP3342", "SCP3342", "Septum\nsnRNA-seq", "24", "19 HFpEF", "Pools of 5–6 donors;\nboth groups in every pool"),
    ("Hahn", "Hahn biopsy", "RV septum\nbulk RNA-seq", "24", "41 HFpEF", "Same program; patients\noverlap SCP3342"),
    ("GSE290367", "GSE290367", "LV\nsnRNA-seq", "139", "67 DCM", "Atlas of 10 studies;\ncases from 3 (panel C)"),
    ("MAGNet", "MAGNet", "LV\nbulk RNA-seq", "102", "128 DCM", "Samples overlapping\nGSE290367 excluded"),
    ("GTEx", "GTEx v11", "LV, AA\nbulk RNA-seq", f"LV {450 - gtex_n['GTEx_LV']}\nAA {461 - gtex_n['GTEx_AA']}",
     f"LV {gtex_n['GTEx_LV']}\nAA {gtex_n['GTEx_AA']}", "Ventilator-managed deaths\nas controls; not blinded"),
    (None, "In vitro ×4", "Monocyte, MDM;\narray, RNA-seq", "vehicle", "GC", "6, 4, 4 and 3 donors\n(datasets in Fig 2A)"),
]
assert 450 - gtex_n["GTEx_LV"] == 253 and 461 - gtex_n["GTEx_AA"] == 215
row_h = (cB["h"] - (hy0 + 2.4)) / len(trows)
rec = []
for i, (coh, name, tissue, ctrl, case, constr) in enumerate(trows):
    yy = hy0 + 2.4 + (i + 0.5) * row_h
    if i % 2 == 1:
        b.add_patch(Rectangle((0.4, yy - row_h / 2), cB["w"] - 0.4, row_h, fc="#F4F4F4", ec="none", zorder=0))
    if coh:
        mk(b, CX["name"] - 1.5, yy, coh, s=4.0)
    else:
        b.add_patch(Circle((CX["name"] - 1.5, yy), 1.1, fc="white", ec=TEXT, lw=0.6, zorder=3))
        b.add_patch(Circle((CX["name"] - 1.5, yy), 0.42, fc=GC_CLASS["induced"], ec="none", zorder=4))
    t8(b, CX["name"] + 0.3, yy, name, fontweight="bold")
    for k, s in [("tissue", tissue), ("constr", constr)]:
        t8(b, CX[k], yy, s, linespacing=1.0)
        assert CX[k] + max(tw(z) for z in s.split("\n")) < (CX["ctrl"] - 7 if k == "tissue" else cB["w"]), s
    for k, s in [("ctrl", ctrl), ("case", case)]:
        t8(b, CX[k], yy, s, ha="right", linespacing=1.0)
    if coh == "GTEx":     # controls = ventilator-managed deaths (donor-like); cases = other deaths
        ctrl, case = f"ventilator-managed: {ctrl}", f"other deaths: {case}"
    rec.append(dict(cohort=name, tissue_assay=tissue.replace("\n", ", "), controls=ctrl.replace("\n", "; "),
                    cases=case.replace("\n", "; "), design_constraint=constr.replace("\n", " ")))
y_a = hy0 + 2.4 + 0.1 * row_h; y_b = hy0 + 2.4 + 1.9 * row_h
b.plot([0.9, 0.45, 0.45, 0.9], [y_a, y_a, y_b, y_b], color=COHORT["SCP3342"], lw=0.9, solid_capstyle="butt")
for txt, val, src, comp, m in [
        ("SCP3342 24 / 19 HFpEF", "24; 19", "S2_datasets", "as tabulated", "yes"),
        ("SCP3342 pools of 5–6 donors", "5-6", "S2_datasets", "as tabulated", "yes"),
        ("Hahn biopsy 24 / 41 HFpEF", "24; 41", "S2_datasets", "as tabulated", "yes"),
        ("GSE290367 139 / 67 DCM", "139; 67", "pseudobulk_L1_meta.csv + duplicate_individuals.csv",
         "LV, ND/DCM, age>=20, postnatal, drop duplicate records, unique donors", "yes"),
        ("GSE290367 atlas of 10 studies; cases from 3", "10; 3", "pseudobulk_L1_meta.csv",
         "unique studies; studies with DCM>0", "yes"),
        ("MAGNet 102 / 128 DCM", "102; 128", "S2_datasets", "as tabulated", "yes"),
        ("GTEx controls (ventilator-managed) LV 253, AA 215; cases (other deaths) LV 197, AA 246", "197; 253", "S2_datasets (450, 253); G2_G5_decisions.json G5 n=197",
         "other = total - ventilator", "new (S: 450, 253, 461, 215 and 197 are in the text; AA other 246 is from S7b G5 n)"),
        ("In vitro 6, 4, 4 and 3 donors", "6; 4; 4; 3", "S2_datasets", "as tabulated", "new (S)")]:
    reg("B", txt, val, src, comp, m)
assert gse_ct.index.nunique() == 10 and (gse_ct.DCM > 0).sum() == 3
SRC["1B_cohorts"] = pd.DataFrame(rec)

# ======================================================================= C  GSE290367 study composition
cC = cell(1, 7, 12)
header(cC, "C", "GSE290367: DCM from 3 of 10 studies")
both = gse_ct[(gse_ct.ND > 0) & (gse_ct.DCM > 0)].sort_values("DCM", ascending=False)
only = gse_ct[gse_ct.DCM == 0].sort_values("ND", ascending=False)
order = list(both.index) + list(only.index)
short = {s: ("ENCODE v4" if s.startswith("ENCODE") else s.split(" ")[0]) for s in order}
axC = axmm(cC["x"] + 15.6, cC["y"] + 9.4, cC["w"] - 15.6 - 1.6, cC["h"] - 9.4 - 7.4)
yy = np.arange(len(order))
nd = gse_ct.loc[order, "ND"].values; dc = gse_ct.loc[order, "DCM"].values
axC.axhspan(-0.5, 2.5, color="#E3F1EF", lw=0, zorder=0)
axC.barh(yy, nd, color=CONTROL, height=0.7, lw=0, zorder=2)
axC.barh(yy, dc, left=nd, color=COHORT["GSE290367"], height=0.7, lw=0, zorder=2)
for i in range(len(order)):
    lab = f"{nd[i]} + {dc[i]}" if dc[i] else f"{nd[i]}"
    axC.text(nd[i] + dc[i] + 1.0, i, lab, fontsize=FS_S, va="center", ha="left", zorder=3)
    reg("C", f"{short[order[i]]} {lab}", f"{nd[i]};{dc[i]}", "pseudobulk_L1_meta.csv + duplicate_individuals.csv",
        "unique adult LV donors per study after de-duplication (ND;DCM)", "new (N)")
axC.set_yticks(yy); axC.set_yticklabels([short[s] for s in order], fontsize=FS_S)
axC.set_ylim(len(order) - 0.5, -0.5)
axC.set_xlim(0, 66); axC.set_xticks([0, 20, 40, 60])
axC.set_xlabel("Donors (n)", labelpad=1.0)
axC.tick_params(axis="y", length=0, pad=1.5)
axC.text(64.5, 6.6, "Shaded: studies with\nboth groups, used for\nDCM contrasts:\n45 vs 67",
         fontsize=FS_S, ha="right", va="center", linespacing=1.05, color="#1F6F66")
assert both.ND.sum() == 45 and both.DCM.sum() == 67
reg("C", "DCM contrasts: 45 vs 67", "45;67", "pseudobulk_L1_meta.csv", "sum over studies with both groups",
    "new (S: S8a n/n_case)")
kc = canvas(cC)
for j, (col, lab) in enumerate([(CONTROL, "Non-failing (139)"), (COHORT["GSE290367"], "DCM (67)")]):
    xx = 6.0 + j * 30.0
    kc.add_patch(Rectangle((xx, 5.9), 2.4, 2.2, fc=col, ec="none"))
    t8(kc, xx + 3.1, 7.0, lab)
reg("C", "Non-failing (139); DCM (67)", "139;67", "pseudobulk_L1_meta.csv", "column totals", "yes")
SRC["1C_GSE290367_studies"] = pd.DataFrame({"study": order, "nonfailing_n": nd, "DCM_n": dc,
                                            "used_for_DCM_contrasts": [s in both.index for s in order]})

# ======================================================================= D  SCP3342 pools
cD = cell(2, 0, 3)
header(cD, "D", "SCP3342: both groups\nin every pool")
axD = axmm(cD["x"] + 14.0, cD["y"] + 11.4, cD["w"] - 14.4, cD["h"] - 11.4 - 8.2)
pools = sorted(scp.pool_n.unique())
recD = []
for j, p in enumerate(pools):
    for r_, (cv, col) in enumerate([(0, CONTROL), (1, COHORT["SCP3342"])]):
        sub = scp[(scp.pool_n == p) & (scp.case == cv)].sort_values("unit_id")
        k = len(sub)
        for q in range(k):
            axD.plot(j, r_ + (q - (k - 1) / 2) * 0.28, "o", ms=3.6, mfc=col, mec=TEXT, mew=0.35)
        recD += [dict(pool=f"pool{p}", group="HFpEF" if cv else "control", unit_id=u) for u in sub.unit_id]
axD.set_xlim(-0.6, len(pools) - 0.4); axD.set_ylim(1.55, -0.55)
axD.set_xticks(range(len(pools))); axD.set_xticklabels([str(p) for p in pools])
axD.set_yticks([0, 1]); axD.set_yticklabels(["Control\n(n=24)", "HFpEF\n(n=19)"], linespacing=1.0)
axD.tick_params(axis="y", length=0, pad=1.0); axD.tick_params(axis="x", length=0, pad=1.5)
for sp in ["left", "bottom"]:
    axD.spines[sp].set_visible(False)
for j in range(len(pools) - 1):
    axD.axvline(j + 0.5, color=GRID, lw=0.5, zorder=0)
axD.set_xlabel("Pool", labelpad=1.0)
cnt = pd.crosstab(scp.pool_n, scp.case)
assert (cnt[0.0] == 3).all() and cnt[1.0].between(2, 3).all() and len(pools) == 8
reg("D", "Control (n=24); HFpEF (n=19)", "24;19", "perdonor_scores_SCP3342_HFpEF.csv", "unique unit_id by case", "yes")
reg("D", "Pools 2, 3, 4, 5, 7, 8, 9, 10 (8 pools)", 8, "perdonor_scores_SCP3342_HFpEF.csv", "unique pool", "yes")
SRC["1D_SCP3342_pools"] = pd.DataFrame(recD)

# ======================================================================= E  score genes x in vitro class
cE = cell(2, 3, 9)
header(cE, "E", "Score genes and their in vitro GC class")
e = canvas(cE)
cls = S5.set_index("gene")["class"]
cmap_cls = {"GC_induced": "induced", "GC_repressed": "repressed", "GC_insensitive": "insensitive",
            "ambiguous": "ambiguous"}
gs3 = S3.groupby("gene_set").gene.apply(list)
gc3 = gs3["GC_specific3"]; gcr = [x for x in gs3["GC_response"] if x not in gc3]
erows = [("GC score", "v1.3", gc3), ("8-gene GC\nresponse", "v1.2", gcr),
         ("Morning +", None, gs3["CLOCK_morning"]), ("Evening −", None, gs3["CLOCK_evening"]),
         ("cAMP", None, gs3["CAMP_IEG"])]
ys = [7.9, 12.6, 20.4, 25.1, 29.8]
xg0 = 17.0
recE = []


def tile(x, y, gname, k):
    wbox = tw(gname, FS_S, style="italic") + 1.2
    rr = Rectangle((x, y - 1.75), wbox, 3.5, fc=GC_CLASS[k], ec=(AMBIGUOUS_EDGE if k == "ambiguous" else
                   ("#9A9A9A" if k == "not_measured" else "none")), lw=0.6,
                   hatch="////" if k == "not_measured" else None, zorder=2)
    e.add_patch(rr)
    e.text(x + wbox / 2, y + 0.05, gname, fontsize=FS_S, ha="center", va="center", fontstyle="italic", zorder=3)
    return wbox


for i, (lab, tag, genes) in enumerate(erows):
    yy = ys[i]
    t8(e, 0.6, yy, lab, linespacing=0.95)
    xx = xg0
    if i == 1:
        t8(e, xx, yy, "GC score +", color="#555555"); xx += tw("GC score +", FS_S) + 1.2
    for gname in genes:
        c_ = cls.get(gname, None); k = cmap_cls.get(c_, "not_measured")
        wb = tile(xx, yy, gname, k)
        recE.append(dict(score=lab.replace("\n", " "), plan=tag or "", gene=gname,
                         in_vitro_class=c_ if c_ else "not measured in >=2 in vitro datasets (not in S5)"))
        if i == 1 and gname == "PER1": x_beg_ex = xx
        if i == 1 and gname == "DDIT4": x_end_ex = xx + wb
        xx += wb + 0.45
    if tag:
        t8(e, xx + 0.6, yy, tag, color="#555555", fontstyle="italic")
    assert xx < cE["w"], (lab, xx)
yb = ys[1] + 2.35
e.plot([x_beg_ex, x_beg_ex, x_end_ex, x_end_ex], [yb - 0.4, yb, yb, yb - 0.4], color=TEXT, lw=0.6)
t8(e, x_beg_ex, yb + 1.55, "not in GC score: also clock- or cAMP-responsive", color="#444444")
ky = cE["h"] - 1.6; xx = 0.6
assert not any(r_["in_vitro_class"] == "GC_repressed" for r_ in recE)
for k, lab in [("induced", "Induced"), ("insensitive", "Insensitive"), ("ambiguous", "Ambiguous"),
               ("not_measured", "Not measured in ≥2 datasets")]:
    rr = Rectangle((xx, ky - 1.1), 2.2, 2.2, fc=GC_CLASS[k], ec=AMBIGUOUS_EDGE if k == "ambiguous" else
                   ("#9A9A9A" if k == "not_measured" else "none"), lw=0.6, hatch="////" if k == "not_measured" else None)
    e.add_patch(rr); t8(e, xx + 2.8, ky, lab)
    xx += 2.8 + tw(lab, FS_S) + 1.6
    assert xx - 1.6 < cE["w"] - 0.3, xx
SRC["1E_score_genes"] = pd.DataFrame(recE)

# ======================================================================= F  analysis history
cF = cell(2, 9, 12)
header(cF, "F", "Analysis plans,\nfrozen internally")
f = canvas(cF)
nodes = [("v1.0", True, False, "Before any data"),
         ("v1.1", True, False, "Modules after discovery"),
         ("v1.2", True, True, "G1–G5"),
         ("v1.3", True, True, "C1–C3, K1, K2, A1–A4"),
         ("Post hoc", False, False, "Labeled in text")]
xt = 2.0; ytop, ybot = 11.2, cF["h"] - 0.9
f.plot([xt, xt], [ytop, ybot], color="#7F7F7F", lw=0.9, zorder=1)
for i, (ver, locked, sha, desc) in enumerate(nodes):
    yy = ytop + i * (ybot - ytop) / (len(nodes) - 1)
    if locked:
        f.add_patch(Rectangle((xt - 1.35, yy - 0.75), 2.7, 2.1, fc=TEXT, ec="none", zorder=3))
        th = np.linspace(np.pi, 2 * np.pi, 30)
        f.plot(xt + 0.85 * np.cos(th), yy - 0.75 + 1.2 * np.sin(th), color=TEXT, lw=0.9, zorder=3)
        f.add_patch(Rectangle((xt - 1.6, yy - 2.25), 3.2, 3.9, fc="#FFFFFF", ec="none", zorder=2))
    else:
        f.add_patch(Circle((xt, yy), 1.35, fc="white", ec=TEXT, lw=0.7, ls=(0, (1.4, 0.9)), zorder=3))
    xx = xt + 2.6
    f.text(xx, yy - 1.2, ver, fontsize=FS_S, fontweight="bold", va="center", ha="left")
    xx += tw(ver, FS_S, "bold") + 1.1
    if sha:
        f.text(xx + 0.3, yy - 1.2, "SHA-256", fontsize=FS_S, va="center", ha="left", color="white",
               fontweight="bold", bbox=dict(boxstyle="round,pad=0.12,rounding_size=0.4", fc="#7F7F7F", ec="none"))
    f.text(xt + 2.6, yy + 1.45, desc, fontsize=FS_S, va="center", ha="left", color="#333333")
    assert xt + 2.6 + tw(desc, FS_S) < cF["w"], ver
SRC["1F_history"] = pd.DataFrame([dict(plan=n[0], frozen=n[1], sha256_hash_stated=n[2], content=n[3])
                                  for n in nodes])


# ======================================================================= G  scoreboard
cG = cell(3, 0, 12)
header(cG, "G", "Outcomes of the pre-specified decision tests (Table 1 order)")
gg = canvas(cG)
tests = [  # id, outcome glyph, outcome text, where shown, value check
    ("v1.0", "notmet", "Not met", "Suppl. S-R1"),
    ("G1", "fail", "Failed", "Fig 3A"),
    ("G2", "grey", "No decision", "Fig 7E"),
    ("G3", "fail", "Failed", "Fig 7C"),
    ("G4a", "dash", "Not evaluated", "Fig 7D"),
    ("G4b", "notmet", "Not met", "Fig 7C"),
    ("G5", "pass", "Passed", "Fig 7H"),
    ("C2", "pass", "Passed", "Fig 4E, G"),
    ("K1", "pass", "Passed", "Fig 4E, G"),
    ("A3", "pass", "Passed", "Fig 5H"),
]
# outcome checks against the result files
assert dec_g1["G1_PASS"] is False
assert dec_g25["G3"]["PASS"] is False
assert dec_g25["G4"]["prediction_b_attenuation_ge50_both"] is False
assert dec_g25["G5"]["PASS_LV_rawTOST"] is True
assert dec_v13["C2"]["PASS"] and dec_v13["K1"]["PASS"] and dec_v13["A3"]["PASS"]
nT = len(tests); tw_ = (cG["w"] - 1.0) / nT; pw, ph = tw_ - 3.6, 4.4; py = 5.3
recG = []
for i, (tid, kind, out, where) in enumerate(tests):
    x0 = 0.6 + i * tw_ + (tw_ - pw) / 2
    style = dict(boxstyle="round,pad=0,rounding_size=2.4")
    if kind == "pass":
        gg.add_patch(FancyBboxPatch((x0, py), pw, ph, fc=TEXT, ec=TEXT, lw=0.8, **style)); tc = "white"
        lab = tid
    elif kind == "notmet":    # same weight as 'failed'; en dash instead of x
        gg.add_patch(FancyBboxPatch((x0, py), pw, ph, fc="white", ec=TEXT, lw=1.8, **style)); tc = TEXT
        lab = tid + "  –"
    elif kind == "fail":
        gg.add_patch(FancyBboxPatch((x0, py), pw, ph, fc="white", ec=TEXT, lw=1.8, **style)); tc = TEXT
        lab = tid + "  ×"
    elif kind == "dash":
        gg.add_patch(FancyBboxPatch((x0, py), pw, ph, fc="white", ec=TEXT, lw=0.8, ls=(0, (2.0, 1.2)), **style))
        tc = TEXT; lab = tid
    elif kind == "grey":
        gg.add_patch(FancyBboxPatch((x0, py), pw, ph, fc="#D9D9D9", ec="#7F7F7F", lw=0.8, **style)); tc = TEXT
        lab = tid
    else:   # half-filled
        outline = FancyBboxPatch((x0, py), pw, ph, fc="white", ec=TEXT, lw=0.8, **style)
        gg.add_patch(outline)
        half = Rectangle((x0, py), pw / 2, ph, fc=TEXT, ec="none"); gg.add_patch(half); half.set_clip_path(outline)
        tc = None; lab = tid
    if kind == "half":
        gg.text(x0 + pw / 4, py + ph / 2 + 0.1, lab, color="white", fontsize=FS_A, fontweight="bold",
                ha="center", va="center", zorder=5)
    else:
        gg.text(x0 + pw / 2, py + ph / 2 + 0.1, lab, color=tc, fontsize=FS_A, fontweight="bold",
                ha="center", va="center", zorder=5)
    gg.text(x0 + pw / 2, py + ph + 1.9, out, fontsize=FS_S, ha="center", va="center")
    gg.text(x0 + pw / 2, py + ph + 4.75, where, fontsize=FS_S, ha="center", va="center", color="#555555")
    recG.append(dict(test=tid, outcome=out, glyph=kind, shown_in=where))
    reg("G", f"{tid}: {out}", out, "G1_decision.json; G2_G5_decisions.json; v13_decisions.json; Table 1",
        "outcome as declared in Table 1", "yes")
SRC["1G_tests"] = pd.DataFrame(recG)

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

# ------------------------------------------------------------------ save
save_figure(fig, "Figure_1")
# JACC: RGB (no alpha) TIFF, LZW, 600 dpi
from PIL import Image as _Im
_t = _Im.open("Figure_1.tif"); _t.convert("RGB").save("Figure_1.tif", compression="tiff_lzw", dpi=(600, 600))

fig.savefig("Figure_1_preview_full.png", dpi=300)
with pd.ExcelWriter("Figure_1_source_data.xlsx") as xw:
    for k, v in SRC.items():
        v.to_excel(xw, sheet_name=k[:31], index=False)
_reg = pd.DataFrame(REG)
_reg["status_note"] = _reg.matches_locked_text
_reg["matches_locked_text"] = _reg.matches_locked_text.str.split(" ").str[0].str.strip()
assert set(_reg.matches_locked_text) <= {"yes", "no", "new"}, set(_reg.matches_locked_text)
_reg.to_csv("Figure_1_number_registry.csv", index=False)
print("saved Figure_1", fig.get_size_inches() / MM)
