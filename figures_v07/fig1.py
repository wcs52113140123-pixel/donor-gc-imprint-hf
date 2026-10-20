#!/usr/bin/env python
"""Figure 1 (JACC: BTS v0.7). Study design: an in vitro glucocorticoid reference applied to
heart failure cohorts with donor controls.  183 mm wide (height set in code, <=240 mm), 7 panels (A-G);
every text element >=10 pt.  v0.7b: panel A redesigned on one grid (FIG1_A_VARIANT=1|2).  v0.7: reader-facing content only (analysis-plan history and decision-test
outcome panels removed; plan tags removed from E); new F (donor age) and G (sex composition) from the
per-donor bundle.

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

import os  # public package: inputs resolved relative to this script (override with JBTS_INPUTS)
HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.environ.get("JBTS_INPUTS", os.path.join(HERE, "..", "figures_v05", "inputs"))
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

FS_T, FS_A, FS_S = SIZES["panel_title"], SIZES["tick"], SIZES["min_any"]   # 11 / 10 / 10 pt (v0.6: >=10 pt)
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

dec_g25 = json.load(open(PD + "v12_rerun/G2_G5_decisions.json"))

# ------------------------------------------------------------------ canvas (v0.6: every text >=10 pt)
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties

W = WIDTH_MM
GUT = 3.0
LN = 4.2                                   # 10-pt line pitch (mm)
A_H, BC_H, R3_H, R4_H = 78.6, 66.0, 43.0, 43.2
Y_A, Y_BC = 0.0, 78.6 + GUT
Y_R3 = Y_BC + BC_H + GUT
Y_R4 = Y_R3 + R3_H + GUT
H = Y_R4 + R4_H
assert H <= 240.0, H
fig = plt.figure(figsize=(W * MM, H * MM))
fig.patches.append(Rectangle((0, 0), 1, 1, transform=fig.transFigure, fc="white", ec="none", zorder=-10))


def box(x, y, w, h):
    return dict(x=x, y=y, w=w, h=h)


def axmm(x, y, w, h, **kw):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **kw)


def canvas(c):
    ax = axmm(c["x"], c["y"], c["w"], c["h"])
    ax.set_xlim(0, c["w"]); ax.set_ylim(c["h"], 0); ax.axis("off")
    return ax


def tw(s, size=FS_S, weight="normal", style="normal"):
    """Rendered width of a string in mm."""
    fp = FontProperties(family=mpl.rcParams["font.sans-serif"][0], size=size, weight=weight, style=style)
    return TextPath((0, 0), s, prop=fp).get_extents().width / 72 * 25.4


def wrap(text, width, size=FS_S, weight="normal", style="normal"):
    lines, cur = [], ""
    for w_ in text.split():
        trial = (cur + " " + w_).strip()
        if tw(trial, size, weight, style) <= width or not cur:
            cur = trial
        else:
            lines.append(cur); cur = w_
    lines.append(cur)
    assert all(tw(l, size, weight, style) <= width + 1e-6 for l in lines), (text, width)
    return lines


def header(c, letter, title, title_dx=5.8):
    """Panel letter (14 pt) + bold title (11 pt), wrapped to the panel width. Returns title height (mm)."""
    fr = axmm(c["x"], c["y"], c["w"], c["h"]); fr.axis("off")
    panel_label(fr, letter, x=0, y=1 - 5.0 / c["h"])
    lines = wrap(title, c["w"] - title_dx - 0.3, FS_T, "bold")
    fr.text(title_dx / c["w"], 1 - 0.75 / c["h"], "\n".join(lines), transform=fr.transAxes, fontsize=FS_T,
            fontweight="bold", va="top", ha="left", linespacing=1.05)
    return 0.75 + len(lines) * 4.6


def rbox(ax, x, y, w, h, fc="white", ec="#7F7F7F", lw=0.6, ls="-", r=1.0, z=1):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw,
                       ls=ls, zorder=z)
    ax.add_patch(p); return p


def arrow(ax, x0, y0, x1, y1, color="#555555", lw=0.9, ms=6):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                mutation_scale=ms, shrinkA=0, shrinkB=0), zorder=3)


def num(ax, x, y, k, r=1.95):
    ax.add_patch(Circle((x, y), r, fc=TEXT, ec="none", zorder=4))
    ax.text(x, y + 0.05, str(k), color="white", fontsize=FS_S, fontweight="bold", ha="center", va="center",
            zorder=5)


def heart_path(cx, cy, s):
    t = np.linspace(0, 2 * np.pi, 160)
    x = 16 * np.sin(t) ** 3
    y = 13 * np.cos(t) - 5 * np.cos(2 * t) - 2 * np.cos(3 * t) - np.cos(4 * t)
    return Path(np.c_[cx + x / 34 * s, cy - (y + 2.5) / 34 * s], closed=True)


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


def mk(ax, x, y, cohort, s=4.0, z=4):
    ax.plot([x], [y], marker=COHORT_MARKER[cohort], ms=s, mfc=COHORT[cohort], mec=TEXT, mew=0.4, ls="none",
            zorder=z)


def t8(ax, x, y, s, **kw):
    kw.setdefault("fontsize", FS_S); kw.setdefault("va", "center"); kw.setdefault("ha", "left")
    kw.setdefault("linespacing", 1.1)
    return ax.text(x, y, s, **kw)


def chip(ax, x, y, s, ha="left"):
    """Figure-reference chip (figure panels that show this step)."""
    w = tw(s, FS_S, "bold") + 1.8
    x0 = x if ha == "left" else x - w
    ax.add_patch(FancyBboxPatch((x0, y - 2.0), w, 4.0, boxstyle="round,pad=0,rounding_size=0.9", fc="#555555",
                                ec="none", zorder=6))
    ax.text(x0 + w / 2, y + 0.05, s, fontsize=FS_S, color="white", fontweight="bold", ha="center", va="center",
            zorder=7)
    return x0, x0 + w


def fits(x0, s, right, **kw):
    assert x0 + tw(s, FS_S, **kw) <= right + 1e-6, (s, round(x0 + tw(s, FS_S, **kw) - right, 2))


# ======================================================================= A  design schematic (v0.7b: one grid)
# Rows: context strip | inputs (in vitro, hearts) | 4 equal step cards | outcomes.  One card style (rounded,
# thin mid-gray border, white fill); step cards share a dark header band with a white number badge; chips
# (panels showing each step) bottom-right; uniform 1.2-pt solid arrows, orthogonal, no crossings.
A_VARIANT = int(os.environ.get("FIG1_A_VARIANT", "1"))
cA = box(0, Y_A, W, A_H)
header(cA, "A", "Study design and analysis logic")
a = canvas(cA)
BORDER, BAND, CTX = "#8C8C8C", "#3A3A3A", "#EDEDED"
PADX, PIT = 1.2, 4.2                                  # inner padding, 10-pt line pitch (mm)
CW, GAP = 42.75, 4.0                                  # step-card width, gap (4 cards + 3 gaps = 183 mm)
SX = [i * (CW + GAP) for i in range(4)]
assert abs(SX[-1] + CW - W) < 1e-6
Y_CTX, H_CTX = 6.2, 9.6
Y_IN, H_IN = Y_CTX + H_CTX + 2.2, 10.0
Y_ST, H_ST = Y_IN + H_IN + 3.8, 32.0
Y_OUT = Y_ST + H_ST + 4.4
H_OUT = 10.2
assert Y_OUT + H_OUT <= A_H - 0.15, Y_OUT + H_OUT
BAND_H = 5.4


def card(x, y, w, h, ls="-"):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2", fc="white", ec=BORDER, lw=0.6,
                       ls=ls, zorder=2)
    a.add_patch(p)
    return p


def arrowA(x0, y0, x1, y1):
    a.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color="#444444", lw=1.2,
               mutation_scale=7.5, shrinkA=0, shrinkB=0, joinstyle="miter"), zorder=6)


def polyarrow(pts):
    """Orthogonal arrow through pts (mm); head on the last segment."""
    xs_, ys_ = zip(*pts[:-1])
    a.plot(list(xs_) + [pts[-2][0]], list(ys_) + [pts[-2][1]], color="#444444", lw=1.2, solid_joinstyle="miter",
           solid_capstyle="butt", zorder=6)
    arrowA(*pts[-2], *pts[-1])


def step_card(i, title, ls="-", accent=None):
    x, y = SX[i], Y_ST
    p = card(x, y, CW, H_ST, ls=ls)
    band = Rectangle((x, y), CW, BAND_H, fc=BAND, ec="none", zorder=3)
    a.add_patch(band); band.set_clip_path(p)
    if accent:
        a.add_patch(Rectangle((x + 0.3, y + BAND_H), CW - 0.6, 0.7, fc=accent, ec="none", zorder=3))
    a.add_patch(Circle((x + PADX + 1.9, y + BAND_H / 2), 1.85, fc="white", ec="none", zorder=4))
    a.text(x + PADX + 1.9, y + BAND_H / 2 + 0.05, str(i + 1), color=BAND, fontsize=FS_S, fontweight="bold",
           ha="center", va="center", zorder=5)
    a.text(x + PADX + 4.9, y + BAND_H / 2 + 0.05, title, color="white", fontsize=FS_S, fontweight="bold",
           ha="left", va="center", zorder=5)
    fits(x + PADX + 4.9, title, x + CW - PADX, weight="bold")
    return x


def body_y(k):                                        # centre of body line k (0-based) in a step card
    return Y_ST + BAND_H + 3.0 + k * PIT


def body_lines(x, lines, k0=0, **kw):
    for j, s in enumerate(lines):
        t8(a, x + PADX, body_y(k0 + j), s, zorder=5, **kw)
        fits(x + PADX, s, x + CW - PADX, **{k_.replace("font", ""): v for k_, v in kw.items() if k_.startswith("font")})
    return k0 + len(lines)


CHIP_Y = body_y(5)                                   # 6th body line, bottom-right
assert CHIP_Y + 2.0 <= Y_ST + H_ST - 0.5


def chips(x_right, y, labels):
    xr_ = x_right - 1.0; spans = []
    for s in labels[::-1]:
        x0_, x1_ = chip(a, xr_, y, s, ha="right"); spans.append((x0_, x1_)); xr_ = x0_ - 0.7
    return min(s_[0] for s_ in spans)


# --- 1. context strip (full width, light gray)
a.add_patch(FancyBboxPatch((0.0, Y_CTX), W, H_CTX, boxstyle="round,pad=0,rounding_size=1.2", fc=CTX, ec="none",
                           zorder=1))
hy = Y_CTX + H_CTX / 2
hx_c, hx_f = 4.6, 15.8
heart(a, hx_c, hy, 6.4, CONTROL)
heart(a, hx_f, hy, 6.4, None, stripes=[COHORT[k] for k in ["SCP3342", "Hahn", "GSE290367", "MAGNet"]])
t8(a, (hx_c + hx_f) / 2, hy + 0.2, "vs", ha="center", color="#555555")
XT = hx_f + 4.6
lead = "Donor controls:"
ctx_text = "ICU management often includes corticosteroids (refs 10–12); exposure not recorded in any cohort"
wl = tw(lead, FS_S, "bold") + 1.8
words, l1, l2 = ctx_text.split(), "", ""
for w_ in words:
    if not l2 and tw((l1 + " " + w_).strip()) <= W - 3.5 - XT - wl:
        l1 = (l1 + " " + w_).strip()
    else:
        l2 = (l2 + " " + w_).strip()
assert tw(l2) <= W - 3.5 - XT, l2
t8(a, XT, hy - PIT / 2, lead, fontweight="bold", zorder=5)
t8(a, XT + wl, hy - PIT / 2, l1, zorder=5)
t8(a, XT, hy + PIT / 2, l2, zorder=5)

# --- 2. inputs: in vitro datasets | hearts
lead_iv = "4 in vitro datasets:"
iv_l2 = "human monocytes/macrophages"
WIV = max(tw(lead_iv, FS_S, "bold") + 1.2 + tw("GC vs vehicle"), tw(iv_l2)) + 2 * PADX + 0.4
card(0.0, Y_IN, WIV, H_IN)
y1_, y2_ = Y_IN + H_IN / 2 - PIT / 2, Y_IN + H_IN / 2 + PIT / 2
t8(a, PADX, y1_, lead_iv, fontweight="bold", zorder=5)
t8(a, PADX + tw(lead_iv, FS_S, "bold") + 1.2, y1_, "GC vs vehicle", zorder=5)
t8(a, PADX, y2_, iv_l2, zorder=5)
XH = WIV + GAP; WH = W - XH
card(XH, Y_IN, WH, H_IN)
t8(a, XH + PADX, y1_, "Hearts", fontweight="bold", zorder=5)
c_lab = XH + PADX + tw("Hearts", FS_S, "bold") + 3.0
c1 = c_lab + max(tw("HFpEF"), tw("DCM")) + 1.6
c2 = c1 + 3.0 + max(tw("SCP3342"), tw("GSE290367")) + 2.6
c3 = c2 + 3.0 + max(tw("Hahn biopsy"), tw("MAGNet")) + 4.0
for yy_, grp, (k1, n1), (k2, n2) in [(y1_, "HFpEF", ("SCP3342", "SCP3342"), ("Hahn", "Hahn biopsy")),
                                     (y2_, "DCM", ("GSE290367", "GSE290367"), ("MAGNet", "MAGNet"))]:
    t8(a, c_lab, yy_, grp, color="#555555", zorder=5)
    mk(a, c1 + 1.1, yy_, k1, z=5); t8(a, c1 + 3.0, yy_, n1, zorder=5)
    mk(a, c2 + 1.1, yy_, k2, z=5); t8(a, c2 + 3.0, yy_, n2, zorder=5)
a.plot([c3 - 2.0, c3 - 2.0], [Y_IN + 1.6, Y_IN + H_IN - 1.6], color="#D0D0D0", lw=0.6, zorder=4)
mk(a, c3 + 1.1, y1_, "GTEx", z=5); t8(a, c3 + 3.0, y1_, "GTEx v11", zorder=5)
t8(a, c3 + 3.0, y2_, "(ventilator proxy)", zorder=5)
fits(c3 + 3.0, "(ventilator proxy)", W - PADX)
# inputs -> steps 1 and 2
arrowA(SX[0] + CW / 2, Y_IN + H_IN, SX[0] + CW / 2, Y_ST)
assert XH + 2.0 < SX[1] + CW / 2
arrowA(SX[1] + CW / 2, Y_IN + H_IN, SX[1] + CW / 2, Y_ST)

# --- 3. analysis row: 4 equal step cards
x = step_card(0, "GC reference", accent=GC_CLASS["induced"])
k = body_lines(x, ["Random-effects", "meta-analysis"])
t8(a, x + PADX, body_y(k), "13,368 genes classified", fontweight="bold", zorder=5)
fits(x + PADX, "13,368 genes classified", x + CW - PADX, weight="bold")
SW = 2.0
cls_items = [("induced", "Induced"), ("repressed", "Repressed"), ("insensitive", "Insensitive"),
             ("ambiguous", "Ambiguous")]
colx = [x + 1.0, x + 1.0 + SW + 0.5 + max(tw("Induced"), tw("Insensitive")) + 1.8]
for j, (kk, lab) in enumerate(cls_items):
    xx_ = colx[j % 2]; yy_ = body_y(k + 1 + j // 2)
    a.add_patch(Rectangle((xx_, yy_ - SW / 2), SW, SW, fc=GC_CLASS[kk], lw=0.6, zorder=5,
                          ec=AMBIGUOUS_EDGE if kk == "ambiguous" else "none"))
    t8(a, xx_ + SW + 0.5, yy_, lab, zorder=5)
    assert xx_ + SW + 0.5 + tw(lab) <= x + CW - 0.8, lab
chips(x + CW, CHIP_Y, ["Fig 2", "Fig 6"])
assert body_y(k + 2) + PIT <= CHIP_Y + 1e-6
reg("A", "4 in vitro datasets", 4, "S2_datasets", "rows with role == 'in vitro GC reference'", "yes")
reg("A", "13,368 genes classified", 13368, "S5_GC_classification", "row count", "yes")

x = step_card(1, "Scores")
body_lines(x, ["GC score", "8-gene GC response", "Circadian, cAMP", "Macrophage modules"])
chips(x + CW, CHIP_Y, ["Fig 1E", "Fig 7"])

x = step_card(2, "Contrast")
k = body_lines(x, wrap("Control − case, per cohort and cell type", CW - 2 * PADX))
k = body_lines(x, wrap("± circadian or cAMP score", CW - 2 * PADX), k)
assert body_y(k - 1) + PIT <= CHIP_Y + 1e-6
chips(x + CW, CHIP_Y, ["Fig 3", "Fig 4"])

x = step_card(3, "Identifiable?", ls=(0, (2.4, 1.3)))
k = body_lines(x, wrap("AUC and variance inflation of the disease term", CW - 2 * PADX))
t8(a, x + PADX, body_y(k), "Added post hoc", fontstyle="italic", color="#555555", zorder=5)
if A_VARIANT == 2:                                    # decision diamond in the chip slot
    dx_, dy_ = x + CW - 4.6, CHIP_Y - 0.6
    a.add_patch(Polygon([[dx_ - 3.4, dy_], [dx_, dy_ - 2.4], [dx_ + 3.4, dy_], [dx_, dy_ + 2.4]], closed=True,
                        fc="white", ec=TEXT, lw=0.7, ls=(0, (1.6, 0.9)), zorder=5))
for i in range(3):                                    # step arrows (mid-height)
    arrowA(SX[i] + CW, Y_ST + H_ST / 2, SX[i + 1], Y_ST + H_ST / 2)

# --- 4. outcome row
yes_t1, yes_t2 = "adjust for GC score;", "audit genes"
no_t1, no_t2 = "report as non-identifiable; illustrative, post hoc", "AUC ≥0.9 or ≤0.1, or VIF ≥5"
if A_VARIANT == 1:        # Yes under step 4 (straight down); No under steps 2-3 (short elbow); key under step 1
    YX, YW = SX[3], CW
    NX, NW_ = SX[1], 2 * CW + GAP
else:                     # grid halves: Yes under steps 3-4, No under steps 1-2 (bus to the left)
    YX, YW = SX[2], 2 * CW + GAP
    NX, NW_ = SX[0], 2 * CW + GAP
oy1, oy2 = Y_OUT + H_OUT / 2 - PIT / 2, Y_OUT + H_OUT / 2 + PIT / 2
OCHIP_Y = oy2 + 0.1
if YW > 60:                                           # wide Yes card: one text line
    yes_t1, yes_t2 = "adjust for GC score; audit genes", ""
for (xo, wo, lab, t1, t2, ch) in [(YX, YW, "Yes", yes_t1, yes_t2, "Fig 5E–I"), (NX, NW_, "No", no_t1, no_t2, "Fig 5A–D")]:
    card(xo, Y_OUT, wo, H_OUT)
    t8(a, xo + PADX, oy1, lab, fontweight="bold", zorder=5)
    xl_ = xo + PADX + tw(lab, FS_S, "bold") + 1.4
    t8(a, xl_, oy1, t1, zorder=5); fits(xl_, t1, xo + wo - PADX)
    t8(a, xo + PADX, oy2, t2, zorder=5)
    c0_ = chips(xo + wo, OCHIP_Y, [ch])
    assert xo + PADX + (tw(t2) if t2 else 0) < c0_ - 0.8, (lab, t2)
ybus = Y_ST + H_ST + 2.2
arrowA(SX[3] + CW / 2, Y_ST + H_ST, SX[3] + CW / 2, Y_OUT)                     # Yes: straight down
if A_VARIANT == 1:
    xs_no = SX[3] + 4.0; xe_no = SX[2] + CW - 12.0
else:
    xs_no = SX[3] + 4.0; xe_no = NX + NW_ - 12.0
polyarrow([(xs_no, Y_ST + H_ST), (xs_no, ybus), (xe_no, ybus), (xe_no, Y_OUT)])
assert xe_no < NX + NW_ - 1.0 and xe_no > NX + 1.0
if A_VARIANT == 1:      # key for the dashed outline, in the free slot under step 1
    kx, ky_ = SX[0] + 1.0, Y_OUT + H_OUT / 2
    a.add_patch(FancyBboxPatch((kx, ky_ - 1.6), 5.0, 3.2, boxstyle="round,pad=0,rounding_size=0.6", fc="white",
                               ec=BORDER, lw=0.6, ls=(0, (2.4, 1.3)), zorder=4))
    t8(a, kx + 6.2, ky_, "Dashed: post hoc", color="#555555", zorder=5)
    fits(kx + 6.2, "Dashed: post hoc", SX[0] + CW - 1.0)
for s_, v_ in [("AUC ≥0.9 or ≤0.1", "0.9; 0.1"), ("VIF ≥5", 5)]:
    reg("A", s_, v_, "Manuscript v0.4 Methods (Identifiability)", "threshold as stated", "yes")
reg("A", "refs 10–12", "10-12", "Manuscript v0.4 Introduction", "reference numbers as cited", "yes")
SRC["1A_schematic"] = pd.DataFrame({
    "element": ["context: donor controls", "input: in vitro datasets", "input: hearts", "step 1 GC reference",
                "step 2 scores", "step 3 contrast", "step 4 identifiability (added post hoc)", "outcome: yes",
                "outcome: no", "identifiability AUC threshold", "identifiability VIF threshold"],
    "value": ["ICU management often includes corticosteroids (refs 10-12); exposure not recorded in any cohort",
              "4 datasets; human monocytes/macrophages; GC vs vehicle",
              "HFpEF: SCP3342, Hahn biopsy; DCM: GSE290367, MAGNet; GTEx v11 (ventilator proxy)",
              "random-effects meta-analysis; 13,368 genes classified (induced, repressed, insensitive, ambiguous); "
              "Fig 2, Fig 6",
              "GC score; 8-gene GC response; circadian; cAMP; macrophage modules; Fig 1E, Fig 7",
              "control - case per cohort and cell type; +/- circadian or cAMP score; Fig 3, Fig 4",
              "AUC and variance inflation of the disease term",
              "adjust for GC score; audit genes; Fig 5E-I",
              "report as non-identifiable; illustrative, post hoc; Fig 5A-D",
              "AUC >=0.9 or <=0.1", "VIF >=5"],
    "source": ["Manuscript v0.4 Introduction", "S2_datasets", "S2_datasets", "S5_GC_classification (n rows)",
               "S3_gene_sets", "Methods", "Methods (Identifiability)", "Methods", "Methods",
               "Manuscript v0.4 Methods", "Manuscript v0.4 Methods"]})

# ======================================================================= B  cohort table
BW_ = 125.0
cB = box(0, Y_BC, BW_, BC_H)
header(cB, "B", "Cohorts and their design constraints")
b = canvas(cB)
# GTEx n from the per-donor files (analysis set; vent = Hardy 0); cross-checked against the stored decision output
GT_LV = pd.read_csv(PD + "v12_rerun/perdonor_scores_GTEx_LV.csv"); GT_LV = GT_LV[GT_LV.analysis_ok]
GT_AA = pd.read_csv(PD + "v12_rerun/perdonor_scores_GTEx_AA.csv"); GT_AA = GT_AA[GT_AA.analysis_ok]
gtex_n = {"GTEx_LV": int((GT_LV.vent == 0).sum()), "GTEx_AA": int((GT_AA.vent == 0).sum())}   # other deaths
_stored = {d["dataset"]: d["n"] for d in dec_g25["G5"]["rows"]}
assert all(_stored[k] == v for k, v in gtex_n.items()), (_stored, gtex_n)
assert len(GT_LV) == 450 and len(GT_AA) == 461
assert 450 - gtex_n["GTEx_LV"] == 253 and 461 - gtex_n["GTEx_AA"] == 215
trows = [
    ("SCP3342", "SCP3342", "Septum\nsnRNA-seq", "24", "19 HFpEF", "Pools of 5–6 donors; both groups per pool"),
    ("Hahn", "Hahn biopsy", "RV septum\nbulk RNA-seq", "24", "41 HFpEF", "Same program; overlaps SCP3342"),
    ("GSE290367", "GSE290367", "LV\nsnRNA-seq", "139", "67 DCM", "10 studies; DCM from 3 (panel C)"),
    ("MAGNet", "MAGNet", "LV\nbulk RNA-seq", "102", "128 DCM", "GSE290367 overlap excluded"),
    ("GTEx", "GTEx v11", "LV, AA\nbulk RNA-seq", f"LV {450 - gtex_n['GTEx_LV']}\nAA {461 - gtex_n['GTEx_AA']}",
     f"LV {gtex_n['GTEx_LV']}\nAA {gtex_n['GTEx_AA']}", "Ventilator deaths as controls; not blinded"),
    (None, "In vitro ×4", "Monocyte, MDM;\narray, RNA-seq", "vehicle", "GC", "6, 4, 4 and 3 donors (Fig 2A)"),
]
heads = ["Cohort", "Tissue / assay", "Controls", "Cases", "Design constraint"]
def mw(vals, bold=False): return max(tw(z, FS_S, "bold" if bold else "normal") for v in vals for z in v.split("\n"))
cx_name = 4.8
cx_tis = cx_name + max(mw([r_[1] for r_ in trows], True), tw(heads[0], FS_S, "bold")) + 2.2
cx_ctl = cx_tis + max(mw([r_[2] for r_ in trows]), tw(heads[1], FS_S, "bold")) + 2.2 + \
    max(mw([r_[3] for r_ in trows]), tw(heads[2], FS_S, "bold"))
cx_cas = cx_ctl + 2.2 + max(mw([r_[4] for r_ in trows]), tw(heads[3], FS_S, "bold"))
cx_con = cx_cas + 2.2
con_w = BW_ - cx_con - 0.3
hy0 = 9.6
for x_, lab, ha in [(cx_name, heads[0], "left"), (cx_tis, heads[1], "left"), (cx_ctl, heads[2], "right"),
                    (cx_cas, heads[3], "right"), (cx_con, heads[4], "left")]:
    t8(b, x_, hy0, lab, fontweight="bold", ha=ha, color="#444444")
fits(cx_con, heads[4], BW_ - 0.3, weight="bold")
b.plot([0.4, BW_], [hy0 + 2.4, hy0 + 2.4], color=TEXT, lw=0.6)
row0 = hy0 + 3.0
row_h = (BC_H - row0 - 0.2) / len(trows)
assert row_h >= 2 * LN + 0.2, row_h
rec = []
for i, (coh, name, tissue, ctrl, case, constr) in enumerate(trows):
    yy = row0 + (i + 0.5) * row_h
    if i % 2 == 1:
        b.add_patch(Rectangle((0.4, yy - row_h / 2), BW_ - 0.4, row_h, fc="#F4F4F4", ec="none", zorder=0))
    if coh:
        mk(b, cx_name - 1.9, yy, coh, s=4.4)
    else:
        b.add_patch(Circle((cx_name - 1.9, yy), 1.3, fc="white", ec=TEXT, lw=0.6, zorder=3))
        b.add_patch(Circle((cx_name - 1.9, yy), 0.5, fc=GC_CLASS["induced"], ec="none", zorder=4))
    t8(b, cx_name, yy, name, fontweight="bold")
    t8(b, cx_tis, yy, tissue, linespacing=1.0)
    t8(b, cx_ctl, yy, ctrl, ha="right", linespacing=1.0)
    t8(b, cx_cas, yy, case, ha="right", linespacing=1.0)
    cl_ = wrap(constr, con_w)
    assert len(cl_) <= 2, cl_
    t8(b, cx_con, yy, "\n".join(cl_), linespacing=1.0)
    if coh == "GTEx":
        ctrl, case = f"ventilator-managed: {ctrl}", f"other deaths: {case}"
    rec.append(dict(cohort=name, tissue_assay=tissue.replace("\n", ", "), controls=ctrl.replace("\n", "; "),
                    cases=case.replace("\n", "; "), design_constraint=constr))
y_a = row0 + 0.1 * row_h; y_b = row0 + 1.9 * row_h
b.plot([0.9, 0.45, 0.45, 0.9], [y_a, y_a, y_b, y_b], color=COHORT["SCP3342"], lw=0.9, solid_capstyle="butt")
for txt, val, src, comp, m_ in [
        ("SCP3342 24 / 19 HFpEF", "24; 19", "S2_datasets", "as tabulated", "yes"),
        ("SCP3342 pools of 5–6 donors", "5-6", "S2_datasets", "as tabulated", "yes"),
        ("Hahn biopsy 24 / 41 HFpEF", "24; 41", "S2_datasets", "as tabulated", "yes"),
        ("GSE290367 139 / 67 DCM", "139; 67", "pseudobulk_L1_meta.csv + duplicate_individuals.csv",
         "LV, ND/DCM, age>=20, postnatal, drop duplicate records, unique donors", "yes"),
        ("GSE290367 10 studies; DCM from 3", "10; 3", "pseudobulk_L1_meta.csv", "unique studies; studies with DCM>0",
         "yes"),
        ("MAGNet 102 / 128 DCM", "102; 128", "S2_datasets", "as tabulated", "yes"),
        ("GTEx controls (ventilator-managed) LV 253, AA 215; cases (other deaths) LV 197, AA 246", "253; 215; 197; 246",
         "S2_datasets (450, 253); perdonor_scores_GTEx_LV.csv / _AA.csv (analysis set, vent = Hardy 0)",
         "other = total - ventilator", "new (S: 450, 253, 461, 215 and 197 are in the text; AA other 246 = AA analysis set "
         "minus ventilator deaths)"),
        ("In vitro 6, 4, 4 and 3 donors", "6; 4; 4; 3", "S2_datasets", "as tabulated", "new (S)")]:
    reg("B", txt, val, src, comp, m_)
assert gse_ct.index.nunique() == 10 and (gse_ct.DCM > 0).sum() == 3
SRC["1B_cohorts"] = pd.DataFrame(rec)

# ======================================================================= C  GSE290367 study composition
cC = box(BW_ + GUT, Y_BC, W - BW_ - GUT, BC_H)
th = header(cC, "C", "GSE290367: DCM from 3 of 10 studies")
both = gse_ct[(gse_ct.ND > 0) & (gse_ct.DCM > 0)].sort_values("DCM", ascending=False)
only = gse_ct[gse_ct.DCM == 0].sort_values("ND", ascending=False)
order = list(both.index) + list(only.index)
short = {s: ("ENCODE v4" if s.startswith("ENCODE") else s.split(" ")[0]) for s in order}
kc = canvas(cC)
ky = th + 2.6
xx = 0.6
for col, lab in [(CONTROL, "Non-failing (139)"), (COHORT["GSE290367"], "DCM (67)")]:
    kc.add_patch(Rectangle((xx, ky - 1.3), 2.6, 2.6, fc=col, ec="none"))
    t8(kc, xx + 3.3, ky, lab)
    xx += 3.3 + tw(lab) + 2.6
assert xx - 2.6 < cC["w"]
reg("C", "Non-failing (139); DCM (67)", "139;67", "pseudobulk_L1_meta.csv", "column totals", "yes")
labw = max(tw(short[s]) for s in order) + 1.8
ay0 = ky + 3.4
axC = axmm(cC["x"] + labw, cC["y"] + ay0, cC["w"] - labw - 0.8, cC["h"] - ay0 - 9.4)
yy = np.arange(len(order))
nd = gse_ct.loc[order, "ND"].values; dc = gse_ct.loc[order, "DCM"].values
axC.axhspan(-0.5, len(both) - 0.5, color="#E3F1EF", lw=0, zorder=0)
axC.barh(yy, nd, color=CONTROL, height=0.7, lw=0, zorder=2)
axC.barh(yy, dc, left=nd, color=COHORT["GSE290367"], height=0.7, lw=0, zorder=2)
for i in range(len(order)):
    lab = f"{nd[i]} + {dc[i]}" if dc[i] else f"{nd[i]}"
    axC.text(nd[i] + dc[i] + 1.2, i, lab, fontsize=FS_S, va="center", ha="left", zorder=3)
    reg("C", f"{short[order[i]]} {lab}", f"{nd[i]};{dc[i]}", "pseudobulk_L1_meta.csv + duplicate_individuals.csv",
        "unique adult LV donors per study after de-duplication (ND;DCM)", "new (N)")
axC.set_yticks(yy); axC.set_yticklabels([short[s] for s in order], fontsize=FS_S)
axC.set_ylim(len(order) - 0.5, -0.5)
axC.set_xlim(0, 78); axC.set_xticks([0, 25, 50, 75])
axC.set_xlabel("Donors (n)", labelpad=1.0)
axC.tick_params(axis="y", length=0, pad=1.5)
assert both.ND.sum() == 45 and both.DCM.sum() == 67
reg("C (legend)", "shaded studies used for DCM contrasts: 45 non-failing vs 67 DCM", "45;67", "pseudobulk_L1_meta.csv",
    "sum over studies with both groups", "new (S: S8a n/n_case)")
SRC["1C_GSE290367_studies"] = pd.DataFrame({"study": order, "nonfailing_n": nd, "DCM_n": dc,
                                            "used_for_DCM_contrasts": [s in both.index for s in order]})

# ======================================================================= D  SCP3342 pools
DW = 46.0
EW = W - DW - GUT
cD = box(0, Y_R3, DW, R3_H)
th = header(cD, "D", "SCP3342: both groups in every pool")
ylab_w = max(tw("Control"), tw("(n=24)")) + 1.6
axD = axmm(cD["x"] + ylab_w, cD["y"] + th + 1.6, DW - ylab_w - 0.8, R3_H - th - 1.6 - 10.6)
pools = sorted(scp.pool_n.unique())
recD = []
for j, p in enumerate(pools):
    for r_, (cv, col) in enumerate([(0, CONTROL), (1, COHORT["SCP3342"])]):
        sub = scp[(scp.pool_n == p) & (scp.case == cv)].sort_values("unit_id")
        k = len(sub)
        for q in range(k):
            axD.plot(j, r_ + (q - (k - 1) / 2) * 0.27, "o", ms=4.0, mfc=col, mec=TEXT, mew=0.35)
        recD += [dict(pool=f"pool{p}", group="HFpEF" if cv else "control", unit_id=u) for u in sub.unit_id]
axD.set_xlim(-0.6, len(pools) - 0.4); axD.set_ylim(1.6, -0.6)
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
cE = box(DW + GUT, Y_R3, EW, R3_H)
th = header(cE, "E", "Score genes and their in vitro GC class")
e = canvas(cE)
cls = S5.set_index("gene")["class"]
cmap_cls = {"GC_induced": "induced", "GC_repressed": "repressed", "GC_insensitive": "insensitive",
            "ambiguous": "ambiguous"}
gs3 = S3.groupby("gene_set").gene.apply(list)
gc3 = gs3["GC_specific3"]; gcr = [x for x in gs3["GC_response"] if x not in gc3]
xr = EW - 0.3
TH_, TP = 4.4, 5.0                                  # tile height, row pitch (mm)
recE = []


def tile(x, y, gname, k):
    wbox = tw(gname, FS_S, style="italic") + 1.2
    e.add_patch(Rectangle((x, y - TH_ / 2), wbox, TH_, fc=GC_CLASS[k],
                          ec=(AMBIGUOUS_EDGE if k == "ambiguous" else ("#9A9A9A" if k == "not_measured" else "none")),
                          lw=0.6, hatch="////" if k == "not_measured" else None, zorder=2))
    e.text(x + wbox / 2, y + 0.05, gname, fontsize=FS_S, ha="center", va="center", fontstyle="italic", zorder=3)
    return wbox


def tile_row(y, x, genes, score):
    xs_ = {}
    for gname in genes:
        c_ = cls.get(gname, None); k = cmap_cls.get(c_, "not_measured")
        wb_ = tile(x, y, gname, k); xs_[gname] = (x, x + wb_)
        recE.append(dict(score=score, gene=gname,
                         in_vitro_class=c_ if c_ else "not measured in >=2 in vitro datasets (not in S5)"))
        x += wb_ + 0.45
    assert x - 0.45 <= xr, (score, x)
    return x, xs_


labw = max(tw("Morning +"), tw("Evening −"), tw("cAMP", FS_S, "bold"), tw("GC score", FS_S, "bold")) + 1.6
yy = th + 3.4
# GC score (3 genes)
t8(e, 0.4, yy, "GC score", fontweight="bold")
tile_row(yy, 0.4 + labw, gc3, "GC score")
# 8-gene GC response = GC score + 5 genes, one line + bracket
yy += TP + 1.0
t8(e, 0.4, yy, "8-gene GC response", fontweight="bold")
x0_ = 0.4 + tw("8-gene GC response", FS_S, "bold") + 1.2
t8(e, x0_, yy, "= GC score +", color="#555555")
x_, xs_ = tile_row(yy, x0_ + tw("= GC score +") + 1.2, gcr, "8-gene GC response")
x_beg_ex, x_end_ex = xs_["PER1"][0], xs_["DDIT4"][1]
yb = yy + TH_ / 2 + 0.5
e.plot([x_beg_ex, x_beg_ex, x_end_ex, x_end_ex], [yb - 0.4, yb, yb, yb - 0.4], color=TEXT, lw=0.6)
t8(e, x_beg_ex, yb + 2.3, "not in GC score", color="#444444")
# circadian arms (score = morning minus evening)
yy = yb + 2.3 + LN + 0.6
for lab, genes in [("Morning +", gs3["CLOCK_morning"]), ("Evening −", gs3["CLOCK_evening"])]:
    t8(e, 0.4, yy, lab)
    tile_row(yy, 0.4 + labw, genes, "Circadian " + lab.split()[0].lower() + " arm")
    yy += TP
# cAMP (one line)
yy += 0.2
t8(e, 0.4, yy, "cAMP", fontweight="bold")
tile_row(yy, 0.4 + labw, gs3["CAMP_IEG"], "cAMP")
# class key (flowed)
ky = yy + TP + 0.4; xx = 0.4
for k, lab in [("induced", "Induced"), ("insensitive", "Insensitive"), ("ambiguous", "Ambiguous"),
               ("not_measured", "Not measured in ≥2 datasets")]:
    wd = 3.3 + tw(lab) + 2.4
    if xx + wd - 2.4 > xr:
        xx = 0.4; ky += LN
    e.add_patch(Rectangle((xx, ky - 1.3), 2.6, 2.6, fc=GC_CLASS[k], lw=0.6,
                          ec=AMBIGUOUS_EDGE if k == "ambiguous" else ("#9A9A9A" if k == "not_measured" else "none"),
                          hatch="////" if k == "not_measured" else None))
    t8(e, xx + 3.3, ky, lab)
    xx += wd
assert ky + 2.2 <= R3_H, ky
assert not any(r_["in_vitro_class"] == "GC_repressed" for r_ in recE)
SRC["1E_score_genes"] = pd.DataFrame(recE)

# ======================================================================= F / G  donor demographics (per-donor bundle)
# One row per donor: the score files hold one row per donor x cell type (SCP3342, GSE290367) -> de-duplicate.
B_N = {"SCP3342": (24, 19), "Hahn": (24, 41), "Chaffin": (16, 11), "Koenig": (23, 10), "Reichart": (6, 46),
       "GSE290367": (45, 67), "MAGNet": (102, 128), "GTEx LV": (253, 197)}     # panel B / C counts (control, case)
d_scp = pd.read_csv(PD + "v13_rerun/perdonor_scores_SCP3342_HFpEF.csv")
d_hahn = pd.read_csv(PD + "v13_rerun/perdonor_scores_Hahn2021_RVS_HFpEF.csv")
d_mag = pd.read_csv(PD + "v13_rerun/perdonor_scores_MAGNet_exclOverlap_DCM.csv")
d_gse = pd.read_csv(PD + "v13_rerun/perdonor_scores_GSE290367_DCM.csv")
for d_, k_ in [(d_scp, ["sex", "case"]), (d_hahn, ["age", "sex", "case"]), (d_mag, ["age", "gender", "case"]),
               (d_gse, ["age", "sex", "case", "study"])]:
    assert (d_.groupby("unit_id")[k_].nunique(dropna=False) <= 1).all().all()   # demographics constant within donor
d_scp = d_scp.drop_duplicates("unit_id"); d_hahn = d_hahn.drop_duplicates("unit_id")
d_mag = d_mag.drop_duplicates("unit_id"); d_gse = d_gse.drop_duplicates("unit_id")
assert GT_LV.SUBJID.is_unique
don = []
def _add(df, cohort, study, group, fem, age, note):
    for u, f_, a_ in zip(df.unit_id if "unit_id" in df else df.SUBJID, fem, age):
        don.append(dict(cohort=cohort, study=study, group=group, donor_id=str(u), female=int(f_),
                        age_years=a_, age_note=note))
for cv, gname in [(0, "control"), (1, "HFpEF")]:
    s_ = d_scp[d_scp.case == cv]; _add(s_, "SCP3342", "", gname, s_.female, [np.nan] * len(s_), "per-donor age not available")
    s_ = d_hahn[d_hahn.case == cv]; _add(s_, "Hahn", "", gname, 1 - s_.male, s_.age, "years")
for cv, gname in [(0, "control"), (1, "DCM")]:
    for st in ["Chaffin", "Koenig", "Reichart"]:
        s_ = d_gse[(d_gse.case == cv) & d_gse.study.str.startswith(st)]
        _add(s_, "GSE290367", st, gname, s_.female, s_.age, "years as deposited")
    s_ = d_mag[d_mag.case == cv]; _add(s_, "MAGNet", "", gname, 1 - s_.male, s_.age, "years")
for v_, gname in [(1, "ventilator death"), (0, "other death")]:
    s_ = GT_LV[GT_LV.vent == v_]; _add(s_, "GTEx LV", "", gname, 1 - s_.male, s_.age, "decade-bin midpoint (" + s_.AGE + ")")
DON = pd.DataFrame(don)
DON["role"] = np.where(DON.group.isin(["control", "ventilator death"]), "control", "case")
# n check against panel B / C (report, do not force)
NCHK = []
for key, (nc, nk) in B_N.items():
    if key in ("Chaffin", "Koenig", "Reichart"):
        sub = DON[DON.study == key]
    else:
        sub = DON[DON.cohort == key]
    oc, ok_ = int((sub.role == "control").sum()), int((sub.role == "case").sum())
    NCHK.append(dict(cohort=key, panel_B_controls=nc, panel_B_cases=nk, perdonor_controls=oc, perdonor_cases=ok_,
                     match=(oc == nc and ok_ == nk)))
NCHK = pd.DataFrame(NCHK)
print(NCHK.to_string(index=False))
assert NCHK.match.all(), "n mismatch vs panel B - report"
SRC["1F_1G_n_check"] = NCHK

R4_TOP = Y_R4
# ---- F layout (mm)
F_AXL = 10.4; SCP_W = 16.6; PP = 15.2; COLD = 3.7
pairs = [("Hahn", "Hahn", ""), ("Chaffin", "GSE290367", "Chaffin"), ("Koenig", "GSE290367", "Koenig"),
         ("Reichart", "GSE290367", "Reichart"), ("MAGNet", "MAGNet", ""), ("GTEx LV", "GTEx LV", "")]
AXW = SCP_W + len(pairs) * PP
FW = F_AXL + AXW + 0.4
GW = W - FW - GUT
cF = box(0, Y_R4, FW, R4_H)
th = header(cF, "F", "Donor age, controls vs cases")
fc_ = canvas(cF)
AX_T, AX_B = th + 1.4, 13.6                       # top offset, bottom margin (n / name / study tiers)
axF = axmm(F_AXL, Y_R4 + AX_T, AXW, R4_H - AX_T - AX_B)
axF.set_xlim(0, AXW); axF.set_ylim(10, 92)
axF.set_yticks([20, 40, 60, 80]); axF.set_ylabel("Age (years)", labelpad=1.0)
axF.set_xticks([]); axF.spines["bottom"].set_visible(True)
rng = np.random.default_rng(20261020)
ph_mm = R4_H - AX_T - AX_B
recFs = []
def ymm(v):                                        # data y -> canvas mm (panel F)
    return AX_T + (92 - v) / 82 * ph_mm
# SCP3342 block
axF.text(SCP_W / 2, 51, "per-donor\nage not\navailable", ha="center", va="center", fontsize=FS_S, fontstyle="italic",
         color="#555555", linespacing=1.05)
axF.axvline(SCP_W, color=GRID, lw=0.6, zorder=0)
t8(fc_, F_AXL + SCP_W / 2, R4_H - AX_B + 7.9, "SCP3342", ha="center", va="baseline")
t8(fc_, F_AXL - 1.2, R4_H - AX_B + 3.7, "n", ha="right", va="baseline", fontstyle="italic")
reg("F", "SCP3342: per-donor age not available", "n/a", "perdonor_scores_SCP3342_HFpEF.csv",
    "no age column in the per-donor file (ages not public)", "new")
for j, (lab, coh, st) in enumerate(pairs):
    xc = SCP_W + (j + 0.5) * PP
    if j:
        axF.axvline(SCP_W + j * PP, color=GRID, lw=0.6, zorder=0)
    sub = DON[(DON.study == st)] if st else DON[DON.cohort == coh]
    for side, role in [(-1, "control"), (1, "case")]:
        ss = sub[sub.role == role]
        x0 = xc + side * COLD
        big = len(ss) > 80
        if coh == "GTEx LV":
            col = COHORT["GTEx"] if role == "control" else CONTROL
            mkr = COHORT_MARKER["GTEx"]
        else:
            col = CONTROL if role == "control" else COHORT[coh]
            mkr = COHORT_MARKER[coh]
        jit = rng.uniform(-1, 1, len(ss)) * (2.9 if big else 2.3)
        axF.plot(x0 + jit, ss.age_years, ls="none", marker=mkr, ms=1.9 if big else 2.8, mfc=col,
                 mec="none" if big else TEXT, mew=0.3, alpha=0.75 if big else 1.0, zorder=2)
        med = float(np.median(ss.age_years))
        axF.plot([x0 - 2.6, x0 + 2.6], [med, med], color=TEXT, lw=1.4, solid_capstyle="butt", zorder=4)
        q1, q3 = np.percentile(ss.age_years, [25, 75])
        recFs.append(dict(cohort=coh, study=st, group=ss.group.iloc[0], role=role, n=len(ss), median_age=med,
                          q1=q1, q3=q3, min=ss.age_years.min(), max=ss.age_years.max(),
                          note=("GTEx ages are decade bins; points at bin midpoints" if coh == "GTEx LV" else "")))
        t8(fc_, F_AXL + x0, R4_H - AX_B + 3.7, str(len(ss)), ha="center", va="baseline")
        reg("F", f"{lab} {ss.group.iloc[0]} n {len(ss)}", len(ss), "per-donor bundle (" +
            {"Hahn": "perdonor_scores_Hahn2021_RVS_HFpEF.csv", "GSE290367": "perdonor_scores_GSE290367_DCM.csv",
             "MAGNet": "perdonor_scores_MAGNet_exclOverlap_DCM.csv", "GTEx LV": "perdonor_scores_GTEx_LV.csv"}[coh] + ")",
            "unique donors" + (" (analysis set; vent = Hardy 0)" if coh == "GTEx LV" else ""), "yes (= panel B/C)")
    t8(fc_, F_AXL + xc, R4_H - AX_B + 7.9, lab, ha="center", va="baseline")
    assert tw(lab) < PP - 0.6, lab
# GSE290367 study tier
x_a = F_AXL + SCP_W + 1 * PP + 0.8; x_b = F_AXL + SCP_W + 4 * PP - 0.8
yt = R4_H - AX_B + 9.2
fc_.plot([x_a, x_a, x_b, x_b], [yt - 0.5, yt, yt, yt - 0.5], color=TEXT, lw=0.6)
t8(fc_, (x_a + x_b) / 2, yt + 3.6, "GSE290367", ha="center", va="baseline")
assert yt + 2.3 + 2.0 <= R4_H + 0.05
SRC["1F_donor_age"] = DON[DON.cohort != "SCP3342"].drop(columns="female")
SRC["1F_age_summary"] = pd.DataFrame(recFs)

# ---- G  sex composition
cG = box(FW + GUT, Y_R4, GW, R4_H)
th = header(cG, "G", "Sex composition, controls vs cases")
gx_ = canvas(cG)
rowsG = [("SCP3342", "SCP3342", None), ("Hahn", "Hahn", None), ("GSE290367", "GSE290367", None),
         ("MAGNet", "MAGNet", None), ("GTEx LV", "GTEx LV", None)]
LABW = max(tw(r_[0]) for r_ in rowsG) + 1.4
NW = tw("253") + 0.8
BW = (GW - 0.3 - LABW - 2 * (NW + 1.6)) / 2
FEM, MAL = "#9DB2C6", "#E6E6E6"
xcol = [LABW + NW, LABW + NW + BW + 1.6 + NW]
yh = th + 2.6
for xb_, lab in zip(xcol, ["Control", "Case"]):
    t8(gx_, xb_ + (BW - NW) / 2, yh, lab, ha="center", fontweight="bold", color="#444444")   # over n + bar
recG = []
RP = 4.8; BH = 3.4
for i, (lab, coh, _) in enumerate(rowsG):
    yy = yh + 1.4 + (i + 0.5) * RP
    t8(gx_, 0.2, yy, lab)
    sub = DON[DON.cohort == coh]
    for xb_, role in zip(xcol, ["control", "case"]):
        ss = sub[sub.role == role]
        n_, nf_ = len(ss), int(ss.female.sum()); p_ = 100 * nf_ / n_
        fw_ = BW * nf_ / n_
        gx_.add_patch(Rectangle((xb_, yy - BH / 2), fw_, BH, fc=FEM, ec="none", zorder=2))
        gx_.add_patch(Rectangle((xb_ + fw_, yy - BH / 2), BW - fw_, BH, fc=MAL, ec="none", zorder=2))
        s_ = f"{p_:.0f}%"
        assert tw(s_) + 1.0 <= BW, (lab, role)       # label centred on the bar; both fills light enough for dark text
        t8(gx_, xb_ + BW / 2, yy + 0.05, s_, ha="center", color=TEXT, zorder=3)
        t8(gx_, xb_ - 0.6, yy, str(n_), ha="right")
        grp = ss.group.iloc[0]
        recG.append(dict(cohort=coh, group=grp, role=role, n=n_, n_female=nf_, n_male=n_ - nf_, pct_female=p_,
                         shown_in_figure=True))
        reg("G", f"{lab} {grp}: {s_} female (n={n_})", f"{p_:.4f}; {n_}", "per-donor bundle",
            "female / unique donors (GSE290367: the 3 contrast studies pooled)", "new")
    assert tw(lab) <= LABW - 1.0
# per-study GSE290367 rows (source data only)
for st in ["Chaffin", "Koenig", "Reichart"]:
    for role in ["control", "case"]:
        ss = DON[(DON.study == st) & (DON.role == role)]
        recG.append(dict(cohort="GSE290367 " + st, group=ss.group.iloc[0], role=role, n=len(ss),
                         n_female=int(ss.female.sum()), n_male=int(len(ss) - ss.female.sum()),
                         pct_female=100 * ss.female.mean(), shown_in_figure=False))
ky = yh + 1.4 + len(rowsG) * RP + 2.6; xx = LABW + NW
for col, lab in [(FEM, "Female"), (MAL, "Male")]:
    gx_.add_patch(Rectangle((xx, ky - 1.3), 2.6, 2.6, fc=col, ec="none"))
    t8(gx_, xx + 3.3, ky, lab); xx += 3.3 + tw(lab) + 2.6
assert xx - 2.6 <= GW and ky + 2.2 <= R4_H, (xx, ky)
SRC["1G_sex"] = pd.DataFrame(recG)

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
print("font sizes used (pt):", _sizes)

_bad = min_font_audit(fig, 10); print("min_font_audit:", _bad); assert _bad == []

# ------------------------------------------------------------------ save
save_figure(fig, "Figure_1", tight=False)   # canvas drawn at exactly 183 mm
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
