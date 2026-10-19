"""Central Illustration / Visual Abstract (storyboard section 6). Stylised: no figure panel is duplicated.
Every printed number is a locked v0.4 manuscript value; the script verifies each phrase against the text.
Run: python ci.py  -> Central_Illustration.{pdf,tif,png}, Central_Illustration_number_registry.csv
"""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch, Polygon, Wedge, Circle

import os  # v0.5 package: inputs resolved relative to this script (override with JBTS_INPUTS)
HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.environ.get("JBTS_INPUTS", os.path.join(HERE, "inputs"))
STYLE = os.path.join(HERE, "jbts_style.py")
MANU = os.path.join(INP, "Manuscript_JACC-BTS_v0.4_text.md")
exec(open(STYLE).read())
apply_jbts_style()
MS = open(MANU, encoding="utf-8").read().replace("\u2212", "-") if os.path.exists(MANU) else None  # not distributed in the public release

REG = []
def reg(text, lock, comp="locked v0.4 value"):
    ok = (lock.replace("\u2212", "-") in MS) if MS is not None else None
    if ok is False: print("LOCK MISMATCH", lock)
    text = text.replace("\n", " / ")
    REG.append(dict(figure="Central Illustration", panel="CI", printed_text=text, value=text, source_file="Manuscript v0.4 text",
                    computation=comp, matches_locked_text="yes", locked_text_phrase=lock,
                    check=("not checked (manuscript text not distributed)" if ok is None else ("verbatim in v0.4 text" if ok else "MISMATCH"))))

W, H = 183.0, 170.0
fig = plt.figure(figsize=(W * MM, H * MM))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.set_axis_off()
DARK, MID, LIGHT = "#4D4D4D", "#7F7F7F", "#F2F2F2"
GCI, GCR = GC_CLASS["induced"], GC_CLASS["repressed"]

def T(x, y, s, fs=9, **kw):
    kw.setdefault("va", "center"); kw.setdefault("ha", "left")
    return ax.text(x, y, s, fontsize=fs, linespacing=1.12, **kw)
def rbox(x, y, w, h, fc="white", ec="#BBBBBB", lw=0.7, r=1.5, ls="-", z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw, ls=ls, zorder=z))
def arrow(x0, y0, x1, y1, c=DARK, lw=1.0, ms=8):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=ms, lw=lw, color=c, zorder=3))
def band(y, h, num, ttl, tw=72):
    rbox(1, y, W - 2, h, fc="white", ec="#9A9A9A", lw=0.8, r=2)
    ax.add_patch(FancyBboxPatch((1, y), tw, 5.6, boxstyle="round,pad=0,rounding_size=2", fc=DARK, ec="none", zorder=2))
    ax.add_patch(Rectangle((1, y + 3), tw, 2.6, fc=DARK, ec="none", zorder=2))
    T(3.5, y + 2.9, f"{num}   {ttl}", fs=10, color="white", fontweight="bold", zorder=3)
def heart(cx, cy, s, fc, ec="#333333"):
    t = np.linspace(0, 2 * np.pi, 200)
    x = 16 * np.sin(t) ** 3; yy = 13 * np.cos(t) - 5 * np.cos(2 * t) - 2 * np.cos(3 * t) - np.cos(4 * t)
    ax.add_patch(Polygon(np.c_[cx + s * x / 17, cy - s * yy / 17], closed=True, fc=fc, ec=ec, lw=0.8, zorder=3))
def drop(cx, cy, s, fc):
    t = np.linspace(0, 2 * np.pi, 120)
    x = s * np.sin(t) * (1 - np.cos(t)) / 2 * 1.2; yy = -s * np.cos(t)
    ax.add_patch(Polygon(np.c_[cx + x, cy + yy], closed=True, fc=fc, ec="none", zorder=4))

# ---------------------------------------------------------------- title
T(2, 4.0, "A glucocorticoid-like imprint in donor control hearts can confound heart failure transcriptomics",
  fs=10.5, fontweight="bold")

# ---------------------------------------------------------------- band 1
y1 = 8.5; band(y1, 30.0, "1", "Standard design: donor control hearts versus failing hearts", 112)
heart(14, y1 + 17.5, 6.8, CONTROL)
T(3.0, y1 + 27.3, "Donor control heart", fs=9, ha="left", fontweight="bold")
for i, s in enumerate(["Organ donor, not transplanted", "Brain or circulatory death", "Drug exposure not recorded\nin any cohort"]):
    yy = [y1 + 10.2, y1 + 14.8, y1 + 20.6][i]
    ax.add_patch(Circle((24.5, yy - (1.6 if i == 2 else 0)), 0.7, fc=MID, ec="none", zorder=3))
    T(26.3, yy, s, fs=9)
reg("Drug exposure not recorded in any cohort", "not recorded in any cohort")
T(74.0, y1 + 17.5, "vs", fs=10, ha="center", style="italic", color=MID)
heart(89, y1 + 15.6, 7.2, "#6E6E6E")
T(89, y1 + 26.2, "Failing heart\n(explant or biopsy)", fs=9, ha="center", fontweight="bold", va="center")
arrow(101.5, y1 + 17.5, 110.0, y1 + 17.5, lw=1.2)
rbox(111, y1 + 7.5, 70, 20.5, fc=LIGHT, ec="none")
T(146, y1 + 11.3, "Case\u2013control transcriptome contrast", fs=9.5, ha="center", fontweight="bold")
for k, (hc, hf) in enumerate([(5.0, 2.0), (4.4, 1.7), (3.6, 1.6)]):
    xb = 115 + k * 6.0
    ax.add_patch(Rectangle((xb, y1 + 25.0 - hc), 2.2, hc, fc=CONTROL, ec="none", zorder=3))
    ax.add_patch(Rectangle((xb + 2.5, y1 + 25.0 - hf), 2.2, hf, fc="#6E6E6E", ec="none", zorder=3))
ax.plot([114, 133.5], [y1 + 25.0] * 2, color=DARK, lw=0.6, zorder=3)
T(158.0, y1 + 20.0, "A GC-like program in controls\ncan make its genes\nappear \u2018lost\u2019 in HF", fs=9, ha="center")

# ---------------------------------------------------------------- band 2
y2 = 41.5; band(y2, 80.0, "2", "What we found: 4 cohorts and an in vitro reference", 100)
rbox(3.5, y2 + 8.5, 30.5, 46.0, fc=LIGHT, ec="none")
T(18.75, y2 + 13.0, "In vitro GC\nreference", fs=9.5, ha="center", fontweight="bold")
T(18.75, y2 + 22.2, "4 human monocyte/\nmacrophage datasets", fs=8.5, ha="center")
for k, (val, col, lab) in enumerate([(394, GCI, "induced"), (508, GCR, "repressed")]):
    yy = y2 + 29.5 + k * 8.6
    ax.add_patch(Rectangle((5.5, yy - 1.4), 26.0 * val / 508, 2.8, fc=col, ec="none", zorder=3))
    T(5.5, yy + 3.8, f"{val} GC-{lab}", fs=9, color=col)
reg("394 GC-induced", "394 GC-induced"); reg("508 GC-repressed", "508 GC-repressed")
T(18.75, y2 + 50.2, "genes, applied\nto each cohort", fs=8, ha="center", color=MID, style="italic")
arrow(34.3, y2 + 31.0, 36.6, y2 + 31.0, lw=1.0)

cols = [("SCP3342", "snRNA-seq, HFpEF", "24/19", "19 vs 24"), ("Hahn", "bulk biopsy, HFpEF", "24/41", "41 vs 24"),
        ("MAGNet", "bulk LV, DCM", "102/128", "128 vs 102"), ("GSE290367", "snRNA-seq, DCM", "139/67", "67 vs 139")]
cx0, cw = 71.5, 27.4
mx = lambda j: cx0 + j * cw + cw / 2
for j, (nm, d, n, lk) in enumerate(cols):
    T(mx(j) + 1.6, y2 + 9.4, nm, fs=9, ha="center", fontweight="bold", color=DARK)
    _w = len(nm) * 1.0 + 1.6
    ax.plot(mx(j) + 1.6 - _w - 0.3, y2 + 9.4, marker=COHORT_MARKER[nm], ms=5.2, color=COHORT[nm], mec="#333333", mew=0.4, ls="none", zorder=4)
    T(mx(j), y2 + 13.0, d, fs=8, ha="center", color=DARK)
    T(mx(j), y2 + 16.4, n, fs=8, ha="center", color=DARK)
    reg(f"{nm} n={n} (ctrl/case)", lk, "control/case order; text gives cases vs controls")
T(mx(3), y2 + 19.6, "(45 ctrl in contrasts)", fs=7.5, ha="center", color=MID)
reg("45 ctrl in contrasts", "64 DCM and 45 controls", "GSE290367 controls from studies contributing both groups")
T(cx0 - 1.5, y2 + 16.4, "n, ctrl/case", fs=8, ha="right", color=MID, style="italic")
rows = [(y2 + 28.0, 6.4, "3-gene GC score higher\nin controls?\n(control-group SD)"),
        (y2 + 39.6, 5.2, "Explained by circadian\nor cAMP score?"),
        (y2 + 50.0, 5.2, "Disease separable from\nthe imprint? (post hoc)")]
for i, (yy, hh, rl) in enumerate(rows):
    ax.add_patch(Rectangle((37.0, yy - hh), cx0 + 4 * cw - 37.0, 2 * hh, fc=LIGHT if i % 2 == 0 else "white", ec="none", zorder=1))
    T(38.0, yy, rl, fs=8.5)
def varrow(x, y, c=GCI):
    arrow(x, y + 2.0, x, y - 2.2, c=c, lw=1.4, ms=7)
y_r1 = rows[0][0]
for j, (s, lk) in enumerate([("+3.4\u20136.5", "3.4 to 6.5"), ("+4.65", "+4.65 SD"), ("+1.67", "1.67 SD")]):
    varrow(mx(j), y_r1 - 3.0)
    T(mx(j), y_r1 + 1.5, s, fs=8.5, ha="center")
    reg(s, lk)
T(mx(1), y_r1 + 4.9, "(post hoc)", fs=7.5, ha="center", color=MID, style="italic")
reg("+4.65 (post hoc)", "(+4.65 SD; 95% CI: 4.03 to 5.27; post hoc)")
T(mx(3), y_r1, "Chaffin \u2191, Koenig =,\nReichart \u2193 (by study)", fs=8, ha="center")
reg("Chaffin up, Koenig =, Reichart down (by study)", "opposite study-level effects")
for j, g in enumerate(["No*", "\u2013", "No", "n.i."]):
    T(mx(j), rows[1][0] - (1.4 if j == 1 else 0), g, fs=9.5, ha="center", fontweight="bold")
T(mx(1), rows[1][0] + 2.5, "not tested", fs=7.5, ha="center", color=MID)
for j, g in enumerate(["No", "No", "Partly", "Yes (VIF only)"]):
    T(mx(j), rows[2][0], g, fs=9.5, ha="center", fontweight="bold")
reg("No* (SCP3342 circadian)", "only weakly identifiable", "footnote caveat")
T(37.5, y2 + 57.6, "*Circadian adjustment weakly identifiable in SCP3342 (the circadian score also separates groups).",
  fs=8, color=DARK, style="italic")
cy = y2 + 61.0
chips = [("HFpEF myeloid cells (discovery contrast;\nhypothesis-generating): 60% of the top\n100 downregulated genes are GC-induced\n(18% in bulk DCM)",
          ["60% of the top 100", "18% in bulk DCM", "(hypothesis-generating)"]),
         ("8 of 9 genes of a published LVAD\n\u2018GC receptor signature\u2019 are\nGC-induced in vitro", ["Eight of 9"]),
         ("Macrophages: GC-independent loss not\nshown (G3 failed); GC-induced loss\nattenuated 37%/34%, below the\npre-specified 50% (G4b not met)",
          ["(G3) failed", "attenuation of 37% and 34%", "pre-specified 50% threshold", "prediction G4b was not met"])]
for k, (s, lks) in enumerate(chips):
    x0 = 3.5 + k * 59.0
    rbox(x0, cy, 57.0, 17.0, fc="white", ec=GCI, lw=1.0)
    T(x0 + 28.5, cy + 8.5, s, fs=8, ha="center")
    for lk in lks: reg(s.replace("\n", " "), lk)

# ---------------------------------------------------------------- band 3
y3 = 124.5; band(y3, 30.0, "3", "What to do", 66)
todo = ["Record donor and recipient\ncorticosteroid exposure\n(agent, dose, timing)",
        "Compute a GC score\nin control and case\nsamples",
        "Check AUC and variance\ninflation before\nadjusting for the score",
        "Flag GC-responsive\ngenes (FKBP5, ZBTB16,\nPDK4, CD163)"]
for k, s in enumerate(todo):
    x0 = 3.5 + k * 44.5
    ax.add_patch(Circle((x0 + 3.2, y3 + 15.5), 3.0, fc=DARK, ec="none", zorder=3))
    T(x0 + 3.2, y3 + 15.6, str(k + 1), fs=10, ha="center", color="white", fontweight="bold", zorder=4)
    T(x0 + 7.3, y3 + 15.5, s, fs=8.5)
reg("FKBP5, ZBTB16, PDK4", "FKBP5, ZBTB16 and PDK4", "named GC-target genes (Highlights)")
T(91.5, y3 + 26.4, "GC-like: the pre-specified systemic-exposure test failed; donor drug exposure was not recorded in any cohort.",
  fs=8.5, ha="center", style="italic", color=DARK)
reg("systemic-exposure test failed", "systemic-exposure test failed")
T(1.5, H - 6.3, "AUC = area under the curve; cAMP = cyclic adenosine monophosphate; ctrl = control; DCM = dilated cardiomyopathy;\n"
  "GC = glucocorticoid; HFpEF = heart failure with preserved ejection fraction; LV = left ventricle; LVAD = left ventricular assist device;\n"
  "n.i. = non-informative; snRNA-seq = single-nucleus RNA sequencing; VIF = variance inflation factor.", fs=8, color=MID)

# ---- verify text extents
fig.canvas.draw(); _r = fig.canvas.get_renderer()
import matplotlib as _m
_t = [(t, t.get_window_extent(_r)) for t in fig.findobj(_m.text.Text) if t.get_text().strip()]
_ov = [(a.get_text()[:25], b.get_text()[:25]) for i, (a, ba) in enumerate(_t) for b, bb in _t[i + 1:] if ba.overlaps(bb)]
_out = [t.get_text()[:25] for t, bb in _t if bb.x0 < 0 or bb.x1 > fig.bbox.x1 or bb.y0 < 0 or bb.y1 > fig.bbox.y1]
print("OVERLAPS", _ov); print("OUTSIDE", _out)
matplotlib.rcParams["savefig.pad_inches"] = 0.0
save_figure(fig, "Central_Illustration")
pd.DataFrame(REG).to_csv("Central_Illustration_number_registry.csv", index=False)
print("CI saved; registry", len(REG), "mismatch", sum(r["check"] != "verbatim in v0.4 text" for r in REG))
