"""Central Illustration / Visual Abstract (storyboard section 6), v0.6: every text element >= 10 pt.
Stylised: no figure panel is duplicated. Every printed number is a locked v0.4 manuscript value;
the script verifies each phrase against the text.
Run: python ci.py  -> Central_Illustration.{pdf,tif,png}, Central_Illustration_number_registry.csv
"""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch, Polygon, Circle

import os  # public package: inputs resolved relative to this script (override with JBTS_INPUTS)
HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.environ.get("JBTS_INPUTS", os.path.join(HERE, "..", "figures_v05", "inputs"))
STYLE = os.path.join(HERE, "jbts_style.py")
MANU = os.path.join(INP, "Manuscript_JACC-BTS_v0.4_text.md")
exec(open(STYLE).read())
apply_jbts_style()
matplotlib.rcParams["font.family"] = [matplotlib.rcParams["font.sans-serif"][0], "FreeSans"]   # glyph fallback
MS = open(MANU, encoding="utf-8").read().replace("\u2212", "-") if os.path.exists(MANU) else None  # not distributed in the public release
FS, FTT = 10, SIZES["panel_title"]

REG = []
def reg(text, lock, comp="locked v0.4 value"):
    ok = (lock.replace("\u2212", "-") in MS) if MS is not None else None
    if ok is False: print("LOCK MISMATCH", lock)
    text = text.replace("\n", " / ")
    REG.append(dict(figure="Central Illustration", panel="CI", printed_text=text, value=text, source_file="Manuscript v0.4 text",
                    computation=comp, matches_locked_text="yes", locked_text_phrase=lock,
                    check=("not checked (manuscript text not distributed)" if ok is None else ("verbatim in v0.4 text" if ok else "MISMATCH"))))

W, H = 183.0, 214.0
fig = plt.figure(figsize=(W * MM, H * MM))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.set_axis_off()
DARK, MID, LIGHT = "#4D4D4D", "#6F6F6F", "#F2F2F2"
GCI, GCR = GC_CLASS["induced"], GC_CLASS["repressed"]

def T(x, y, s, fs=FS, **kw):
    kw.setdefault("va", "center"); kw.setdefault("ha", "left")
    return ax.text(x, y, s, fontsize=fs, linespacing=1.12, **kw)
def rbox(x, y, w, h, fc="white", ec="#BBBBBB", lw=0.7, r=1.5, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw, zorder=z))
def arrow(x0, y0, x1, y1, c=DARK, lw=1.0, ms=8):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=ms, lw=lw, color=c, zorder=3))
def band(y, h, num, ttl, tw):
    rbox(1, y, W - 2, h, fc="white", ec="#9A9A9A", lw=0.8, r=2)
    ax.add_patch(FancyBboxPatch((1, y), tw, 6.4, boxstyle="round,pad=0,rounding_size=2", fc=DARK, ec="none", zorder=2))
    ax.add_patch(Rectangle((1, y + 3.4), tw, 3.0, fc=DARK, ec="none", zorder=2))
    T(3.5, y + 3.3, f"{num}   {ttl}", fs=FTT, color="white", fontweight="bold", zorder=3)
def heart(cx, cy, s, fc, ec="#333333"):
    t = np.linspace(0, 2 * np.pi, 200)
    x = 16 * np.sin(t) ** 3; yy = 13 * np.cos(t) - 5 * np.cos(2 * t) - 2 * np.cos(3 * t) - np.cos(4 * t)
    ax.add_patch(Polygon(np.c_[cx + s * x / 17, cy - s * yy / 17], closed=True, fc=fc, ec=ec, lw=0.8, zorder=3))

# ---------------------------------------------------------------- title
T(2, 6.0, "A glucocorticoid-like imprint in donor control hearts can confound\nheart failure transcriptomics",
  fs=FTT, fontweight="bold")

# ---------------------------------------------------------------- band 1
y1 = 13.0; band(y1, 35.0, "1", "Standard design: donor control hearts versus failing hearts", 132)
heart(13, y1 + 19.5, 7.0, CONTROL)
T(2.5, y1 + 31.0, "Donor control heart", fontweight="bold")
for i, s in enumerate(["Organ donor, not transplanted", "Brain or circulatory death", "Drug exposure not recorded\nin any cohort"]):
    yy = [y1 + 11.0, y1 + 16.0, y1 + 22.8][i]
    ax.add_patch(Circle((25.0, yy - (2.0 if i == 2 else 0)), 0.75, fc=MID, ec="none", zorder=3))
    T(27.0, yy, s)
reg("Drug exposure not recorded in any cohort", "not recorded in any cohort")
T(91.0, y1 + 19.0, "vs", fs=FTT, ha="center", style="italic", color=MID)
heart(106, y1 + 18.0, 7.0, "#6E6E6E")
T(106, y1 + 29.6, "Failing heart\n(explant or biopsy)", ha="center", fontweight="bold")
arrow(115.0, y1 + 18.0, 121.0, y1 + 18.0, lw=1.2)
rbox(122, y1 + 8.5, 59, 25.0, fc=LIGHT, ec="none")
T(151.5, y1 + 13.6, "Case\u2013control\ntranscriptome contrast", ha="center", fontweight="bold")
for k, (hc, hf) in enumerate([(5.0, 2.0), (4.4, 1.7)]):
    xb = 124.5 + k * 6.0
    ax.add_patch(Rectangle((xb, y1 + 30.5 - hc), 2.2, hc, fc=CONTROL, ec="none", zorder=3))
    ax.add_patch(Rectangle((xb + 2.5, y1 + 30.5 - hf), 2.2, hf, fc="#6E6E6E", ec="none", zorder=3))
ax.plot([123.5, 137.0], [y1 + 30.5] * 2, color=DARK, lw=0.6, zorder=3)
T(159.5, y1 + 26.0, "A GC-like program\nin controls can make\nits genes appear\n\u2018lost\u2019 in HF", ha="center")

# ---------------------------------------------------------------- band 2
y2 = 51.0; band(y2, 103.0, "2", "What we found: 4 cohorts and an in vitro reference", 112)
rbox(3.0, y2 + 9.5, 33.5, 52.0, fc=LIGHT, ec="none")
T(19.75, y2 + 14.2, "In vitro GC\nreference", ha="center", fontweight="bold")
T(19.75, y2 + 26.6, "4 human\nmonocyte/\nmacrophage\ndatasets", ha="center")
for k, (val, col, lab) in enumerate([(394, GCI, "induced"), (508, GCR, "repressed")]):
    yy = y2 + 37.6 + k * 8.4
    ax.add_patch(Rectangle((4.5, yy - 1.3), 30.0 * val / 508, 2.6, fc=col, ec="none", zorder=3))
    T(4.5, yy + 3.7, f"{val} GC-{lab}", color=col)
reg("394 GC-induced", "394 GC-induced"); reg("508 GC-repressed", "508 GC-repressed")
T(19.75, y2 + 56.6, "genes, applied\nto each cohort", ha="center", color=MID, style="italic")
arrow(36.8, y2 + 38.0, 39.0, y2 + 38.0, lw=1.0)

cols = [("SCP3342", "snRNA-seq\nHFpEF", "24/19", "19 vs 24"), ("Hahn", "bulk biopsy\nHFpEF", "24/41", "41 vs 24"),
        ("MAGNet", "bulk LV\nDCM", "102/128", "128 vs 102"), ("GSE290367", "snRNA-seq\nDCM", "139/67\u2021", "67 vs 139")]
cx0, cw = 79.5, 25.4
mx = lambda j: cx0 + j * cw + cw / 2
for j, (nm, d, n, lk) in enumerate(cols):
    T(mx(j) + 1.8, y2 + 10.4, nm, ha="center", fontweight="bold", color=DARK)
    _w = len(nm) * 1.06 + 1.8
    ax.plot(mx(j) + 1.8 - _w - 0.4, y2 + 10.4, marker=COHORT_MARKER[nm], ms=5.4, color=COHORT[nm], mec="#333333", mew=0.4, ls="none", zorder=4)
    T(mx(j), y2 + 16.4, d, ha="center", color=DARK)
    T(mx(j), y2 + 22.4, n, ha="center", color=DARK)
    reg(f"{nm} n={n} (ctrl/case)", lk, "control/case order; text gives cases vs controls")
T(cx0 - 1.5, y2 + 22.4, "n, ctrl/case", ha="right", color=MID, style="italic")
reg("\u202145 controls in contrasts", "64 DCM and 45 controls", "GSE290367 controls from studies contributing both groups")
rows = [(y2 + 32.3, 7.4, "3-gene GC score\nhigher in controls?\n(ctrl-group SD)"),
        (y2 + 44.4, 4.7, "Explained by\ncircadian or cAMP?"),
        (y2 + 54.2, 5.1, "Disease separable?\n(post hoc)")]
for i, (yy, hh, rl) in enumerate(rows):
    ax.add_patch(Rectangle((39.5, yy - hh), cx0 + 4 * cw - 39.5, 2 * hh, fc=LIGHT if i % 2 == 0 else "white", ec="none", zorder=1))
    T(40.5, yy, rl)
y_r1 = rows[0][0]
for j, (s, lk) in enumerate([("+3.4\u20136.5", "3.4 to 6.5"), ("+4.65", "+4.65 SD"), ("+1.67", "1.67 SD")]):
    arrow(mx(j), y_r1 - 1.0, mx(j), y_r1 - 5.6, c=GCI, lw=1.4, ms=7)
    T(mx(j), y_r1 + 1.3, s, ha="center")
    reg(s, lk)
T(mx(1), y_r1 + 5.5, "(post hoc)", ha="center", color=MID, style="italic")
reg("+4.65 (post hoc)", "(+4.65 SD; 95% CI: 4.03 to 5.27; post hoc)")
T(mx(3), y_r1, "Chaffin \u2191\nKoenig =\nReichart \u2193", ha="center")
reg("Chaffin up, Koenig =, Reichart down (by study)", "opposite study-level effects")
for j, g in enumerate(["No*", "\u2013", "No", "n.i."]):
    T(mx(j), rows[1][0] - (1.8 if j == 1 else 0), g, fs=FTT, ha="center", fontweight="bold")
T(mx(1), rows[1][0] + 2.2, "not tested", ha="center", color=MID)
reg("No* (SCP3342 circadian)", "only weakly identifiable", "footnote caveat")
for j, g in enumerate(["No", "No", "Partly", "Yes"]):
    T(mx(j), rows[2][0] - (1.8 if j == 3 else 0), g, fs=FTT, ha="center", fontweight="bold")
T(mx(3), rows[2][0] + 2.3, "(VIF only)", ha="center")
T(39.5, y2 + 63.4, "*Circadian adjustment weakly identifiable in SCP3342 (the circadian score\nalso separates groups). \u202145 controls in contrasts (studies with both groups).",
  color=DARK, style="italic")
chips = [("HFpEF myeloid cells (discovery contrast; hypothesis-generating): 60% of the top 100\ndownregulated genes are GC-induced (18% in bulk DCM)",
          ["60% of the top 100", "18% in bulk DCM", "(hypothesis-generating)"], 9.6),
         ("8 of 9 genes of a published LVAD \u2018GC receptor signature\u2019 are GC-induced in vitro", ["Eight of 9"], 6.0),
         ("Macrophages: GC-independent loss not shown (G3 failed); GC-induced loss attenuated\n37%/34%, below the pre-specified 50% (G4b not met)",
          ["(G3) failed", "attenuation of 37% and 34%", "pre-specified 50% threshold", "prediction G4b was not met"], 9.6)]
cy = y2 + 69.0
for s, lks, hh in chips:
    rbox(3.0, cy, 177.0, hh, fc="white", ec=GCI, lw=1.0)
    T(91.5, cy + hh / 2, s, ha="center")
    for lk in lks: reg(s.replace("\n", " "), lk)
    cy += hh + 1.8

# ---------------------------------------------------------------- band 3
y3 = 157.0; band(y3, 39.0, "3", "What to do", 40)
todo = ["Record donor and recipient corticosteroid\nexposure (agent, dose, timing)",
        "Compute a GC score in control\nand case samples",
        "Check AUC and variance inflation\nbefore adjusting for the score",
        "Flag GC-responsive genes\n(FKBP5, ZBTB16, PDK4, CD163)"]
for k, s in enumerate(todo):
    x0 = 3.5 + (k % 2) * 89.0; yc = y3 + 13.0 + (k // 2) * 10.5
    ax.add_patch(Circle((x0 + 3.2, yc), 3.0, fc=DARK, ec="none", zorder=3))
    T(x0 + 3.2, yc + 0.1, str(k + 1), fs=FTT, ha="center", color="white", fontweight="bold", zorder=4)
    T(x0 + 7.5, yc, s)
reg("FKBP5, ZBTB16, PDK4", "FKBP5, ZBTB16 and PDK4", "named GC-target genes (Highlights)")
T(91.5, y3 + 33.4, "GC-like: the pre-specified systemic-exposure test failed;\ndonor drug exposure was not recorded in any cohort.",
  ha="center", style="italic", color=DARK)
reg("systemic-exposure test failed", "systemic-exposure test failed")
T(1.5, H - 8.6, "AUC = area under the curve; cAMP = cyclic adenosine monophosphate; ctrl = control;\n"
  "DCM = dilated cardiomyopathy; GC = glucocorticoid; HFpEF = heart failure with preserved ejection\n"
  "fraction; LV = left ventricle; LVAD = left ventricular assist device; n.i. = non-informative;\n"
  "snRNA-seq = single-nucleus RNA sequencing; VIF = variance inflation factor.", color=MID)

# ---- verify text extents and fonts
fig.canvas.draw(); _r = fig.canvas.get_renderer()
import matplotlib as _m
_t = [(t, t.get_window_extent(_r)) for t in fig.findobj(_m.text.Text) if t.get_text().strip()]
_ov = [(a.get_text()[:25], b.get_text()[:25]) for i, (a, ba) in enumerate(_t) for b, bb in _t[i + 1:] if ba.overlaps(bb)]
_out = [t.get_text()[:25] for t, bb in _t if bb.x0 < 0 or bb.x1 > fig.bbox.x1 or bb.y0 < 0 or bb.y1 > fig.bbox.y1]
print("OVERLAPS", _ov); print("OUTSIDE", _out)
_small = min_font_audit(fig, 10); print("MIN FONT AUDIT", _small); assert _small == []
save_figure(fig, "Central_Illustration")
pd.DataFrame(REG).to_csv("Central_Illustration_number_registry.csv", index=False)
print("CI saved; registry", len(REG), "mismatch", sum(r["check"] != "verbatim in v0.4 text" for r in REG))
