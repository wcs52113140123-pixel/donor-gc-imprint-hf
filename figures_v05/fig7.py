"""Figure 7 (R6 + R7): macrophage programs after GC adjustment, GTEx ventilator proxy, aging.
Standalone: reproduces Figure_7.{pdf,tif,png}, Figure_7_source_data.xlsx and Figure_7_number_registry.csv
from the input artifacts below (read-only).  Run: python fig7.py
"""
import os, re, tarfile, tempfile
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch, Polygon
from matplotlib.lines import Line2D

import os  # v0.5 package: inputs resolved relative to this script (override with JBTS_INPUTS)
HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.environ.get("JBTS_INPUTS", os.path.join(HERE, "inputs"))
STYLE = os.path.join(HERE, "jbts_style.py")
TAR_FI = os.path.join(INP, "figure_inputs_v05.tar.gz")
TAR_PD = os.path.join(INP, "figure_inputs_v05_perdonor.tar.gz")
SUPP = os.path.join(INP, "Supplemental_Tables_S1-S14.xlsx")
MANU = os.path.join(INP, "Manuscript_JACC-BTS_v0.4_text.md")

exec(open(STYLE).read())
apply_jbts_style()

# ---------------------------------------------------------------- data
TMP = tempfile.mkdtemp(prefix="fig7_")
for t in (TAR_FI, TAR_PD):
    with tarfile.open(t) as tf:
        tf.extractall(TMP)
PD_ = f"{TMP}/figure_inputs_v05_perdonor/v12_rerun/"
FI_ = f"{TMP}/figure_inputs_v05/04_results/"
SRC = {
    "gse": "PD:v12_rerun/GSE290367_Resident_DCM_scores.csv",
    "mag": "PD:v12_rerun/MAGNet_exclOverlap_DCM_scores.csv",
    "lv": "PD:v12_rerun/perdonor_scores_GTEx_LV.csv",
    "s7b": "XL:S7b_G2_G5",
    "r2": "FI:04_results/validation_v1.3/revision_posthoc/rev_posthoc_models.csv",
    "s11a": "XL:S11a_substates_dedup",
    "sub": "FI:04_results/validation_v1.3/substates_dedup/GSE290367_subset_donor_scores.csv",
    "s3": "XL:S3_gene_sets",
}
gse = pd.read_csv(PD_ + "GSE290367_Resident_DCM_scores.csv")
mag = pd.read_csv(PD_ + "MAGNet_exclOverlap_DCM_scores.csv")
lv = pd.read_csv(PD_ + "perdonor_scores_GTEx_LV.csv")
lv = lv[lv["analysis_ok"] == True].copy()
s7b = pd.read_excel(SUPP, "S7b_G2_G5")
s11a = pd.read_excel(SUPP, "S11a_substates_dedup")
s3 = pd.read_excel(SUPP, "S3_gene_sets")
r2 = pd.read_csv(FI_ + "validation_v1.3/revision_posthoc/rev_posthoc_models.csv")
sub = pd.read_csv(FI_ + "validation_v1.3/substates_dedup/GSE290367_subset_donor_scores.csv")
manuscript = open(MANU, encoding="utf-8").read() if os.path.exists(MANU) else None  # manuscript text is not distributed in the public release

assert (gse.case == 1).sum() == 64 and (gse.case == 0).sum() == 45
assert (mag.case == 1).sum() == 128 and (mag.case == 0).sum() == 102
assert len(lv) == 450 and (lv.vent == 1).sum() == 253

def s7(dataset, contrast, module, adjustment=None, test=None):
    q = (s7b.dataset == dataset) & (s7b.contrast == contrast) & (s7b.module == module)
    if adjustment is not None: q &= s7b.adjustment == adjustment
    if test is not None: q &= s7b.test == test
    r = s7b[q]; assert len(r) == 1, (dataset, contrast, module, adjustment, len(r)); return r.iloc[0]

def r2row(dataset, module, adjustment):
    r = r2[(r2.test == "R2_macrophage_GC3adj") & (r2.dataset == dataset) & (r2.module == module) & (r2.adjustment == adjustment)]
    assert len(r) == 1; return r.iloc[0]

# ---------------------------------------------------------------- registry / formatting
FAM = matplotlib.rcParams["font.sans-serif"][0]
matplotlib.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": FAM, "mathtext.it": f"{FAM}:italic",
                            "mathtext.bf": f"{FAM}:bold", "mathtext.default": "regular"})
REG, SHEETS = [], {}
MS = manuscript
def _n(s): return s.replace("\u2212", "-")
def reg(panel, text, value, src, comp, status, lock=None, where="figure"):
    """status: yes (locked v0.4 text; `lock` = verbatim phrase that must occur in the text), S (supplement), new."""
    note = ""
    if status == "yes":
        assert lock is not None, text
        if MS is None:
            note = "not checked (manuscript text not distributed)"
        else:
            ok = _n(lock) in _n(MS)
            note = "verbatim in v0.4 text" if ok else "MISMATCH: phrase not found in v0.4 text"
            if not ok: print("LOCK MISMATCH", panel, lock)
    REG.append(dict(figure="Figure 7", panel=panel, location=where, printed_text=text, value=value, source_file=src,
                    computation=comp, matches_locked_text=status, locked_text_phrase=lock or "", check=note))

def f(x, nd=2, sign=True):
    s = f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"
    return s.replace("-", "\u2212")
def fp(p):
    """2 significant figures as in the text; mathtext power of ten below 1e-4."""
    if p < 1e-4:
        m, e = f"{p:.1e}".split("e"); return f"{m}$\\times$10$^{{{int(e)}}}$"
    return f"{p:.2g}"
def ptext(ax, x, y, base, expo, fs=9, transform=None, **kw):
    """Draw base text then a separate superscript (>=7 pt) anchored to its right edge."""
    t = ax.text(x, y, base, fontsize=fs, transform=transform if transform is not None else ax.transData, **kw)
    ax.annotate(expo, xy=(1.0, 0.28), xycoords=t, fontsize=7.2, ha="left", va="bottom", color=kw.get("color", TEXT))
    return t
def psplit(p):
    m, e = f"{p:.1e}".split("e"); return f"{m}\u00d710", str(int(e)).replace("-", "\u2212")
def plain(s):  # registry form of a mathtext string
    return re.sub(r"\$\\times\$10\$\^\{(-?\d+)\}\$", lambda m: "\u00d710^" + m.group(1), s)

# ---------------------------------------------------------------- canvas helpers (mm, top-left origin)
W, H = 183.0, 230.0
fig = plt.figure(figsize=(W * MM, H * MM))
def AX(x, y, w, h, **kw):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **kw)
def FT(x, y, s, **kw):
    kw.setdefault("fontsize", 9); kw.setdefault("va", "top"); kw.setdefault("ha", "left")
    return fig.text(x / W, 1 - y / H, s, **kw)
def letter(x, y, L): FT(x, y, L, fontsize=SIZES["panel_letter"], fontweight="bold")
def title(x, y, t): FT(x, y + 0.9, t, fontsize=SIZES["panel_title"], fontweight="bold", linespacing=1.05)
def banner(y, t):
    a = AX(0, y, W, 4.2); a.set_axis_off()
    a.add_patch(Rectangle((0, 0), 1, 1, color="#666666", lw=0, transform=a.transAxes))
    a.text(0.006, 0.5, t, color="white", fontsize=9, va="center", ha="left", fontweight="bold", transform=a.transAxes)
def tag(x, y, s="post hoc"): FT(x, y, s, ha="right", **POSTHOC_TAG)

GSEc, MAGc, GTXc = COHORT["GSE290367"], COHORT["MAGNet"], COHORT["GTEx"]
GSEm, MAGm = COHORT_MARKER["GSE290367"], COHORT_MARKER["MAGNet"]
MINUS = "\u2212"
def mticks(ax, ticks, axis="x", nd=None):
    lab = [(f"{t:g}").replace("-", MINUS) for t in ticks]
    (ax.set_xticks if axis == "x" else ax.set_yticks)(ticks); (ax.set_xticklabels if axis == "x" else ax.set_yticklabels)(lab)
rng = np.random.default_rng(7)   # jitter only (display)

# schematic primitives
def box(ax, x, y, w, h, txt, fc="white", ec="#555555", fs=8, tc=TEXT, lw=0.7, ls="-"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2", fc=fc, ec=ec, lw=lw, ls=ls))
    ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs, color=tc, linespacing=1.12)
def arrow(ax, x0, y0, x1, y1, c="#555555"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=6, lw=0.7, color=c))
def pill(ax, x, y, w, h, txt, kind, fs=8):
    """Test badges (paper-wide): passed = black; failed = white, outlined, x; notmet = half-filled; noteval = dashed."""
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={h/2}", fc="black" if kind == "passed" else "white",
                       ec="black", lw=0.8, ls="--" if kind == "noteval" else "-")
    ax.add_patch(p)
    if kind == "notmet":
        hb = Rectangle((x, y), w * 0.2, h, fc="black", ec="none"); ax.add_patch(hb); hb.set_clip_path(p)
        ax.text(x + w * 0.6, y + h / 2, txt, ha="center", va="center", fontsize=fs)
    elif kind == "failed":
        ax.text(x + 2.4, y + h / 2, "\u00d7", ha="center", va="center", fontsize=fs + 1, fontweight="bold")
        ax.text(x + w / 2 + 1.4, y + h / 2, txt, ha="center", va="center", fontsize=fs)
    else:
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs, color="white" if kind == "passed" else TEXT)

# ================================================================= banners
banner(0, "Macrophage modules after GC adjustment (GSE290367 resident macrophages; MAGNet)")
banner(112.0, "GTEx: ventilator-managed vs other deaths (supportive test G2, not blinded)")
banner(170.5, "Aging (pre-specified test G5)")

# ================================================================= A schematic
letter(0, 5.5, "A"); title(6, 5.5, "What each module tests")
a = AX(0, 12.0, 52, 47.0); a.set_xlim(0, 52); a.set_ylim(47.0, 0); a.set_axis_off()
n_g3 = int((s3.gene_set == "G3_GC_insensitive").sum()); n_g4 = int((s3.gene_set == "G4_GC_induced").sum())
assert (n_g3, n_g4) == (4, 14)
box(a, 0.5, 0.3, 51, 5.0, "Is the DCM loss GC-independent?", fc="#F0F0F0", ec="#F0F0F0")
arrow(a, 13, 5.5, 13, 7.8); arrow(a, 39, 5.5, 39, 7.8)
box(a, 0.5, 8.0, 25, 13.8, f"G3 ({n_g3} genes)\nGC-insensitive\nITGAV, ITGB5,\nLRP1, SCARB1", ec=GC_CLASS["insensitive"], lw=1.5)
box(a, 26.5, 8.0, 25, 13.8, f"G4 ({n_g4} genes)\nGC-induced\nC1QA\u2013C, CD163,\nMERTK, LYVE1 \u2026", ec=GC_CLASS["induced"], lw=1.5)
arrow(a, 13, 22.0, 13, 24.3); arrow(a, 39, 22.0, 39, 24.3)
box(a, 0.5, 24.5, 25, 7.6, "Predicted: declines\nafter adjustment", ec="#BBBBBB")
box(a, 26.5, 24.5, 25, 7.6, "Predicted: decline\nattenuates \u226550%", ec="#BBBBBB")
pill(a, 1.5, 34.2, 23, 4.6, "G3 failed", "failed")
a.text(13, 40.0, "no decline in\neither cohort", ha="center", va="top", fontsize=8, linespacing=1.1)
pill(a, 27.5, 34.2, 23, 4.6, "G4b not met", "notmet")
a.text(39, 39.6, "37% and 34%", ha="center", va="top", fontsize=8)
pill(a, 26.7, 42.6, 24.6, 4.2, "G4a not evaluated", "noteval")
reg("A", "G3 failed", "failed", "Table 1 (v0.4)", "pre-specified outcome", "yes", lock="(G3) failed")
reg("A", "G4b not met", "not met", "Table 1 (v0.4)", "pre-specified outcome", "yes", lock="prediction G4b was not met")
reg("A", "37% and 34%", "37; 34", SRC["s7b"], "G4 attenuation (base - +GC)/base, GSE290367; MAGNet", "yes", lock="attenuation of 37% and 34%")
reg("A", "G4a not evaluated", "not evaluated", "Table 1 (v0.4)", "pre-specified outcome", "yes", lock="prediction G4a was not evaluated")
reg("A", f"G3 ({n_g3} genes)", n_g3, SRC["s3"], "count G3_GC_insensitive", "yes", lock="G3, GC-insensitive receptor set (ITGAV, ITGB5, LRP1, SCARB1)")
reg("A", f"G4 ({n_g4} genes)", n_g4, SRC["s3"], "count G4_GC_induced", "yes", lock="G4, GC-induced efferocytic subset (14 genes)")
SHEETS["A_gene_sets"] = s3[s3.gene_set.isin(["G3_GC_insensitive", "G4_GC_induced"])][["version", "gene_set", "gene", "meta_log2FC", "GC_class"]]

# ================================================================= B strip + box (2 x 2, shared y per row)
letter(57, 5.5, "B"); title(63, 5.5, "G3 and G4 per donor")
dsets = [("GSE290367", gse, GSEc, GSEm, "GSE290367\nresident mac."), ("MAGNet", mag, MAGc, MAGm, "MAGNet\nbulk LV")]
B_rows = []
for i, mod in enumerate(["G3", "G4"]):
    lo = min(gse[mod].min(), mag[mod].min()); hi = max(gse[mod].max(), mag[mod].max()); pad = 0.06 * (hi - lo)
    for j, (nm, d, col, mk, lab) in enumerate(dsets):
        ax = AX(70 + j * 21.5, 20.0 + i * 20.5, 19.5, 14.5)
        for k, cs in enumerate([0, 1]):
            v = d.loc[d.case == cs, mod].dropna().values
            ax.boxplot(v, positions=[k], widths=0.6, showfliers=False, patch_artist=True, zorder=1,
                       boxprops=dict(fc="none", ec="#333333", lw=0.7), medianprops=dict(color="#111111", lw=1.1),
                       whiskerprops=dict(lw=0.6, color="#333333"), capprops=dict(lw=0.6, color="#333333"))
            ax.scatter(k + rng.uniform(-0.22, 0.22, len(v)), v, s=4, marker=mk, color=CONTROL if cs == 0 else col, lw=0, alpha=0.85, zorder=2)
            B_rows += [dict(module=mod, cohort=nm, group="DCM" if cs else "control", score=x) for x in v]
        ax.set_xlim(-0.6, 1.6); ax.set_xticks([0, 1]); ax.set_xticklabels(["Ctrl", "DCM"] if i == 1 else ["", ""])
        ax.axhline(0, color="#7F7F7F", lw=0.6, ls="--", zorder=0)
        ax.set_ylim(lo - pad, hi + pad); mticks(ax, [-3, 0, 3], "y")
        if j == 0: ax.set_ylabel(f"{mod} (SD)", labelpad=1)
        else: ax.set_yticklabels([])
        if i == 0: ax.set_title(lab, fontsize=9, pad=2, color=col, linespacing=1.0)
SHEETS["B_G3_G4_per_donor"] = pd.DataFrame(B_rows)
reg("B", "45 controls vs 64 DCM (GSE290367); 102 vs 128 (MAGNet)", "45/64; 102/128", f"{SRC['gse']}; {SRC['mag']}", "row counts by case", "yes",
    lock="64 DCM and 45 controls", where="legend")

# ================================================================= C forest (module x cohort x adjustment)
letter(117, 5.5, "C"); title(123, 5.5, "Module effects without and\nwith GC adjustment")
cax = AX(120.5, 19.0, 36.5, 47.5)
mods = [("G3", f"G3, GC-insensitive ({n_g3})"), ("G4", f"G4, GC-induced ({n_g4})"),
        ("RES_EFF_core", "Resident efferocytic (9)"), ("MHCII_core", "MHC-II (10)")]
cohorts = [("GSE290367", "GSE290367_Resident", "DCM_vs_ND", "age+sex+study", "age+sex+study+GC", "GSE290367_Resident_DCM", GSEc, GSEm),
           ("MAGNet", "MAGNet_excl_overlap", "DCM_vs_NF", "age+sex", "age+sex+GC", "MAGNet_exclOverlap_DCM", MAGc, MAGm)]
locks = {("G3", "GSE290367"): "unadjusted estimates were −0.04 and +0.31 SD", ("G3", "MAGNet"): "+0.43 to +0.51 SD in MAGNet",
         ("G4", "GSE290367"): "from −1.18 to −0.75 SD after adjustment in GSE290367", ("G4", "MAGNet"): "from −0.97 to −0.65 SD in MAGNet",
         ("RES_EFF_core", "GSE290367"): "(−2.58 to −2.27 SD", ("RES_EFF_core", "MAGNet"): "MAGNet (−1.53 to −1.06 SD)",
         ("MHCII_core", "GSE290367"): "GSE290367 +0.94 to +0.80 SD", ("MHCII_core", "MAGNet"): "MAGNet +1.50 to +1.56 SD"}
C_rows = []; y = 0.0; TX = 1.04
for mod, mlab in mods:
    cax.text(0.0, y, mlab, transform=cax.get_yaxis_transform(), ha="left", va="center", fontsize=9, zorder=5,
             bbox=dict(fc="white", ec="none", pad=0.4))
    if mod == "RES_EFF_core":
        cax.text(1.0 + 26.0 / 36.5, y, "circular (GSE)", transform=cax.get_yaxis_transform(), ha="right", va="center", fontsize=8, **POSTHOC_TAG)
    for ci, (cn, ds, ct, adj0, adj1, r2ds, col, mk) in enumerate(cohorts):
        yy = y + 1.35 + ci * (2.25 if mod == "G4" else 1.3)
        b0 = s7(ds, ct, mod, adj0); b1 = s7(ds, ct, mod, adj1)
        for row, filled, dy in [(b0, True, -0.3), (b1, False, 0.0)]:
            cax.plot([row.ci_lo, row.ci_hi], [yy + dy] * 2, color=col, lw=0.8, zorder=1)
            cax.plot(row.beta, yy + dy, marker=mk, ms=4.4, mfc=col if filled else "white", mec=col, mew=0.9, zorder=3, ls="none")
            C_rows.append(dict(module=mod, cohort=cn, adjustment=row.adjustment, marker="filled" if filled else "open", pipeline="v1.2",
                               beta=row.beta, ci_lo=row.ci_lo, ci_hi=row.ci_hi, p=row.p, source=SRC["s7b"]))
        txt = f"{f(b0.beta)} \u2192 {f(b1.beta)}"
        cax.text(TX, yy - 0.15, txt, transform=cax.get_yaxis_transform(), ha="left", va="center", fontsize=9, color=col)
        reg("C", txt, f"{b0.beta:.4f} -> {b1.beta:.4f}", SRC["s7b"], f"{ds} {ct} {mod}: {adj0} -> {adj1} (DCM - control)", "yes",
            lock=locks[(mod, cn)])
        if mod == "G3" and cn == "GSE290367":
            reg("C", "", f"{b1.ci_lo:.2f} to {b1.ci_hi:.2f}", SRC["s7b"], "G3 +GC 95% CI, GSE290367", "yes",
                lock="−0.37 SD (95% CI: −0.93 to 0.19)", where="legend")
        if mod == "G4":
            h = r2row(r2ds, "G4", "+GC3+MAC")
            cax.plot([h.ci_lo, h.ci_hi], [yy + 0.3] * 2, color=col, lw=0.8, zorder=1)
            cax.plot(h.beta, yy + 0.3, marker=mk, ms=4.4, fillstyle="left", mfc=col, mfcalt="white", mec=col, mew=0.9, zorder=3, ls="none")
            C_rows.append(dict(module=mod, cohort=cn, adjustment="+3-gene GC score + macrophage content (post hoc)", marker="half",
                               pipeline="v1.3", beta=h.beta, ci_lo=h.ci_lo, ci_hi=h.ci_hi, p=h.p, source=SRC["r2"]))
            t2 = f"+mac. {f(h.beta)}"
            cax.text(TX, yy + 0.95, t2, transform=cax.get_yaxis_transform(), ha="left", va="center", fontsize=8, color=col)
            reg("C", t2, f"{h.beta:.4f}", SRC["r2"], f"R2_macrophage_GC3adj {r2ds} G4 +GC3+MAC", "yes",
                lock="+0.04 SD (95% CI: −0.51 to 0.58) in GSE290367" if cn == "GSE290367" else "−0.42 SD (P<10⁻⁴) in MAGNet")
    y += 5.8 if mod == "G4" else 4.2
cax.axvline(0, color="#7F7F7F", lw=0.6, ls="--", zorder=0)
cax.set_ylim(y - 1.05, -0.75); cax.set_yticks([]); cax.spines["left"].set_visible(False)
cax.set_xlim(-3.4, 2.1); mticks(cax, [-3, -2, -1, 0, 1, 2])
cax.set_xlabel(f"DCM {MINUS} control (SD, 95% CI)", labelpad=1)
FT(158.0, 15.5, "base \u2192 +GC", fontsize=9)
SHEETS["C_module_effects"] = pd.DataFrame(C_rows)
# key
kax = AX(118, 76.0, 65, 13.0); kax.set_axis_off(); kax.set_xlim(0, 65); kax.set_ylim(13.0, 0)
def kmk(x, yk, mk, col, fs="full", ms=4.4):
    kax.plot(x, yk, marker=mk, ms=ms, ls="none", mew=0.9, mec=col,
             mfc=("white" if fs == "none" else col), fillstyle=("left" if fs == "left" else "full"), mfcalt="white")
kmk(1.5, 1.6, GSEm, GSEc); kax.text(3.5, 1.6, "GSE290367", va="center", fontsize=9, color=GSEc)
kmk(23.0, 1.6, MAGm, MAGc, ms=4.0); kax.text(25.0, 1.6, "MAGNet", va="center", fontsize=9, color=MAGc)
kmk(1.5, 5.6, "o", "#555555"); kax.text(3.5, 5.6, "base", va="center", fontsize=9)
kmk(14.0, 5.6, "o", "#555555", "none"); kax.text(16.0, 5.6, "+ 8-gene GC response score", va="center", fontsize=9)
kmk(1.5, 9.6, "o", "#555555", "left"); kax.text(3.5, 9.6, "+ 3-gene GC + mac. content", va="center", fontsize=9)
kax.text(46.0, 9.6, "(post hoc)", va="center", ha="left", fontsize=9, **POSTHOC_TAG)
# attenuation inset (G4)
FT(118, 88.3, "G4 attenuation by GC adjustment", fontsize=9)
att = AX(137, 96.0, 24, 10.0)
att_rows = []
for k, (cn, ds, ct, adj0, adj1, r2ds, col, mk) in enumerate(cohorts):
    b0 = s7(ds, ct, "G4", adj0).beta; b1 = s7(ds, ct, "G4", adj1).beta; a8 = 100 * (b0 - b1) / b0
    r0 = r2row(r2ds, "G4", "base").beta; r3 = r2row(r2ds, "G4", "+GC3").beta; a3 = 100 * (r0 - r3) / r0
    for j, (val, fc) in enumerate([(a8, col), (a3, "white")]):
        att.barh(k * 1.0 + j * 0.42, val, height=0.36, color=fc, ec=col, lw=0.8, hatch="//////" if j == 1 else None)
    att.text(1.03, k * 1.0 + 0.21, f"{a8:.0f}% / {a3:.0f}%", transform=att.get_yaxis_transform(), va="center", fontsize=9, color=col)
    att.text(-1.5, k * 1.0 + 0.21, cn, ha="right", va="center", fontsize=9, color=col)
    att_rows += [dict(cohort=cn, score="8-gene GC response (pre-specified, v1.2)", base=b0, adjusted=b1, attenuation_pct=a8),
                 dict(cohort=cn, score="3-gene GC score (post hoc, v1.3)", base=r0, adjusted=r3, attenuation_pct=a3)]
    reg("C", f"{a8:.0f}%", round(a8, 2), SRC["s7b"], f"G4 attenuation 8-gene {cn}", "yes", lock="attenuation of 37% and 34%")
    reg("C", f"{a3:.0f}%", round(a3, 2), SRC["r2"], f"G4 attenuation 3-gene {cn} (post hoc, v1.3)", "yes", lock="similar attenuation (31% and 30%)")
att.axvline(50, color="#222222", lw=0.8, ls=(0, (2, 1.5)))
att.text(50.8, -0.75, "G4b", ha="left", va="center", fontsize=8)
att.set_xlim(0, 55); att.set_ylim(1.7, -1.05); att.set_yticks([]); att.spines["left"].set_visible(False)
att.set_xticks([0, 25, 50]); att.set_xticklabels(["0%", "25%", "50%"])
FT(118, 92.2, "solid, 8-gene; hatched, 3-gene score (post hoc)", fontsize=8)
SHEETS["C_G4_attenuation"] = pd.DataFrame(att_rows)

# ================================================================= D scatter GC vs G4
letter(0, 62.5, "D"); title(6, 62.5, "G4 and the 8-gene GC response score, per donor")
FT(6, 67.6, "Descriptive; G4a not evaluated. Gray, control.", fontsize=9, style="italic", color="#7F7F7F")
D_rows = []
for j, (nm, d, col, mk, lab) in enumerate(dsets):
    ax = AX(12 + j * 52, 76.0, 43, 25.0)
    for cs in [0, 1]:
        dd = d[d.case == cs]
        ax.scatter(dd.GC_response, dd.G4, s=7, marker=mk, color=CONTROL if cs == 0 else col, lw=0, alpha=0.85, zorder=2 if cs else 1)
        D_rows += [dict(cohort=nm, group="DCM" if cs else "control", GC_response=a_, G4=b_) for a_, b_ in zip(dd.GC_response, dd.G4)]
    if j == 0: ax.set_ylabel("G4 score (SD)", labelpad=1)
    ax.margins(0.05)
    mticks(ax, [-5, 0, 5] if j == 0 else [-6, -4, -2, 0, 2]); mticks(ax, [-4, -2, 0, 2], "y")
    FT(12 + j * 52, 72.0, lab.replace("\n", " "), color=col)
FT(59.5, 106.3, "8-gene GC response score (SD)", ha="center", fontsize=10)
SHEETS["D_GC_vs_G4"] = pd.DataFrame(D_rows)

# ================================================================= E forest GTEx
letter(0, 118.0, "E"); title(6, 118.0, "Ventilator-managed deaths: efferocytic up,\nMHC-II down, GC score itself not higher")
eax = AX(21.5, 135.0, 23.0, 25.0)
emods = [("GC_response", "8-gene GC\nresponse\nscore"), ("RES_EFF_core", "Resident\nefferocytic"), ("MHCII_core", "MHC-II")]
E_rows = []
elocks = {("LV", "GC_response"): "(+0.16 SD; P=0.34)", ("LV", "RES_EFF_core"): "LV (+0.30 SD; P=0.0021)",
          ("AA", "RES_EFF_core"): "atrial appendage (+0.14 SD; P=0.0077)", ("LV", "MHCII_core"): "(LV −0.36 SD; P=0.0056)"}
for i, (mod, lab) in enumerate(emods):
    eax.text(-0.06, i * 1.0, lab, transform=eax.get_yaxis_transform(), ha="right", va="center", fontsize=9, linespacing=1.0)
    for k, (ds, mk, tl) in enumerate([("GTEx_LV", "o", "LV"), ("GTEx_AA", "D", "AA")]):
        r = s7(ds, "ventilator_vs_other", mod, test="G2"); yk = i * 1.0 - 0.2 + 0.4 * k
        eax.plot([r.ci_lo, r.ci_hi], [yk, yk], color=GTXc, lw=0.8)
        eax.plot(r.beta, yk, marker=mk, ms=4.6 if mk == "o" else 3.9, mfc=GTXc, mec=GTXc, mew=0.9, ls="none")
        if r.p < 1e-4:
            pb, pe = psplit(r.p); txt = f"{f(r.beta)}, P={pb}^{pe}"
            ptext(eax, 1.04, yk, f"{f(r.beta)}, P={pb}", pe, transform=eax.get_yaxis_transform(), ha="left", va="center")
        else:
            txt = f"{f(r.beta)}, P={fp(r.p)}"
            eax.text(1.04, yk, txt, transform=eax.get_yaxis_transform(), ha="left", va="center", fontsize=9)
        key = (tl, mod)
        reg("E", plain(txt), f"{r.beta:.4f}; P={r.p:.3g}", SRC["s7b"], f"G2 {ds} ventilator_vs_other {mod}",
            "yes" if key in elocks else "S", lock=elocks.get(key))
        E_rows.append(dict(tissue=tl, module=mod, n=r.n, n_vent=r.n_case, beta=r.beta, ci_lo=r.ci_lo, ci_hi=r.ci_hi, p=r.p, adjustment=r.adjustment))
eax.axvline(0, color="#7F7F7F", lw=0.6, ls="--", zorder=0)
eax.set_ylim(2.5, -0.5); eax.set_yticks([]); eax.spines["left"].set_visible(False)
eax.set_xlim(-0.9, 0.6); mticks(eax, [-0.5, 0, 0.5])
eax.set_xlabel(f"Ventilator {MINUS} other (SD)", labelpad=1)
ek = AX(0, 128.0, 78, 4.5); ek.set_axis_off(); ek.set_xlim(0, 78); ek.set_ylim(0, 1)
ek.plot(14.5, 0.5, marker="o", ms=4.6, color=GTXc, ls="none"); ek.text(16.3, 0.5, "LV (n=197/253)", va="center", fontsize=9)
ek.plot(44.5, 0.5, marker="D", ms=3.9, color=GTXc, ls="none"); ek.text(46.3, 0.5, "AA (n=246/215)", va="center", fontsize=9)
reg("E", "LV (n=197/253)", "450 = 197 other / 253 ventilator", SRC["s7b"], "n GTEx_LV G2 (other/ventilator)", "yes", lock="LV n=450 (253 ventilator-managed)")
reg("E", "AA (n=246/215)", "461 = 246 other / 215 ventilator", SRC["s7b"], "n GTEx_AA G2 (other/ventilator)", "yes", lock="atrial appendage n=461 (215 ventilator-managed)")
SHEETS["E_GTEx_ventilator"] = pd.DataFrame(E_rows)

# ================================================================= F violins (raw, unadjusted)
letter(80, 118.0, "F"); title(86, 118.0, "Per-sample LV scores\n(unadjusted)")
fax = AX(91, 135.0, 33, 25.0)
F_rows = []
for j, col_ in enumerate(["GC_response", "RES_EFF_core"]):
    for k, vv in enumerate([0, 1]):
        pos = j * 2.5 + k
        v = lv.loc[lv.vent == vv, col_].values; c = CONTROL if vv == 0 else GTXc
        vp = fax.violinplot(v, positions=[pos], widths=0.85, showextrema=False)
        for b in vp["bodies"]: b.set_facecolor(c); b.set_edgecolor("none"); b.set_alpha(0.35)
        fax.scatter(pos + rng.uniform(-0.13, 0.13, len(v)), v, s=1.5, color=c, lw=0, alpha=0.8, zorder=2)
        fax.plot([pos - 0.28, pos + 0.28], [np.median(v)] * 2, color="#111111", lw=1.2, zorder=3)
        F_rows += [dict(score=col_, group="ventilator" if vv else "other", value=x_) for x_ in v]
fax.set_xticks([0.5, 3.0]); fax.set_xticklabels(["8-gene GC\nresponse", "Resident\nefferocytic"], linespacing=1.0)
fax.set_xlim(-0.6, 4.1); fax.set_ylim(-3.3, 3.6); mticks(fax, [-2, 0, 2], "y")
fax.axhline(0, color="#7F7F7F", lw=0.6, ls="--", zorder=0)
fax.set_ylabel("Score (SD)", labelpad=1)
fk = AX(84, 128.0, 42, 4.5); fk.set_axis_off(); fk.set_xlim(0, 42); fk.set_ylim(0, 1)
fk.plot(1.5, 0.5, "o", ms=4, color=CONTROL); fk.text(3.0, 0.5, "other (197)", va="center", fontsize=9)
fk.plot(22.0, 0.5, "o", ms=4, color=GTXc); fk.text(23.5, 0.5, "ventilator (253)", va="center", fontsize=9)
reg("F", "other (197); ventilator (253)", "197; 253", SRC["lv"], "counts analysis_ok by vent", "yes", lock="LV n=450 (253 ventilator-managed)")
SHEETS["F_GTEx_LV_per_sample"] = pd.DataFrame(F_rows)

# ================================================================= G scatter GC vs RES (GTEx LV)
letter(129, 118.0, "G"); title(135, 118.0, "Resident efferocytic score\nis higher with higher GC\nresponse score")
gax = AX(142.5, 136.0, 40.0, 22.5)
for vv in [0, 1]:
    dd = lv[lv.vent == vv]
    gax.scatter(dd.GC_response, dd.RES_EFF_core, s=3, color=CONTROL if vv == 0 else GTXc, lw=0, alpha=0.75, zorder=1 + vv)
co = np.polyfit(lv.GC_response, lv.RES_EFF_core, 1); xx = np.linspace(lv.GC_response.min(), lv.GC_response.max(), 50)
gax.plot(xx, np.polyval(co, xx), color="#111111", lw=1.0, zorder=4)
r = s7("GTEx_LV", "per_SD_GC_response", "RES_EFF_core", test="G2")
assert f"{r.p:.1e}" == "2.3e-29"
gax.text(0.02, 0.98, f"Adjusted: {r.beta:.2f} SD per SD", transform=gax.transAxes, va="top", fontsize=9)
pb, pe = psplit(r.p)
ptext(gax, 0.02, 0.81, f"({r.ci_lo:.2f} to {r.ci_hi:.2f}); P={pb}", pe, transform=gax.transAxes, va="top", ha="left")
reg("G", f"{r.beta:.2f} SD per SD ({r.ci_lo:.2f} to {r.ci_hi:.2f})", f"{r.beta:.4f} ({r.ci_lo:.4f} to {r.ci_hi:.4f})", SRC["s7b"],
    "G2 GTEx_LV per_SD_GC_response RES_EFF_core (adjusted model)", "yes", lock="0.25 SD higher per SD of the 8-gene GC response score (95% CI: 0.21 to 0.28")
reg("G", f"P={pb}^{pe}", f"{r.p:.3g}", SRC["s7b"], "same row", "yes", lock="P=2.3×10⁻²⁹")
gax.set_ylim(lv.RES_EFF_core.min() - 0.3, lv.RES_EFF_core.max() + 2.1); mticks(gax, [-2, 0, 2], "y")
gax.set_xlabel("8-gene GC response\nscore (SD)", labelpad=1, linespacing=1.0); gax.set_ylabel("Resident\nefferocytic (SD)", labelpad=1, linespacing=1.0)
gax.margins(x=0.04); mticks(gax, [-2, 0, 2])
SHEETS["G_GTEx_LV_scatter"] = lv[["SAMPID", "vent", "GC_response", "RES_EFF_core"]].assign(display_ols_slope=co[0], display_ols_intercept=co[1])

# ================================================================= H aging forest
letter(0, 176.5, "H"); title(6, 176.5, "Aging slopes are small; G5 passed in GTEx LV")
hp = AX(0, 183.0, 110, 9.0); hp.set_axis_off(); hp.set_xlim(0, 110); hp.set_ylim(9.0, 0)
pill(hp, 6, 0.6, 17, 4.3, "G5 passed", "passed")
rG5 = s7("GTEx_LV", "age_per_decade_non_ventilator_DTHHRDY1-4", "RES_EFF_core", test="G5")
l1 = f"LV, other deaths: {f(rG5.beta,3)} SD per decade ({f(rG5.ci_lo,3,False)} to {f(rG5.ci_hi,3,False)});"
l2 = f"equivalence P={rG5.p_tost:.2g}; n={int(rG5.n)}"
hp.text(25, 2.75, l1, va="center", fontsize=9); hp.text(25, 6.8, l2, va="center", fontsize=9)
reg("H", l1, f"{rG5.beta:.4f} ({rG5.ci_lo:.4f} to {rG5.ci_hi:.4f})", SRC["s7b"], "G5 GTEx_LV non-ventilator RES per decade", "yes",
    lock="−0.019 SD per decade; 95% CI: −0.078 to 0.040")
reg("H", l2, f"p_tost={rG5.p_tost:.4g}; n={int(rG5.n)}", SRC["s7b"], "G5 p_tost and n", "yes", lock="equivalence P=0.0037; n=197")
reg("H", "G5 passed", "passed", SRC["s7b"], "pre-specified outcome", "yes", lock="pre-specified test G5 passed")
hrows = [("GTEx_LV", "age_per_decade_non_ventilator_DTHHRDY1-4", "GTEx LV, other", GTXc, "D"),
         ("GTEx_AA", "age_per_decade_non_ventilator_DTHHRDY1-4", "GTEx AA, other", GTXc, "D"),
         ("GTEx_LV", "age_per_decade_ventilator_only", "GTEx LV, ventilator", GTXc, "D"),
         ("GTEx_AA", "age_per_decade_ventilator_only", "GTEx AA, ventilator", GTXc, "D"),
         ("Resident", None, "GSE resident mac.", GSEc, GSEm), ("Monocyte_derived", None, "GSE mono-derived", GSEc, GSEm)]
H_rows = []
for j, (mod, mlab) in enumerate([("RES_EFF_core", "Resident efferocytic"), ("MHCII_core", "MHC-II")]):
    ax = AX(44 + j * 34, 195.5, 30, 19.5)
    ax.axvspan(-0.1, 0.1, color="#E3E3E3", lw=0, zorder=0)
    ax.axvline(0, color="#7F7F7F", lw=0.6, ls="--", zorder=0)
    for i, (ds, ct, lab, col, mk) in enumerate(hrows):
        if ds.startswith("GTEx"):
            r = s7(ds, ct, mod, test="G5"); src = SRC["s7b"]
        else:
            q = s11a[(s11a.subset == ds) & (s11a.module == mod) & (s11a.analysis == "ND_ageing_state") & (s11a.model == "pooled_study_FE")]
            assert len(q) == 1; r = q.iloc[0]; src = SRC["s11a"]
        n = int(r.n); pt = float(r.p_tost); eq = pt < 0.05
        ax.plot([r.ci_lo, r.ci_hi], [i, i], color=col, lw=0.8, zorder=2)
        if eq: ax.plot(r.beta, i, marker="o", ms=7.4, mfc="none", mec="black", mew=0.9, zorder=3, ls="none")
        ax.plot(r.beta, i, marker=mk, ms=3.8 if mk == "D" else 4.4, color=col, ls="none", zorder=4)
        if j == 0:
            ax.text(-0.04, i, f"{lab} ({n})", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=9)
            if ds.startswith("GTEx"):
                lk = {0: "equivalence P=0.0037; n=197", 2: "LV n=450 (253 ventilator-managed)", 3: "atrial appendage n=461 (215 ventilator-managed)"}.get(i)
                reg("H", f"({n})", n, src, f"n {lab}", "yes" if lk else "S", lock=lk)
            else:
                reg("H", f"({n})", n, src, f"n non-failing donors {lab}", "S")
        H_rows.append(dict(module=mod, stratum=lab, n=n, beta=r.beta, ci_lo=r.ci_lo, ci_hi=r.ci_hi, p=r.p, p_tost=pt,
                           equivalence_ring=eq, source=src))
        if mod == "RES_EFF_core" and not ds.startswith("GTEx"):
            reg("H", "", f"{r.beta:.2f} per decade, P={r.p:.2g}", src, f"{ds} ND age per decade, pooled_study_FE", "yes",
                lock="(−0.13 SD per decade; P=0.032)" if ds == "Resident" else "(−0.14 SD per decade; P=0.024)", where="legend")
    ax.set_ylim(5.55, -0.55); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlim(-0.36, 0.26); mticks(ax, [-0.2, 0, 0.2])
    ax.set_title(mlab, fontsize=9, pad=2)
    ax.axhline(3.5, color="#CCCCCC", lw=0.5, zorder=0)
FT(76, 221.2, "Change per decade (SD, 95% CI)", ha="center", fontsize=10)
kk = AX(0, 225.6, 110, 4.3); kk.set_axis_off(); kk.set_xlim(0, 110); kk.set_ylim(0, 1)
kk.add_patch(Rectangle((8, 0.2), 4, 0.6, color="#E3E3E3", lw=0))
kk.text(13, 0.5, "\u00b10.1 SD equivalence bounds", va="center", fontsize=9)
kk.plot(57, 0.5, marker="o", ms=7.4, mfc="none", mec="black", mew=0.9); kk.plot(57, 0.5, marker="D", ms=3.8, color=GTXc)
kk.text(59.5, 0.5, "equivalence P<0.05 (TOST)", va="center", fontsize=9)
SHEETS["H_aging_forest"] = pd.DataFrame(H_rows)

# ================================================================= I per-donor aging
letter(116, 176.5, "I"); title(122, 176.5, "Per-donor scores by age")
iax = AX(127.5, 194.5, 25.5, 20.0)
nv = lv[lv.vent == 0].copy(); assert len(nv) == 197
bins = ["20-29", "30-39", "40-49", "50-59", "60-69", "70-79"]
I_rows = []
for k, b in enumerate(bins):
    v = nv.loc[nv.AGE == b, "RES_EFF_core"].values
    iax.boxplot(v, positions=[k], widths=0.62, showfliers=False, patch_artist=True,
                boxprops=dict(fc="none", ec="#333333", lw=0.7), medianprops=dict(color="#111111", lw=1.0),
                whiskerprops=dict(lw=0.6, color="#333333"), capprops=dict(lw=0.6, color="#333333"))
    iax.scatter(k + rng.uniform(-0.22, 0.22, len(v)), v, s=2.0, color=CONTROL, lw=0, zorder=2)
    iax.text(k, -0.27, str(len(v)), transform=iax.get_xaxis_transform(), ha="center", va="top", fontsize=8)
    reg("I", str(len(v)), len(v), SRC["lv"], f"GTEx LV analysis_ok & vent==0 & AGE=={b}", "new")
    I_rows += [dict(panel="I_left_GTEx_LV_other", age_bin=b, RES_EFF_core=x_) for x_ in v]
iax.set_xticks(range(6)); iax.set_xticklabels(["20", "30", "40", "50", "60", "70"], fontsize=8)
iax.set_xlim(-0.6, 5.6)
iax.text(-0.12, -0.27, "n", transform=iax.transAxes, ha="right", va="top", fontsize=8, style="italic")
iax.set_ylabel("Resident\nefferocytic (SD)", labelpad=1, linespacing=1.0)
FT(140.25, 185.0, "GTEx LV,\nother deaths", ha="center", color=GTXc, linespacing=1.0)
FT(155.0, 224.6, "Age (y)", ha="center", fontsize=10)
jax = AX(158.5, 194.5, 24.0, 20.0)
sr = sub[(sub.subset == "Resident") & (sub.isND == 1)].copy(); assert len(sr) == 138
main3 = ["Chaffin 2022", "Koenig 2022", "Reichart 2022"]
oth = sr[~sr.study.isin(main3)]
jax.scatter(oth.age, oth.RES_EFF_core, s=5, marker="o", fc="none", ec=CONTROL, lw=0.5, zorder=2)
for st in main3:
    dd = sr[sr.study == st]; key = st.split()[0]
    jax.scatter(dd.age, dd.RES_EFF_core, s=7, marker=STUDY_MARKER[key], color=GSE_STUDY[key], lw=0, alpha=0.95, zorder=3)
yl = (min(nv.RES_EFF_core.min(), sr.RES_EFF_core.min()) - 0.3, max(nv.RES_EFF_core.max(), sr.RES_EFF_core.max()) + 0.3)
for ax_ in (iax, jax): ax_.set_ylim(*yl); mticks(ax_, [-2, 0, 2], "y")
jax.set_yticklabels([]); jax.set_xticks([20, 50, 80])
FT(170.5, 185.0, "GSE290367", ha="center", color=GSEc)
I_rows += [dict(panel="I_right_GSE290367_ND_resident", donor_id=a_, study=b_, age=c_, RES_EFF_core=d_)
           for a_, b_, c_, d_ in zip(sr.donor_id, sr.study, sr.age, sr.RES_EFF_core)]
ik = AX(157.5, 188.6, 25.5, 5.6); ik.set_axis_off(); ik.set_xlim(0, 25.5); ik.set_ylim(5.6, 0)
for (st, x0, y0) in [("Chaffin", 0.9, 1.3), ("Koenig", 13.6, 1.3), ("Reichart", 0.9, 4.3)]:
    ik.plot(x0, y0, marker=STUDY_MARKER[st], ms=3.2, color=GSE_STUDY[st], ls="none"); ik.text(x0 + 1.3, y0, st, va="center", fontsize=7.5)
ik.plot(13.6, 4.3, marker="o", ms=3.0, mfc="none", mec=CONTROL, mew=0.6, ls="none")
n_oth = int(oth.study.nunique()); ik.text(14.9, 4.3, f"{n_oth} other", va="center", fontsize=7.5)
reg("I", f"{n_oth} other", n_oth, SRC["sub"], "number of other GSE290367 studies contributing ND resident donors", "new")
SHEETS["I_aging_per_donor"] = pd.DataFrame(I_rows)


# ---- render-then-verify: text overlap check (figure-style 9.1)
import matplotlib as _m
fig.canvas.draw(); _r = fig.canvas.get_renderer()
_t = [(t, t.get_window_extent(_r)) for t in fig.findobj(_m.text.Text) if t.get_text().strip() and t.get_visible()]
_ov = [(a.get_text()[:30], b.get_text()[:30]) for i, (a, ba) in enumerate(_t) for b, bb in _t[i+1:] if ba.overlaps(bb)
       and not (b.get_fontsize() == 7.2 and a.get_text().endswith("\u00d710"))]   # exponent vs its own base is intended
_fb = fig.bbox
_out = [t.get_text()[:30] for t, bb in _t if bb.x0 < _fb.x0 - 1 or bb.x1 > _fb.x1 + 1 or bb.y0 < _fb.y0 - 1 or bb.y1 > _fb.y1 + 1]
print("TEXT OVERLAPS", _ov); print("OUTSIDE", _out)
# ================================================================= save
matplotlib.rcParams["savefig.pad_inches"] = 0.0   # keep the 183-mm canvas; tight bbox only trims
save_figure(fig, "Figure_7")
pd.DataFrame(REG).to_csv("Figure_7_number_registry.csv", index=False)
with pd.ExcelWriter("Figure_7_source_data.xlsx") as xw:
    for k_, v_ in SHEETS.items(): v_.to_excel(xw, sheet_name=k_[:31], index=False)
print("saved Figure_7; registry rows", len(REG), "mismatches", sum(r_["check"].startswith("MISMATCH") for r_ in REG))
