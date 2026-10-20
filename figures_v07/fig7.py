"""Figure 7 (R6 + R7): macrophage programs after GC adjustment, GTEx ventilator proxy, aging.
v0.7: reader-facing labels only (plain test names and outcomes; no internal test codes or plan versions).
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

import os  # public package: inputs resolved relative to this script (override with JBTS_INPUTS)
HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.environ.get("JBTS_INPUTS", os.path.join(HERE, "..", "figures_v05", "inputs"))
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
manuscript = open(MANU, encoding="utf-8").read() if os.path.exists(MANU) else None  # manuscript text is not distributed
MANU6 = os.path.join(INP, "Manuscript_JACC-BTS_v0.6_text.md")
MS6 = open(MANU6, encoding="utf-8").read() if os.path.exists(MANU6) else None
# plain reader-facing names (binding dictionary v0.7)
NAME = {"G3": "GC-insensitive receptor set", "G4": "GC-induced efferocytic genes",
        "RES_EFF_core": "Resident efferocytic", "MHCII_core": "MHC-II", "GC_response": "8-gene GC response score"}
GSNAME = {"G3_GC_insensitive": "GC-insensitive receptor set", "G4_GC_induced": "GC-induced efferocytic genes"}

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

# ---------------------------------------------------------------- registry / formatting (v0.6: every text >= 10 pt)
FS = 10            # every non-title text element
matplotlib.rcParams["font.family"] = [matplotlib.rcParams["font.sans-serif"][0], "FreeSans"]   # glyph fallback for ⁻ superscripts
FT_TICK = SIZES["tick"]; FT_AX = SIZES["axis_label"]; FT_TTL = SIZES["panel_title"]
REG, SHEETS = [], {}
MS = manuscript
def _n(s): return s.replace("\u2212", "-")
def reg(panel, text, value, src, comp, status, lock=None, where="figure"):
    note = ""
    if status == "yes":
        assert lock is not None, text
        if MS is None:
            note = "not checked (manuscript text not distributed)"
        else:
            ok = _n(lock) in _n(MS)
            note = "verbatim in v0.4 text" if ok else "MISMATCH: phrase not found in v0.4 text"
            if not ok: print("LOCK MISMATCH", panel, lock)
    REG.append(dict(figure="Figure 7", panel=panel, location=where, printed_text=text.replace("\n", " / "), value=value,
                    source_file=src, computation=comp, matches_locked_text=status, locked_text_phrase=lock or "", check=note,
                    phrase_in_v06_text=("" if (lock is None or MS6 is None) else ("yes" if _n(lock) in _n(MS6) else "no"))))
def f(x, nd=2, sign=True):
    s = f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"
    return s.replace("-", "\u2212")
SUPC = str.maketrans("-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079")
def fp(p):
    """P in the text style; below 1e-4 inline Unicode superscripts at full size (e.g. 5.3×10⁻⁵)."""
    if p < 1e-4:
        m, e = f"{p:.1e}".split("e"); return f"{m}\u00d710{str(int(e)).translate(SUPC)}"
    return f"{p:.2g}"

# ---------------------------------------------------------------- canvas helpers (mm, top-left origin)
W, H = 183.0, 240.0
fig = plt.figure(figsize=(W * MM, H * MM))
def AX(x, y, w, h, **kw):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **kw)
def FT(x, y, s, **kw):
    kw.setdefault("fontsize", FS); kw.setdefault("va", "top"); kw.setdefault("ha", "left")
    return fig.text(x / W, 1 - y / H, s, **kw)
def letter(x, y, L): FT(x, y, L, fontsize=SIZES["panel_letter"], fontweight="bold")
def title(x, y, t): FT(x, y + 0.7, t, fontsize=FT_TTL, fontweight="bold", linespacing=1.05)
def banner(y, t):
    a = AX(0, y, W, 4.8); a.set_axis_off()
    a.add_patch(Rectangle((0, 0), 1, 1, color="#666666", lw=0, transform=a.transAxes))
    a.text(0.005, 0.5, t, color="white", fontsize=FS, va="center", ha="left", fontweight="bold", transform=a.transAxes)
GSEc, MAGc, GTXc = COHORT["GSE290367"], COHORT["MAGNet"], COHORT["GTEx"]
GSEm, MAGm = COHORT_MARKER["GSE290367"], COHORT_MARKER["MAGNet"]
MINUS = "\u2212"
def mticks(ax, ticks, axis="x"):
    lab = [(f"{t:g}").replace("-", MINUS) for t in ticks]
    if axis == "x": ax.set_xticks(ticks); ax.set_xticklabels(lab)
    else: ax.set_yticks(ticks); ax.set_yticklabels(lab)
rng = np.random.default_rng(7)   # jitter only (display)
def box(ax, x, y, w, h, txt, fc="white", ec="#555555", lw=0.7):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2", fc=fc, ec=ec, lw=lw))
    ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=FS, color=TEXT, linespacing=1.1)
def arrow(ax, x0, y0, x1, y1, c="#555555"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7, lw=0.8, color=c))
def pill(ax, x, y, w, h, txt, kind):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={h/2}", fc="black" if kind == "passed" else "white",
                       ec="black", lw=0.8, ls="--" if kind == "noteval" else "-")
    ax.add_patch(p)
    if kind == "notmet":
        hb = Rectangle((x, y), w * 0.16, h, fc="black", ec="none"); ax.add_patch(hb); hb.set_clip_path(p)
        ax.text(x + w * 0.58, y + h / 2, txt, ha="center", va="center", fontsize=FS)
    elif kind == "failed":
        ax.text(x + 2.6, y + h / 2, "\u00d7", ha="center", va="center", fontsize=FS + 1, fontweight="bold")
        ax.text(x + w / 2 + 1.6, y + h / 2, txt, ha="center", va="center", fontsize=FS)
    else:
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=FS, color="white" if kind == "passed" else TEXT)

def outcome(ax, x, y, w, h, txt, kind="notmet"):
    """Plain outcome box: half-black tab = planned test not met; black = met; dashed = not evaluated."""
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.6", fc="black" if kind == "met" else "white",
                       ec="black", lw=0.8, ls="--" if kind == "noteval" else "-")
    ax.add_patch(p)
    if kind == "notmet":
        hb = Rectangle((x, y), 4.0, h, fc="black", ec="none"); ax.add_patch(hb); hb.set_clip_path(p)
        ax.text(x + 4.0 + (w - 4.0) / 2, y + h / 2, txt, ha="center", va="center", fontsize=FS, linespacing=1.1)
    else:
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=FS, linespacing=1.1,
                color="white" if kind == "met" else TEXT)
# ================================================================= banners
banner(0, "Macrophage modules after GC adjustment (GSE290367 resident macrophages; MAGNet)")
banner(124.0, "GTEx: ventilator-managed vs other deaths (supportive comparison, not blinded)")
banner(181.0, "Aging (planned equivalence test)")

# ================================================================= A schematic (full-width logic strip)
letter(0, 6.0, "A"); title(6, 6.0, "What each module tests")
a = AX(0, 13.0, W, 24.0); a.set_xlim(0, W); a.set_ylim(24.0, 0); a.set_axis_off()
n_g3 = int((s3.gene_set == "G3_GC_insensitive").sum()); n_g4 = int((s3.gene_set == "G4_GC_induced").sum())
assert (n_g3, n_g4) == (4, 14)
g3_genes = ", ".join(s3[s3.gene_set == "G3_GC_insensitive"].gene)
rowsA = [(f"{NAME['G3']}\n({g3_genes})", "Prediction: decline persists\nafter adjustment", GC_CLASS["insensitive"],
          "Planned test not met:\nno decline in either cohort"),
         (f"{NAME['G4']} ({n_g4})\n(C1QA\u2013C, CD163, MERTK, \u2026)", "Prediction: decline\nattenuates by \u226550%", GC_CLASS["induced"],
          "Planned test not met:\nattenuation 37%, 34%")]
for r_, (t1, t2, ec, t3) in enumerate(rowsA):
    y0 = 0.5 + r_ * 12.0
    box(a, 0.5, y0, 57.0, 10.0, t1, ec=ec, lw=1.5)
    arrow(a, 58.0, y0 + 5, 61.5, y0 + 5)
    box(a, 62.0, y0, 57.0, 10.0, t2, ec="#BBBBBB")
    arrow(a, 120.0, y0 + 5, 123.5, y0 + 5)
    outcome(a, 124.0, y0, 58.5, 10.0, t3, "notmet")
reg("A", "Planned test not met: no decline in either cohort", "not met", "Table 1 (v0.4)", "planned outcome, GC-independent loss test", "yes",
    lock="because the GC-insensitive receptor set did not decline in DCM")
reg("A", "Planned test not met: attenuation 37%, 34%", "not met", "Table 1 (v0.4)", "planned outcome, attenuation test", "yes",
    lock="This was below the pre-specified 50% threshold, so prediction")
reg("A", "37%, 34%", "37; 34", SRC["s7b"], "attenuation of GC-induced efferocytic genes (base - +GC)/base, GSE290367; MAGNet (models adjusted for the 8-gene GC response score)", "yes", lock="attenuation of 37% and 34%")
reg("A", f"{NAME['G3']} ({g3_genes})", n_g3, SRC["s3"], "count of genes in GC-insensitive receptor set", "yes", lock="GC-insensitive receptor set (ITGAV, ITGB5, LRP1, SCARB1)")
reg("A", f"{NAME['G4']} ({n_g4})", n_g4, SRC["s3"], "count of GC-induced efferocytic genes", "yes", lock="GC-induced efferocytic subset (14 genes)")
reg("A", "\u226550%", 50, "Table 1 (v0.4)", "planned threshold", "yes", lock="pre-specified 50% threshold")
SHEETS["A_gene_sets"] = (s3[s3.gene_set.isin(["G3_GC_insensitive", "G4_GC_induced"])][["gene_set", "gene", "meta_log2FC", "GC_class"]]
                         .assign(gene_set=lambda d: d.gene_set.map(GSNAME)))

# ================================================================= B strip + box (1 x 4, shared y)
letter(0, 39.5, "B"); title(6, 39.5, "Gene-set scores per donor")
dsets = [("GSE290367", gse, GSEc, GSEm), ("MAGNet", mag, MAGc, MAGm)]
lo = min(gse[["G3", "G4"]].min().min(), mag[["G3", "G4"]].min().min()); hi = max(gse[["G3", "G4"]].max().max(), mag[["G3", "G4"]].max().max())
pad = 0.05 * (hi - lo)
B_rows = []; bx = [13.5, 37.0, 63.5, 87.0]; bw = 22.0
for i, mod in enumerate(["G3", "G4"]):
    FT(bx[2 * i] + bw + 1.75, 45.6, {"G3": "GC-insensitive receptors", "G4": "GC-induced efferocytic"}[mod], ha="center", fontweight="bold")
    for j, (nm, d, col, mk) in enumerate(dsets):
        ax = AX(bx[2 * i + j], 55.0, bw, 17.0)
        for k, cs in enumerate([0, 1]):
            v = d.loc[d.case == cs, mod].dropna().values
            ax.boxplot(v, positions=[k], widths=0.6, showfliers=False, patch_artist=True, zorder=1,
                       boxprops=dict(fc="none", ec="#333333", lw=0.7), medianprops=dict(color="#111111", lw=1.1),
                       whiskerprops=dict(lw=0.6, color="#333333"), capprops=dict(lw=0.6, color="#333333"))
            ax.scatter(k + rng.uniform(-0.22, 0.22, len(v)), v, s=4, marker=mk, color=CONTROL if cs == 0 else col, lw=0, alpha=0.85, zorder=2)
            B_rows += [dict(gene_set=NAME[mod], cohort=nm, group="DCM" if cs else "control", score=x) for x in v]
        ax.set_xlim(-0.6, 1.6); ax.set_xticks([0, 1]); ax.set_xticklabels(["Ctrl", "DCM"])
        ax.axhline(0, color="#7F7F7F", lw=0.6, ls="--", zorder=0)
        ax.set_ylim(lo - pad, hi + pad); mticks(ax, [-4, 0, 4], "y")
        if i == 0 and j == 0: ax.set_ylabel("Score (SD)", labelpad=1, fontsize=FT_AX)
        else: ax.set_yticklabels([])
        FT(bx[2 * i + j] + bw / 2, 50.0, nm if nm == "MAGNet" else "GSE290367", ha="center", color=col)
SHEETS["B_gene_sets_per_donor"] = pd.DataFrame(B_rows)
reg("B", "", "45/64; 102/128", f"{SRC['gse']}; {SRC['mag']}", "row counts by case (control/case)", "yes", lock="64 DCM and 45 controls", where="legend")

# ================================================================= D scatter GC vs G4
letter(0, 80.0, "D"); title(6, 80.0, "GC-induced efferocytic genes vs GC response score")
FT(6, 85.6, "Descriptive; planned correlation test not evaluated. Gray, control.", style="italic", color="#7F7F7F")
D_rows = []
for j, (nm, d, col, mk) in enumerate(dsets):
    ax = AX(13.5 + j * 51, 94.5, 45.0, 18.3)
    for cs in [0, 1]:
        dd = d[d.case == cs]
        ax.scatter(dd.GC_response, dd.G4, s=7, marker=mk, color=CONTROL if cs == 0 else col, lw=0, alpha=0.85, zorder=2 if cs else 1)
        D_rows += [dict(cohort=nm, group="DCM" if cs else "control", GC_response_8gene=a_, GC_induced_efferocytic=b_) for a_, b_ in zip(dd.GC_response, dd.G4)]
    if j == 0: ax.set_ylabel("GC-induced\ngenes (SD)", labelpad=1, linespacing=1.0, fontsize=FT_AX)
    ax.margins(0.05)
    mticks(ax, [-4, 0, 4] if j == 0 else [-4, -2, 0], "x"); mticks(ax, [-4, 0], "y")
    FT(13.5 + j * 51, 90.6, "GSE290367 resident mac." if j == 0 else "MAGNet bulk LV", color=col)
FT(61.0, 119.4, "8-gene GC response score (SD)", ha="center", fontsize=FT_AX)
SHEETS["D_GC_vs_GC_induced"] = pd.DataFrame(D_rows)
reg("D", "planned correlation test not evaluated", "not evaluated", "Table 1 (v0.4)", "planned outcome, correlation test", "yes", lock="was not evaluated (Table 1)")

# ================================================================= C forest (module x cohort x adjustment)
letter(117, 39.5, "C"); title(123, 39.5, "Module effects without and\nwith GC adjustment")
cax = AX(118.5, 54.0, 25.5, 45.0)
mods = [("G3", f"GC-insensitive receptors ({n_g3})"), ("G4", f"GC-induced efferocytic ({n_g4})"), ("RES_EFF_core", "Resident efferocytic (9)"), ("MHCII_core", "MHC-II (10)")]
cohorts = [("GSE290367", "GSE290367_Resident", "DCM_vs_ND", "age+sex+study", "age+sex+study+GC", "GSE290367_Resident_DCM", GSEc, GSEm),
           ("MAGNet", "MAGNet_excl_overlap", "DCM_vs_NF", "age+sex", "age+sex+GC", "MAGNet_exclOverlap_DCM", MAGc, MAGm)]
locks = {("G3", "GSE290367"): "unadjusted estimates were −0.04 and +0.31 SD", ("G3", "MAGNet"): "+0.43 to +0.51 SD in MAGNet",
         ("G4", "GSE290367"): "from −1.18 to −0.75 SD after adjustment in GSE290367", ("G4", "MAGNet"): "from −0.97 to −0.65 SD in MAGNet",
         ("RES_EFF_core", "GSE290367"): "(−2.58 to −2.27 SD", ("RES_EFF_core", "MAGNet"): "MAGNet (−1.53 to −1.06 SD)",
         ("MHCII_core", "GSE290367"): "GSE290367 +0.94 to +0.80 SD", ("MHCII_core", "MAGNet"): "MAGNet +1.50 to +1.56 SD"}
C_rows = []; y = 0.0
TXa, TXb = 145.5, 173.0     # text columns (mm)
def ctext(xmm, yy, s, col):
    # place text in figure mm at the data-y of the forest
    cax.text((xmm - 118.5) / 25.5, yy, s, transform=cax.get_yaxis_transform(), ha="left", va="center", fontsize=FS, color=col, clip_on=False)
for mod, mlab in mods:
    cax.text(0.0, y, mlab, transform=cax.get_yaxis_transform(), ha="left", va="center", fontsize=FS, zorder=5,
             bbox=dict(fc="white", ec="none", pad=0.3))
    for ci, (cn, ds, ct, adj0, adj1, r2ds, col, mk) in enumerate(cohorts):
        yy = y + 1.0 + ci
        b0 = s7(ds, ct, mod, adj0); b1 = s7(ds, ct, mod, adj1)
        dys = [(b0, True, -0.25), (b1, False, 0.0)]
        for row, filled, dy in dys:
            cax.plot([row.ci_lo, row.ci_hi], [yy + dy] * 2, color=col, lw=0.8, zorder=1)
            cax.plot(row.beta, yy + dy, marker=mk, ms=4.2, mfc=col if filled else "white", mec=col, mew=0.9, zorder=3, ls="none")
            C_rows.append(dict(module=NAME[mod], cohort=cn, adjustment=row.adjustment, marker="filled" if filled else "open",
                               model="planned; base or adjusted for the 8-gene GC response score",
                               beta=row.beta, ci_lo=row.ci_lo, ci_hi=row.ci_hi, p=row.p, source="Supplemental Tables, sheet S7b"))
        dag = "\u2020" if (mod == "RES_EFF_core" and cn == "GSE290367") else ""
        txt = f"{f(b0.beta)} \u2192 {f(b1.beta)}{dag}"
        ctext(TXa, yy, txt, col)
        reg("C", txt, f"{b0.beta:.4f} -> {b1.beta:.4f}", SRC["s7b"], f"{ds} {ct} {NAME[mod]}: {adj0} -> {adj1} (DCM - control)", "yes", lock=locks[(mod, cn)])
        if mod == "G3" and cn == "GSE290367":
            reg("C", "", f"{b1.ci_lo:.2f} to {b1.ci_hi:.2f}", SRC["s7b"], "GC-insensitive receptor set, +8-gene GC score, 95% CI, GSE290367", "yes", lock="−0.37 SD (95% CI: −0.93 to 0.19)", where="legend")
        if mod == "G4":
            h_ = r2row(r2ds, "G4", "+GC3+MAC")
            cax.plot([h_.ci_lo, h_.ci_hi], [yy + 0.25] * 2, color=col, lw=0.8, zorder=1)
            cax.plot(h_.beta, yy + 0.25, marker=mk, ms=4.2, fillstyle="left", mfc=col, mfcalt="white", mec=col, mew=0.9, zorder=3, ls="none")
            C_rows.append(dict(module=NAME[mod], cohort=cn, adjustment="+3-gene GC score + macrophage content (post hoc)", marker="half",
                               model="post hoc; adjusted for the 3-gene GC score plus macrophage content", beta=h_.beta, ci_lo=h_.ci_lo, ci_hi=h_.ci_hi, p=h_.p, source="post hoc model results (rev_posthoc_models.csv, archived with the code)"))
            ctext(TXb, yy, f(h_.beta), col)
            reg("C", f(h_.beta), f"{h_.beta:.4f}", SRC["r2"], f"post hoc model {r2ds} GC-induced efferocytic genes, +3-gene GC score + macrophage content", "yes",
                lock="+0.04 SD (95% CI: −0.51 to 0.58) in GSE290367" if cn == "GSE290367" else "−0.42 SD (P<10⁻⁴) in MAGNet")
    y += 3.25
cax.axvline(0, color="#7F7F7F", lw=0.6, ls="--", zorder=0)
cax.set_ylim(y - 0.6, -0.75); cax.set_yticks([]); cax.spines["left"].set_visible(False)
cax.set_xlim(-3.4, 2.1); mticks(cax, [-2, 0, 2])
cax.set_xlabel(f"DCM {MINUS} control (SD)", labelpad=1, fontsize=FT_AX)
FT(TXa, 50.6, "base \u2192 +GC"); FT(TXb, 50.6, "+mac.")
SHEETS["C_module_effects"] = pd.DataFrame(C_rows)
# attenuation of GC-induced genes (3-gene, post hoc) values go to the legend; 8-gene values are printed in A
for cn, ds, ct, adj0, adj1, r2ds, col, mk in cohorts:
    r0 = r2row(r2ds, "G4", "base").beta; r3 = r2row(r2ds, "G4", "+GC3").beta; a3 = 100 * (r0 - r3) / r0
    reg("C", "", f"{a3:.2f}", SRC["r2"], f"attenuation of GC-induced efferocytic genes, 3-gene GC score, {cn} (post hoc)", "yes", lock="similar attenuation (31% and 30%)", where="legend")
kax = AX(117.0, 110.0, 66.0, 13.6); kax.set_axis_off(); kax.set_xlim(0, 66); kax.set_ylim(13.6, 0)
def kmk(x, yk, mk, col, fs="full", ms=4.4):
    kax.plot(x, yk, marker=mk, ms=ms, ls="none", mew=0.9, mec=col,
             mfc=("white" if fs == "none" else col), fillstyle=("left" if fs == "left" else "full"), mfcalt="white")
kmk(1.5, 2.2, GSEm, GSEc); kax.text(3.5, 2.2, "GSE290367", va="center", fontsize=FS, color=GSEc)
kmk(28.0, 2.2, MAGm, MAGc, ms=4.0); kax.text(30.0, 2.2, "MAGNet", va="center", fontsize=FS, color=MAGc)
kmk(1.5, 6.8, "o", "#555555"); kax.text(3.5, 6.8, "base", va="center", fontsize=FS)
kmk(15.0, 6.8, "o", "#555555", "none"); kax.text(17.0, 6.8, "+ 8-gene GC score", va="center", fontsize=FS)
kmk(1.5, 11.3, "o", "#555555", "left"); kax.text(3.5, 11.3, "+ 3-gene GC + mac. (post hoc)", va="center", fontsize=FS)

# ================================================================= E forest GTEx
letter(0, 130.0, "E"); title(6, 130.0, "Ventilator-managed: efferocytic up,\nMHC-II down, GC score not higher")
eax = AX(24.0, 145.5, 22.0, 24.5)
emods = [("GC_response", "GC response\n(8-gene)"), ("RES_EFF_core", "Resident\nefferocytic"), ("MHCII_core", "MHC-II")]
elocks = {("LV", "GC_response"): "(+0.16 SD; P=0.34)", ("LV", "RES_EFF_core"): "LV (+0.30 SD; P=0.0021)",
          ("AA", "RES_EFF_core"): "atrial appendage (+0.14 SD; P=0.0077)", ("LV", "MHCII_core"): "(LV −0.36 SD; P=0.0056)"}
E_rows = []
for i, (mod, lab) in enumerate(emods):
    eax.text(-0.06, i * 1.0, lab, transform=eax.get_yaxis_transform(), ha="right", va="center", fontsize=FS, linespacing=1.0)
    for k, (ds, mk, tl) in enumerate([("GTEx_LV", "o", "LV"), ("GTEx_AA", "D", "AA")]):
        r = s7(ds, "ventilator_vs_other", mod, test="G2"); yk = i * 1.0 - 0.21 + 0.42 * k
        eax.plot([r.ci_lo, r.ci_hi], [yk, yk], color=GTXc, lw=0.8)
        eax.plot(r.beta, yk, marker=mk, ms=4.6 if mk == "o" else 3.9, color=GTXc, ls="none")
        txt = f"{f(r.beta)}, P={fp(r.p)}"
        eax.text(1.05, yk, txt, transform=eax.get_yaxis_transform(), ha="left", va="center", fontsize=FS)
        reg("E", txt, f"{r.beta:.4f}; P={r.p:.3g}", SRC["s7b"], f"GTEx ventilator comparison {ds} ventilator_vs_other {NAME[mod]}", "yes" if (tl, mod) in elocks else "S", lock=elocks.get((tl, mod)))
        E_rows.append(dict(tissue=tl, module=NAME[mod], n=r.n, n_vent=r.n_case, beta=r.beta, ci_lo=r.ci_lo, ci_hi=r.ci_hi, p=r.p, adjustment=r.adjustment))
eax.axvline(0, color="#7F7F7F", lw=0.6, ls="--", zorder=0)
eax.set_ylim(2.5, -0.5); eax.set_yticks([]); eax.spines["left"].set_visible(False)
eax.set_xlim(-0.9, 0.6); mticks(eax, [-0.5, 0, 0.5])
eax.set_xlabel(f"Ventilator {MINUS} other (SD)", labelpad=1, fontsize=FT_AX)
ek = AX(0, 139.5, 80, 5.0); ek.set_axis_off(); ek.set_xlim(0, 80); ek.set_ylim(0, 1)
ek.plot(7.0, 0.5, marker="o", ms=4.6, color=GTXc, ls="none"); ek.text(8.8, 0.5, "LV (n=197/253)", va="center", fontsize=FS)
ek.plot(41.0, 0.5, marker="D", ms=3.9, color=GTXc, ls="none"); ek.text(42.8, 0.5, "AA (n=246/215)", va="center", fontsize=FS)
reg("E", "LV (n=197/253)", "450 = 197 other / 253 ventilator", SRC["s7b"], "n GTEx_LV ventilator comparison (other/ventilator)", "yes", lock="LV n=450 (253 ventilator-managed)")
reg("E", "AA (n=246/215)", "461 = 246 other / 215 ventilator", SRC["s7b"], "n GTEx_AA ventilator comparison (other/ventilator)", "yes", lock="atrial appendage n=461 (215 ventilator-managed)")
SHEETS["E_GTEx_ventilator"] = pd.DataFrame(E_rows)

# ================================================================= F violins (raw, unadjusted)
letter(82, 130.0, "F"); title(88, 130.0, "Per-sample LV\nscores (unadjusted)")
fax = AX(93.0, 148.5, 34.0, 20.5)
F_rows = []
for j, col_ in enumerate(["GC_response", "RES_EFF_core"]):
    for k, vv in enumerate([0, 1]):
        pos = j * 2.5 + k
        v = lv.loc[lv.vent == vv, col_].values; c = CONTROL if vv == 0 else GTXc
        vp = fax.violinplot(v, positions=[pos], widths=0.85, showextrema=False)
        for b_ in vp["bodies"]: b_.set_facecolor(c); b_.set_edgecolor("none"); b_.set_alpha(0.35)
        fax.scatter(pos + rng.uniform(-0.13, 0.13, len(v)), v, s=1.5, color=c, lw=0, alpha=0.8, zorder=2)
        fax.plot([pos - 0.28, pos + 0.28], [np.median(v)] * 2, color="#111111", lw=1.2, zorder=3)
        F_rows += [dict(score=col_, group="ventilator" if vv else "other", value=x_) for x_ in v]
fax.set_xticks([0.5, 3.0]); fax.set_xticklabels(["GC\nresponse", "Resident\nefferocytic"], linespacing=1.0)
fax.set_xlim(-0.6, 4.1); fax.set_ylim(-3.3, 3.6); mticks(fax, [-2, 0, 2], "y")
fax.axhline(0, color="#7F7F7F", lw=0.6, ls="--", zorder=0)
fax.set_ylabel("Score (SD)", labelpad=1, fontsize=FT_AX)
fk = AX(84, 139.0, 44, 9.0); fk.set_axis_off(); fk.set_xlim(0, 44); fk.set_ylim(9, 0)
fk.plot(6.0, 2.2, "o", ms=4.2, color=CONTROL); fk.text(7.8, 2.2, "other deaths (197)", va="center", fontsize=FS)
fk.plot(6.0, 6.8, "o", ms=4.2, color=GTXc); fk.text(7.8, 6.8, "ventilator (253)", va="center", fontsize=FS)
reg("F", "other deaths (197); ventilator (253)", "197; 253", SRC["lv"], "counts analysis_ok by vent", "yes", lock="LV n=450 (253 ventilator-managed)")
SHEETS["F_GTEx_LV_per_sample"] = pd.DataFrame(F_rows)

# ================================================================= G scatter GC vs RES (GTEx LV)
letter(131, 130.0, "G"); title(137, 130.0, "Efferocytic score higher\nwith higher GC score")
gax = AX(145.0, 141.5, 37.5, 25.5)
for vv in [0, 1]:
    dd = lv[lv.vent == vv]
    gax.scatter(dd.GC_response, dd.RES_EFF_core, s=3, color=CONTROL if vv == 0 else GTXc, lw=0, alpha=0.75, zorder=1 + vv)
co = np.polyfit(lv.GC_response, lv.RES_EFF_core, 1); xx = np.linspace(lv.GC_response.min(), lv.GC_response.max(), 50)
gax.plot(xx, np.polyval(co, xx), color="#111111", lw=1.0, zorder=4)
r = s7("GTEx_LV", "per_SD_GC_response", "RES_EFF_core", test="G2")
assert f"{r.p:.1e}" == "2.3e-29"
reg("G", f"0.25 SD per SD (0.21 to 0.28); P={fp(r.p)}", f"{r.beta:.4f} ({r.ci_lo:.4f} to {r.ci_hi:.4f}); p={r.p:.3g}", SRC["s7b"],
    "GTEx_LV resident efferocytic score per SD of the 8-gene GC response score (adjusted model); printed in legend", "yes",
    lock="0.25 SD higher per SD of the 8-gene GC response score (95% CI: 0.21 to 0.28; P=2.3×10⁻²⁹)", where="legend")
mticks(gax, [-2, 0, 2], "y"); mticks(gax, [-2, 0, 2])
gax.set_xlabel("GC response (SD)", labelpad=1, fontsize=FT_AX); gax.set_ylabel("Efferocytic\nscore (SD)", labelpad=1, linespacing=1.0, fontsize=FT_AX)
gax.margins(x=0.04)
SHEETS["G_GTEx_LV_scatter"] = lv[["SAMPID", "vent", "GC_response", "RES_EFF_core"]].assign(display_ols_slope=co[0], display_ols_intercept=co[1])

# ================================================================= H aging forest
letter(0, 187.0, "H"); title(6, 187.0, "Aging slopes are small; equivalence met in GTEx LV")
hp = AX(0, 193.5, 112, 9.0); hp.set_axis_off(); hp.set_xlim(0, 112); hp.set_ylim(9.0, 0)
outcome(hp, 0.6, 0.1, 31.0, 8.2, "Equivalence met\n(\u00b10.1 SD/decade)", "met")
rG5 = s7("GTEx_LV", "age_per_decade_non_ventilator_DTHHRDY1-4", "RES_EFF_core", test="G5")
l1 = f"LV, other: {f(rG5.beta,3)} SD/decade ({f(rG5.ci_lo,3,False)} to {f(rG5.ci_hi,3,False)})"
l2 = f"equivalence P={rG5.p_tost:.2g}, n={int(rG5.n)}"
hp.text(33.6, 2.3, l1, va="center", fontsize=FS); hp.text(33.6, 6.9, l2, va="center", fontsize=FS)
reg("H", l1, f"{rG5.beta:.4f} ({rG5.ci_lo:.4f} to {rG5.ci_hi:.4f})", SRC["s7b"], "aging equivalence test: GTEx_LV other deaths, resident efferocytic per decade", "yes",
    lock="−0.019 SD per decade; 95% CI: −0.078 to 0.040")
reg("H", l2, f"p_tost={rG5.p_tost:.4g}; n={int(rG5.n)}", SRC["s7b"], "aging equivalence test: TOST P and n", "yes", lock="equivalence P=0.0037; n=197")
reg("H", "Equivalence met", "met", SRC["s7b"], "planned outcome, aging equivalence test (p_tost<0.05)", "yes", lock="equivalent to no change with age in the LV")
reg("H", "\u00b10.1 SD/decade", 0.1, "Table 1 (v0.4)", "planned equivalence bound", "yes", lock="Equivalence within \u00b10.1 SD per decade")
hrows = [("GTEx_LV", "age_per_decade_non_ventilator_DTHHRDY1-4", "LV, other", GTXc, "D"),
         ("GTEx_AA", "age_per_decade_non_ventilator_DTHHRDY1-4", "AA, other", GTXc, "D"),
         ("GTEx_LV", "age_per_decade_ventilator_only", "LV, ventilator", GTXc, "D"),
         ("GTEx_AA", "age_per_decade_ventilator_only", "AA, ventilator", GTXc, "D"),
         ("Resident", None, "Resident mac.", GSEc, GSEm), ("Monocyte_derived", None, "Mono-derived", GSEc, GSEm)]
H_rows = []
for j, (mod, mlab) in enumerate([("RES_EFF_core", "Resident efferocytic"), ("MHCII_core", "MHC-II")]):
    ax = AX(45.0 + j * 34.0, 207.0, 30.0, 21.5)
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
            ax.text(-0.04, i, f"{lab} ({n})", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=FS)
            if ds.startswith("GTEx"):
                lk = {0: "equivalence P=0.0037; n=197", 2: "LV n=450 (253 ventilator-managed)", 3: "atrial appendage n=461 (215 ventilator-managed)"}.get(i)
                reg("H", f"({n})", n, src, f"n {lab}", "yes" if lk else "S", lock=lk)
            else:
                reg("H", f"({n})", n, src, f"n non-failing donors {lab}", "S")
        H_rows.append(dict(module=NAME[mod], stratum=lab, n=n, beta=r.beta, ci_lo=r.ci_lo, ci_hi=r.ci_hi, p=r.p, p_tost=pt, equivalence_ring=eq,
                           source="Supplemental Tables, sheet S7b" if src == SRC["s7b"] else "Supplemental Tables, sheet S11a"))
        if mod == "RES_EFF_core" and not ds.startswith("GTEx"):
            reg("H", "", f"{r.beta:.2f} per decade, P={r.p:.2g}", src, f"{ds} ND age per decade, pooled_study_FE", "yes",
                lock="(−0.13 SD per decade; P=0.032)" if ds == "Resident" else "(−0.14 SD per decade; P=0.024)", where="legend")
    ax.set_ylim(5.55, -0.55); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlim(-0.36, 0.26); mticks(ax, [-0.2, 0, 0.2])
    FT(45.0 + j * 34.0 + 15.0, 202.6, mlab, ha="center")
    ax.axhline(3.5, color="#CCCCCC", lw=0.5, zorder=0)
# rotated group labels
FT(1.0, 207.0 + 21.5 * (0.05 + 4.0 / 6.1 / 2) - 6.5, "GTEx", rotation=90, va="center", ha="left")
FT(1.0, 207.0 + 21.5 * (4.6 / 6.1) + 0.5, "GSE", rotation=90, va="center", ha="left")
FT(79.5, 235.4, "Change per decade (SD)", ha="center", fontsize=FT_AX)
kk = AX(0, 230.8, 44, 9.0); kk.set_axis_off(); kk.set_xlim(0, 44); kk.set_ylim(9, 0)
kk.add_patch(Rectangle((1.0, 1.0), 3.6, 2.6, color="#E3E3E3", lw=0)); kk.text(5.6, 2.3, "\u00b10.1 SD bounds", va="center", fontsize=FS)
kk.plot(2.8, 6.8, marker="o", ms=7.4, mfc="none", mec="black", mew=0.9); kk.plot(2.8, 6.8, marker="D", ms=3.8, color=GTXc)
kk.text(5.6, 6.8, "equivalent (TOST)", va="center", fontsize=FS)
SHEETS["H_aging_forest"] = pd.DataFrame(H_rows)

# ================================================================= I per-donor aging
letter(115, 187.0, "I"); title(121, 187.0, "Per-donor scores by age")
iax = AX(125.5, 199.0, 27.0, 27.5)
nv = lv[lv.vent == 0].copy(); assert len(nv) == 197
bins = ["20-29", "30-39", "40-49", "50-59", "60-69", "70-79"]
I_rows = []
for k, b_ in enumerate(bins):
    v = nv.loc[nv.AGE == b_, "RES_EFF_core"].values
    iax.boxplot(v, positions=[k], widths=0.62, showfliers=False, patch_artist=True,
                boxprops=dict(fc="none", ec="#333333", lw=0.7), medianprops=dict(color="#111111", lw=1.0),
                whiskerprops=dict(lw=0.6, color="#333333"), capprops=dict(lw=0.6, color="#333333"))
    iax.scatter(k + rng.uniform(-0.22, 0.22, len(v)), v, s=2.0, color=CONTROL, lw=0, zorder=2)
    reg("I", str(len(v)), len(v), SRC["lv"], f"GTEx LV analysis_ok & vent==0 & AGE=={b_} (n per decade, printed in legend)", "new", where="legend")
    I_rows += [dict(panel="I_left_GTEx_LV_other", age_bin=b_, RES_EFF_core=x_) for x_ in v]
iax.set_xticks([0, 2, 4]); iax.set_xticklabels(["20s", "40s", "60s"]); iax.set_xticks([1, 3, 5], minor=True)
iax.set_xlim(-0.6, 5.6)
iax.set_ylabel("Efferocytic\nscore (SD)", labelpad=1, linespacing=1.0, fontsize=FT_AX)
FT(139.0, 193.6, "GTEx LV, other", ha="center", color=GTXc)
jax = AX(156.0, 199.0, 26.5, 27.5)
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
FT(169.25, 193.6, "GSE290367", ha="center", color=GSEc)
FT(154.0, 234.8, "Age (y)", ha="center", fontsize=FT_AX)
n_oth = int(oth.study.nunique()); assert n_oth == 7
reg("I", "", n_oth, SRC["sub"], "number of other GSE290367 studies (legend)", "new", where="legend")
I_rows += [dict(panel="I_right_GSE290367_ND_resident", donor_id=a_, study=b_, age=c_, RES_EFF_core=d_)
           for a_, b_, c_, d_ in zip(sr.donor_id, sr.study, sr.age, sr.RES_EFF_core)]
SHEETS["I_aging_per_donor"] = pd.DataFrame(I_rows)

# ---- render-then-verify
import matplotlib as _m
fig.canvas.draw(); _r = fig.canvas.get_renderer()
_t = [(t, t.get_window_extent(_r)) for t in fig.findobj(_m.text.Text) if t.get_text().strip() and t.get_visible()]
_ov = [(a_.get_text()[:30], b_.get_text()[:30]) for i, (a_, ba) in enumerate(_t) for b_, bb in _t[i+1:] if ba.overlaps(bb)]
_fb = fig.bbox
_out = [t.get_text()[:30] for t, bb in _t if bb.x0 < _fb.x0 - 1 or bb.x1 > _fb.x1 + 1 or bb.y0 < _fb.y0 - 1 or bb.y1 > _fb.y1 + 1]
print("TEXT OVERLAPS", _ov); print("OUTSIDE", _out)
_small = min_font_audit(fig, 10); print("MIN FONT AUDIT", _small)
assert _small == [], _small
save_figure(fig, "Figure_7")
pd.DataFrame(REG).to_csv("Figure_7_number_registry.csv", index=False)
with pd.ExcelWriter("Figure_7_source_data.xlsx") as xw:
    for k_, v_ in SHEETS.items(): v_.to_excel(xw, sheet_name=k_[:31], index=False)
print("saved Figure_7; registry rows", len(REG), "mismatches", sum(str(r_["check"]).startswith("MISMATCH") for r_ in REG))
