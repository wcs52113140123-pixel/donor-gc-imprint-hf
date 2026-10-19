"""Figure 3 (R2) - Donor control hearts can carry a GC-like imprint whose strength and direction vary by cohort.
Standalone: reproduces Figure_3.{pdf,tif,png}, Figure_3_source_data.xlsx and Figure_3_number_registry.csv
from the input artifacts (paths below). Display-level derivations only; no new tests.
"""
import os, tarfile, numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.colors import to_rgb

import os  # v0.5 package: inputs resolved relative to this script (override with JBTS_INPUTS)
HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.environ.get("JBTS_INPUTS", os.path.join(HERE, "inputs"))
STYLE = os.path.join(HERE, "jbts_style.py")
TAR_FI = os.path.join(INP, "figure_inputs_v05.tar.gz")
TAR_PD = os.path.join(INP, "figure_inputs_v05_perdonor.tar.gz")
XLS = os.path.join(INP, "Supplemental_Tables_S1-S14.xlsx")

exec(open(STYLE).read())
FAM = apply_jbts_style()
mpl.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": FAM, "mathtext.it": FAM + ":italic",
                     "mathtext.bf": FAM + ":bold", "mathtext.default": "regular", "hatch.linewidth": 0.4})
FS_T, FS_S, FS_G = 10, 9, 8          # 3-size ladder: titles/axis labels, ticks/stats/legend, gene/track labels
os.makedirs("_inputs", exist_ok=True)
for t in (TAR_FI, TAR_PD):
    root = os.path.basename(t).split("_", 1)[1].replace(".tar.gz", "")
    if not os.path.isdir(f"_inputs/{root}"):
        with tarfile.open(t) as tf:
            tf.extractall("_inputs")
PD = "_inputs/figure_inputs_v05_perdonor/v13_rerun/"
FI = "_inputs/figure_inputs_v05/"

# ------------------------------------------------------------------ helpers
FIG_W, FIG_H = 183.0, 225.0
REG, SRC = [], {}
CT = ["Cardiomyocyte", "Fibroblast", "Endothelial", "Pericyte", "Myeloid"]
CTS = {"Cardiomyocyte": "CM", "Fibroblast": "Fib", "Endothelial": "EC", "Pericyte": "Peri", "Myeloid": "Mye"}
STUDIES = ["Chaffin 2022", "Koenig 2022", "Reichart 2022"]
SNAME = {s: s.split()[0] for s in STUDIES}


def m(v, nd=2, sign=True):
    s = f"{v:+.{nd}f}" if sign else f"{v:.{nd}f}"
    return s.replace("-", "\u2212")


def fmt_p(p):
    if p < 1e-3:
        e = int(np.floor(np.log10(p))); mant = p / 10 ** e
        if round(mant, 1) >= 10: mant, e = mant / 10, e + 1
        return f"{mant:.1f}\u00d710$^{{{e}}}$".replace("-", "\u2212")
    if p >= 0.1: return f"{p:.2f}"
    d = int(-np.floor(np.log10(p))) + 1          # 2 significant figures, trailing zero kept (0.0030)
    return f"{p:.{d}f}"


from matplotlib.textpath import TextToPath as _TTP
from matplotlib.font_manager import FontProperties as _FP
_ttp = _TTP()
def _tw(s, size, weight="normal"):
    w, h, d = _ttp.get_text_width_height_descent(s, _FP(family=FAM, size=size, weight=weight), ismath=False)
    return w


EXP_FS = 8   # superscript exponents >= 7 pt at final size (deck rule)
def ptext(ax, x, y, p, fs, ha="center", va="center", transform=None, **kw):
    """P value with a manually set exponent (mathtext would shrink it to 0.7 x fs)."""
    tr = transform or ax.transData
    if p >= 1e-3:
        return ax.text(x, y, fmt_p(p), fontsize=fs, ha=ha, va=va, transform=tr, **kw)
    e = int(np.floor(np.log10(p))); mant = p / 10 ** e
    if round(mant, 1) >= 10: mant, e = mant / 10, e + 1
    ms = f"{mant:.1f}\u00d710"; es = str(e).replace("-", "\u2212")
    wm, we = _tw(ms, fs), _tw(es, EXP_FS); tot = wm + we + 0.3
    st = {"center": -tot / 2, "right": -tot, "left": 0.0}[ha]
    ax.annotate(ms, (x, y), xycoords=tr, xytext=(st, 0), textcoords="offset points", fontsize=fs, ha="left", va=va, **kw)
    ax.annotate(es, (x, y), xycoords=tr, xytext=(st + wm + 0.3, fs * 0.34), textcoords="offset points", fontsize=EXP_FS,
                ha="left", va=va, **kw)


def save_exact(fig, stem):
    """Exact canvas (183 mm wide), no tight cropping: PDF, PNG preview (300 dpi), TIFF 600 dpi LZW."""
    fig.savefig(f"{stem}.pdf")
    fig.savefig(f"{stem}.png", dpi=300)
    fig.savefig(f"{stem}.tif", dpi=600, pil_kwargs={"compression": "tiff_lzw"})


def reg(panel, text, value, source, computation, status, lock):
    REG.append(dict(figure="Figure 3", panel=panel, printed_text=text.replace("$", ""), value=value,
                    source_file=source, computation=computation, lock_status=status, matches_locked_text=lock))


def ax_mm(x, y, w, h, **kw):
    """axes from mm coordinates measured from the top-left corner of the canvas"""
    return fig.add_axes([x / FIG_W, 1 - (y + h) / FIG_H, w / FIG_W, h / FIG_H], **kw)


def panel_frame(letter, x, y, w, h, title, tag=None):
    fr = ax_mm(x, y + 6, w, h - 6); fr.axis("off")
    panel_label(fr, letter, x=0, y=1.0)
    fr.text(6.2 / w, 1.0 + 0.4 / (h - 6), title, transform=fr.transAxes, fontsize=FS_T, fontweight="bold",
            va="bottom", ha="left")
    if tag:
        fr.text(1 + 1.4 / w, 1.0 + 0.4 / (h - 6), tag, transform=fr.transAxes, fontsize=FS_G, ha="right", va="bottom",
                **POSTHOC_TAG)
    return fr


def badge(ax, x, y, text, passed, fs=FS_S, transform=None):
    tr = transform or ax.transAxes
    if passed:
        ax.text(x, y, text, transform=tr, fontsize=fs, color="white", fontweight="bold", ha="left", va="center",
                bbox=dict(boxstyle="round,pad=0.25,rounding_size=0.6", fc="black", ec="black", lw=0.6))
    else:
        ax.text(x, y, "\u00d7 " + text, transform=tr, fontsize=fs, color="black", fontweight="bold", ha="left",
                va="center", bbox=dict(boxstyle="round,pad=0.25,rounding_size=0.6", fc="white", ec="black", lw=0.8))


def zero_line(ax, orient="v", v=0):
    (ax.axvline if orient == "v" else ax.axhline)(v, color="#7F7F7F", lw=0.6, ls=(0, (3, 2)), zorder=0)


GCDIV = diverging_cmap(); GCDIV.set_bad(GC_CLASS["not_measured"])
CAP = 3.0
rng = np.random.default_rng(20261004)   # jitter only (display)

# ------------------------------------------------------------------ data
s8 = pd.read_excel(XLS, "S8a_v13_score_tests")
s7 = pd.read_excel(XLS, "S7a_G1")
s14 = pd.read_excel(XLS, "S14a_revision_models")
rp = pd.read_csv(PD + "rev_posthoc_models.csv")
scp = pd.read_csv(PD + "perdonor_scores_SCP3342_HFpEF.csv")
hahn = pd.read_csv(PD + "perdonor_scores_Hahn2021_RVS_HFpEF.csv")
mag = pd.read_csv(PD + "perdonor_scores_MAGNet_exclOverlap_DCM.csv")
gse = pd.read_csv(PD + "perdonor_scores_GSE290367_DCM.csv")
ge = pd.read_csv(PD + "perdonor_gene_expression_GC8_clock10.csv.gz")
_HAHN_CLIN = FI + "02_data/raw/Hahn2021_HFpEF_EMB_bulk/ClinicalCharacteristicsRNAseqSamples.csv"
if not os.path.exists(_HAHN_CLIN):
    _alt = os.path.join(os.environ.get("HAHN2021_DIR", ""), "ClinicalCharacteristicsRNAseqSamples.csv")
    if os.environ.get("HAHN2021_DIR") and os.path.exists(_alt):
        _HAHN_CLIN = _alt
    else:
        raise FileNotFoundError("Figure 3H needs ClinicalCharacteristicsRNAseqSamples.csv from the supplementary files of "
                                "Hahn et al. 2021 (not redistributed here). Download it from the original publication and set "
                                "HAHN2021_DIR to its folder; see PUBLIC_RELEASE_NOTES.md.")
clin = pd.read_csv(_HAHN_CLIN)

fig = plt.figure(figsize=(FIG_W * MM, FIG_H * MM))
CW = (FIG_W - 3 - 11 * 3) / 12.0
def cx(c): return 1.5 + c * (CW + 3)
def cw(n): return n * CW + (n - 1) * 3
ROW_Y = [0, 60, 111, 163]; ROW_H = [60, 51, 52, 62]

# ================================================================== A  forest/table
x0, y0, w, h = cx(0), ROW_Y[0], cw(12), ROW_H[0]
frA = panel_frame("A", x0, y0, w, h, "Control \u2212 case GC score in every contrast")
frA.text(1, 1.0 + 0.4 / (h - 6), "magnitudes are not comparable between cohorts", transform=frA.transAxes, fontsize=FS_G,
         ha="right", va="bottom", style="italic", color="#555555")
LW_, TW_ = 35.0, 60.0
rows = []
def ck(ds, ct, test="C/K"):
    r = s8[(s8.dataset == ds) & (s8.test == test) & (s8.cell_type == ct) & (s8.module == "GC3") & (s8.adjustment == "base")].iloc[0]
    return dict(est=-r.beta, lo=-r.ci_hi, hi=-r.ci_lo, p=r.p, n_case=int(r.n_case), n_ctrl=int(r.n - r.n_case))
def g1(ds, ct):
    r = s7[(s7.dataset == ds) & (s7.cell_type == ct)].iloc[0]
    return dict(g_est=r.ctrl_minus_case, g_lo=r.ci_lo, g_hi=r.ci_hi, g_ph=r.p_holm)
for ct in CT:
    rows.append(dict(cohort="SCP3342", ct=CTS[ct], **ck("SCP3342_HFpEF", ct), **g1("SCP3342_HFpEF", ct),
                     src="XL:S8a C/K GC3 base (-beta); S7a G1"))
r = rp[rp.test == "R2_GC3_case_adj_age"].iloc[0]
rows.append(dict(cohort="Hahn", ct="bulk", est=r.beta, lo=r.ci_lo, hi=r.ci_hi, p=r.p,
                 n_case=int((hahn.case == 1).sum()), n_ctrl=int((hahn.case == 0).sum()),
                 src="rev_posthoc_models R2_GC3_case_adj_age (already ctrl-case)"))
rows.append(dict(cohort="MAGNet", ct="bulk", **ck("MAGNet_exclOverlap_DCM", "bulk", "C/K_bulk"), src="XL:S8a C/K_bulk GC3 base (-beta)"))
for ct in CT:
    rows.append(dict(cohort="GSE290367", ct=CTS[ct], **ck("GSE290367_DCM", ct), **g1("GSE290367_DCM", ct),
                     src="XL:S8a C/K GC3 base (-beta); S7a G1"))
A = pd.DataFrame(rows)
ypos, yy, prev = [], 0.0, None
for c in A.cohort:
    if prev is not None and c != prev: yy += 0.45
    ypos.append(yy); yy += 1.0; prev = c
A["y"] = ypos
ytop, ybot = -0.75, A.y.max() + 0.75
ax = ax_mm(x0 + LW_ + 1, y0 + 11.0, w - LW_ - TW_ - 3, h - 11.0 - 8.5)
axl = ax_mm(x0, y0 + 11.0, LW_, h - 11.0 - 8.5); axl.axis("off")
axt = ax_mm(x0 + w - TW_, y0 + 11.0, TW_, h - 11.0 - 8.5); axt.axis("off")
for a in (ax, axl, axt): a.set_ylim(ybot, ytop)
axl.set_xlim(0, 1); axt.set_xlim(0, 1)
zero_line(ax)
for _, r in A.iterrows():
    col = COHORT[r.cohort]; mk = COHORT_MARKER[r.cohort]
    ax.plot([r.lo, r.hi], [r.y, r.y], color=col, lw=0.9, solid_capstyle="butt", zorder=2)
    ax.plot(r.est, r.y, marker=mk, ms=5.0 if mk != "^" else 5.5, mfc=col, mec="#333333" if r.cohort == "Hahn" else col,
            mew=0.5, ls="none", zorder=3)
    if not np.isnan(r.get("g_est", np.nan)):
        ax.plot([r.g_lo, r.g_hi], [r.y + 0.36] * 2, color="#8C8C8C", lw=0.6, zorder=2)
        ax.plot(r.g_est, r.y + 0.36, marker="D", ms=2.8, mfc="#6E6E6E", mec="#6E6E6E", ls="none", zorder=3)
    axl.text(0.76, r.y, r.ct, fontsize=FS_S, va="center", ha="left")
    axt.text(0.09, r.y, f"{r.n_ctrl}/{r.n_case}", fontsize=FS_S, va="center", ha="center")
    axt.text(0.50, r.y, f"{m(r.est)} ({m(r.lo, sign=False)} to {m(r.hi, sign=False)})", fontsize=FS_S, va="center", ha="center")
    ptext(axt, 0.89, r.y, r.p, FS_S)
    est_T = (r.cohort in ("MAGNet", "Hahn") or (r.cohort == "GSE290367" and r.ct == "CM") or (r.cohort == "SCP3342" and r.ct in ("EC", "Fib")))
    ci_T = r.cohort in ("MAGNet", "Hahn")
    reg("A", m(r.est), round(r.est, 4), r.src, f"control - case GC3 estimate ({r.cohort} {r.ct})", "T" if est_T else "S", "yes" if est_T else "new")
    reg("A", f"({m(r.lo, sign=False)} to {m(r.hi, sign=False)})", f"{r.lo:.4f}..{r.hi:.4f}", r.src, f"95% CI ({r.cohort} {r.ct})",
        "T" if ci_T else "S", "yes" if ci_T else "new")
    reg("A", fmt_p(r.p), r.p, r.src, f"nominal P (HC3) ({r.cohort} {r.ct})",
        "T" if (r.cohort == "GSE290367" and r.ct == "CM") else "S",
        "yes" if (r.cohort == "GSE290367" and r.ct == "CM") else "new")
    _nl = r.cohort in ("MAGNet", "Hahn") or (r.cohort == "SCP3342" and r.ct == "CM")
    reg("A", f"{r.n_ctrl}/{r.n_case}", f"{r.n_ctrl}/{r.n_case}", r.src, "n control / n case", "T" if _nl else "S", "yes" if _nl else "new")
# cohort labels and disease brackets
for coh, lab in [("SCP3342", "SCP3342"), ("Hahn", "Hahn biopsy"), ("MAGNet", "MAGNet"), ("GSE290367", "GSE290367")]:
    ys = A.y[A.cohort == coh]
    axl.text(0.135, ys.mean(), lab, fontsize=FS_S, fontweight="bold", color=COHORT[coh] if coh != "Hahn" else TEXT,
             va="center", ha="left")
for (c0, c1), lab in [(("SCP3342", "Hahn"), "HFpEF"), (("MAGNet", "GSE290367"), "DCM")]:
    ya, yb = A.y[A.cohort == c0].min() - 0.35, A.y[A.cohort == c1].max() + 0.35
    axl.plot([0.10, 0.08, 0.08, 0.10], [ya, ya, yb, yb], color=TEXT, lw=0.6, clip_on=False)
    axl.text(0.035, (ya + yb) / 2, lab, rotation=90, fontsize=FS_S, va="center", ha="center")
yh = A.y[A.cohort == "Hahn"].iloc[0]
rh = A[A.cohort == "Hahn"].iloc[0]
ax.text(rh.hi + 0.25, yh, "post hoc, age + sex", fontsize=FS_G, va="center", ha="left", **POSTHOC_TAG)
# column headers
for xx, t in [(0.09, "ctrl/case"), (0.50, "SD (95% CI)"), (0.89, "P")]:
    axt.text(xx, ytop - 0.05, t, fontsize=FS_S, fontweight="bold", ha="center", va="bottom")
axl.text(0.76, ytop - 0.05, "Cell type", fontsize=FS_S, fontweight="bold", ha="left", va="bottom")
ax.set_xlim(-2.6, 9.2); ax.set_xticks(range(-2, 10, 2))
ax.set_xticklabels([m(v, 0, sign=False) for v in range(-2, 10, 2)])
ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xlabel("Control \u2212 case GC score (control-group SD)", fontsize=FS_T, labelpad=1.5)
# separators between cohorts
for coh in ["Hahn", "MAGNet", "GSE290367"]:
    yb = A.y[A.cohort == coh].min() - 0.72
    for a, x_a, x_b in ((ax, 0, 1), (axt, 0, 1), (axl, 0.13, 1)):
        a.plot([x_a, x_b], [yb, yb], transform=a.get_yaxis_transform(), color=GRID, lw=0.6, clip_on=False)
# key + G1 badge (GSE block, empty right side of the forest)
yg = A.y[A.cohort == "GSE290367"]
badge(ax, 3.55, yg.iloc[0] - 0.2, "G1 failed", False, transform=ax.transData)
ax.plot(3.75, yg.iloc[1] - 0.02, marker="D", ms=2.8, mfc="#6E6E6E", mec="#6E6E6E", ls="none")
ax.text(4.1, yg.iloc[1] - 0.02, "8-gene score (test G1):", fontsize=FS_G, va="center", ha="left")
ax.text(3.55, yg.iloc[2] - 0.38, "SCP3342 +3.6 to +7.7;\nGSE290367 Fib only, +1.68\n(Holm P=0.0087)",
        fontsize=FS_G, va="top", ha="left", linespacing=0.98)
# legend for the two score versions
yk_ = ytop - 0.42
ax.plot(-1.4, yk_, marker="o", ms=4.5, mfc="#555555", mec="#555555", mew=0.6, ls="none", clip_on=False)
ax.text(-1.1, yk_, "3-gene GC score", fontsize=FS_G, va="center", ha="left")
ax.plot(3.2, yk_, marker="D", ms=2.8, mfc="#6E6E6E", mec="#6E6E6E", ls="none", clip_on=False)
ax.text(3.5, yk_, "8-gene GC response score", fontsize=FS_G, va="center", ha="left")
sc = A[A.cohort == "SCP3342"]; gs = A[A.cohort == "GSE290367"]
reg("A", "+3.6 to +7.7", f"{sc.g_est.min():.2f}..{sc.g_est.max():.2f}", "XL:S7a_G1 SCP3342 rows", "range ctrl_minus_case 8-gene", "T", "yes")
gf = gs[gs.ct == "Fib"].iloc[0]
reg("A", "+1.68", round(gf.g_est, 4), "XL:S7a_G1 GSE290367 Fibroblast", "ctrl_minus_case 8-gene", "T", "yes")
reg("A", "Holm P=0.0087", gf.g_ph, "XL:S7a_G1 GSE290367 Fibroblast", "p_holm", "T", "yes")
SRC["3A"] = A.drop(columns=["y"]).rename(columns={"est": "ctrl_minus_case", "lo": "ci_lo", "hi": "ci_hi",
                                                  "g_est": "G1_8gene_ctrl_minus_case", "g_lo": "G1_ci_lo", "g_hi": "G1_ci_hi", "g_ph": "G1_p_holm"})

# ================================================================== B  SCP3342 donor x cell type tiles
x0, y0, w, h = cx(0), ROW_Y[1], cw(8), ROW_H[1]
panel_frame("B", x0, y0, w, h, "SCP3342: every HFpEF donor is below the control mean")
S = scp[scp.cell_type.isin(CT)].copy()
piv = S.pivot(index="unit_id", columns="cell_type", values="score_GC3")[CT]
meta = scp.drop_duplicates("unit_id").set_index("unit_id")[["case", "pool"]]
meta["pool_n"] = meta.pool.str.replace("pool", "").astype(int)
meta["mean_gc"] = piv.mean(axis=1)
meta = meta.sort_values(["case", "pool_n", "mean_gc"], ascending=[True, True, False])
piv = piv.loc[meta.index]
nd = len(meta)
lx, mw = 15.0, w - 15.0 - 17.0
axh = ax_mm(x0 + lx, y0 + 6 + 7.2, mw, 29.0)
axg = ax_mm(x0 + lx, y0 + 6 + 1.0, mw, 2.4)
axp = ax_mm(x0 + lx, y0 + 6 + 3.9, mw, 2.4)
Z = piv.T.values
im = axh.imshow(np.clip(Z, -CAP, CAP), cmap=GCDIV, vmin=-CAP, vmax=CAP, aspect="auto", interpolation="nearest")
for (i, j) in zip(*np.where(np.isnan(Z))):
    axh.add_patch(Rectangle((j - .5, i - .5), 1, 1, fc=GC_CLASS["not_measured"], ec="#9A9A9A", hatch="////", lw=0))
axh.set_yticks(range(5)); axh.set_yticklabels([CTS[c] for c in CT], fontsize=FS_S)
axh.set_xticks([]); [s.set_visible(False) for s in axh.spines.values()]
axh.tick_params(axis="y", length=0, pad=2)
for i in range(1, 5): axh.axhline(i - 0.5, color="white", lw=0.8)
ncon = int((meta.case == 0).sum())
for a in (axh, axg, axp): a.axvline(ncon - 0.5, color="white", lw=1.6)
axg.imshow([[to_rgb(CONTROL) if c == 0 else to_rgb(COHORT["SCP3342"]) for c in meta.case]], aspect="auto", interpolation="nearest")
pshade = {p: ("#D9D9D9" if k % 2 == 0 else "#8C8C8C") for k, p in enumerate(sorted(meta.pool_n.unique()))}
axp.imshow([[to_rgb(pshade[p]) for p in meta.pool_n]], aspect="auto", interpolation="nearest")
for a, lab in ((axg, "Group"), (axp, "Pool")):
    a.set_xticks([]); a.set_yticks([0]); a.set_yticklabels([lab], fontsize=FS_G); a.tick_params(axis="y", length=0, pad=2)
    [s.set_visible(False) for s in a.spines.values()]
# pool numbers once per run
runs, start = [], 0
for j in range(1, nd + 1):
    if j == nd or meta.pool_n.iloc[j] != meta.pool_n.iloc[start] or meta.case.iloc[j] != meta.case.iloc[start]:
        runs.append((start, j - 1, meta.pool_n.iloc[start])); start = j
for a0, a1, p in runs:
    axp.text((a0 + a1) / 2, 0, str(p), fontsize=FS_G, ha="center", va="center",
             color="white" if pshade[p] == "#8C8C8C" else TEXT)
# group brackets under the matrix
for (a0, a1, lab, col) in [(0, ncon - 1, f"{ncon} controls", CONTROL), (ncon, nd - 1, f"{nd - ncon} HFpEF", COHORT["SCP3342"])]:
    axh.plot([a0 - .3, a0 - .3, a1 + .3, a1 + .3], [4.75, 4.95, 4.95, 4.75], color=TEXT, lw=0.6, clip_on=False)
    axh.text((a0 + a1) / 2, 5.05, lab + ("  (8 pools, 3 per pool)" if col == CONTROL else "  (2\u20133 per pool)"),
             fontsize=FS_S, ha="center", va="top", color=TEXT)
reg("B", f"{ncon} controls", ncon, "perdonor_scores_SCP3342_HFpEF.csv", "count unit_id case==0", "T", "yes")
reg("B", f"{nd - ncon} HFpEF", nd - ncon, "perdonor_scores_SCP3342_HFpEF.csv", "count unit_id case==1", "T", "yes")
reg("B", "8 pools, 3 per pool", meta.pool.nunique(), "perdonor_scores_SCP3342_HFpEF.csv", "n pools; controls per pool", "S", "new")
reg("B", "2-3 per pool", "2-3", "perdonor_scores_SCP3342_HFpEF.csv", "HFpEF per pool", "S", "new")
cax = ax_mm(x0 + w - 13.5, y0 + 6 + 7.2, 2.2, 29.0)
cb = fig.colorbar(im, cax=cax, ticks=[-3, 0, 3]); cb.outline.set_linewidth(0.4)
cb.ax.yaxis.set_ticklabels(["\u2264\u22123", "0", "\u2265+3"]); cb.ax.tick_params(labelsize=FS_G, length=2)
cb.ax.set_title("GC score,\nSD", fontsize=FS_G, pad=6, loc="left", linespacing=1.0)
assert (piv[meta.case == 1].max() < 0).all() and (np.abs(piv[meta.case == 0].mean()) < 0.05).all(), "3B title must hold for every plotted row"
SRC["3B"] = piv.reset_index().merge(meta.reset_index()[["unit_id", "case", "pool"]], on="unit_id")

# ================================================================== C  bulk cohorts per sample
x0, y0, w, h = cx(8), ROW_Y[1], cw(4), ROW_H[1]
panel_frame("C", x0, y0, w, h, "Bulk cohorts, per sample")
sub = [("Hahn", hahn, "Ctrl", "HFpEF", rp[rp.test == "R2_GC3_case_adj_age"].iloc[0], False, "rev_posthoc_models R2_GC3_case_adj_age"),
       ("MAGNet", mag, "Ctrl", "DCM", s8[(s8.test == "C/K_bulk") & (s8.adjustment == "base")].iloc[0], True, "XL:S8a C/K_bulk base")]
pw = (w - 11.5 - 3) / 2
srcC = []
for k, (coh, d, l0, l1, r, flip, srcn) in enumerate(sub):
    TOPC = 14.0
    a = ax_mm(x0 + 11.5 + k * (pw + 3), y0 + 6 + TOPC, pw, h - 6 - TOPC - 8.5)
    for g, col in ((0, CONTROL), (1, COHORT[coh])):
        v = d.loc[d.case == g, "score_GC3"].values
        a.boxplot(v, positions=[g], widths=0.55, showfliers=False, whis=1.5,
                  medianprops=dict(color=TEXT, lw=1.0), boxprops=dict(color=TEXT, lw=0.6),
                  whiskerprops=dict(color=TEXT, lw=0.6), capprops=dict(color=TEXT, lw=0.6), zorder=1)
        a.scatter(g + rng.uniform(-0.2, 0.2, len(v)), v, s=6, marker="s", fc=col, ec="#555555" if coh == "Hahn" and g == 1 else "none",
                  lw=0.3, alpha=0.85, zorder=2)
        srcC.append(pd.DataFrame(dict(cohort=coh, unit_id=d.loc[d.case == g, "unit_id"].values, case=g, score_GC3=v)))
    est, lo, hi = (-r.beta, -r.ci_hi, -r.ci_lo) if flip else (r.beta, r.ci_lo, r.ci_hi)
    a.text(0.5, 1.03, f"{m(est)} SD\n({m(lo, sign=False)} to {m(hi, sign=False)})", transform=a.transAxes, ha="center", va="bottom",
           fontsize=FS_S, linespacing=1.05)
    a.text(0.5, 1.0 + 8.4 / (h - 6 - TOPC - 8.5), "Hahn biopsy" if coh == "Hahn" else coh, transform=a.transAxes, ha="center",
           va="bottom", fontsize=FS_S, fontweight="bold", color=COHORT[coh] if coh != "Hahn" else TEXT)
    if coh == "Hahn":
        a.text(-0.35, -10.2, "post hoc", ha="left", va="center", fontsize=FS_G, **POSTHOC_TAG)
    reg("C", f"{m(est)} SD ({m(lo, sign=False)} to {m(hi, sign=False)})", round(est, 4), srcn, "control - case GC3", "T", "yes")
    n0, n1 = int((d.case == 0).sum()), int((d.case == 1).sum())
    a.set_xticks([0, 1]); a.set_xticklabels([f"{l0}\n{n0}", f"{l1}\n{n1}"], fontsize=FS_S, linespacing=1.0)
    reg("C", f"{n0} / {n1}", f"{n0}/{n1}", f"perdonor_scores_{coh}", "n per group", "T", "yes")
    a.set_xlim(-0.6, 1.6); a.set_ylim(-11, 4.5); a.set_yticks([-10, -5, 0])
    a.set_yticklabels(["\u221210", "\u22125", "0"])
    zero_line(a, "h")
    if k == 0: a.set_ylabel("GC score (control SD)", fontsize=FS_T, labelpad=1)
    else: a.set_yticklabels([])
SRC["3C"] = pd.concat(srcC)

# ================================================================== D  GSE290367 per study forest
x0, y0, w, h = cx(0), ROW_Y[2], cw(6), ROW_H[2]
panel_frame("D", x0, y0, w, h, "GSE290367: opposite study-level effects", tag="post hoc")
R5 = rp[(rp.test == "R5_per_study") & (rp.module == "GC3")].copy()
R5[["ct", "study"]] = R5.cell_type.str.split("|", expand=True)
ad = ax_mm(x0 + 13, y0 + 6 + 1.0, 46, h - 6 - 1.0 - 8.5)
off = {"Chaffin 2022": -0.27, "Koenig 2022": 0.0, "Reichart 2022": 0.27}
zero_line(ad)
for i, ct in enumerate(CT):
    pr = ck("GSE290367_DCM", ct)
    ad.add_patch(Rectangle((pr["lo"], i - 0.42), pr["hi"] - pr["lo"], 0.84, fc="#E3E3E3", ec="none", zorder=0))
    ad.plot([pr["est"]] * 2, [i - 0.42, i + 0.42], color="#8C8C8C", lw=1.0, zorder=1)
    for s in STUDIES:
        r = R5[(R5.ct == ct) & (R5.study == s)].iloc[0]
        col = GSE_STUDY[SNAME[s]]; yy = i + off[s]
        ad.plot([r.ci_lo, r.ci_hi], [yy, yy], color=col, lw=0.9, zorder=2)
        ad.plot(r.beta, yy, marker="o", ms=4.6, mfc=col, mec="#1F4F4A" if s == "Reichart 2022" else col,
                mew=0.4, ls="none", zorder=3)
        if s == "Chaffin 2022" and ct in ("Myeloid", "Fibroblast", "Pericyte"):
            ad.text(r.ci_hi + 0.15, yy, m(r.beta), fontsize=FS_G, va="center", ha="left", color=col)
            reg("D", m(r.beta), round(r.beta, 4), "rev_posthoc_models R5_per_study GC3", "control - DCM, age+sex adj", "T", "yes")
ad.set_ylim(4.6, -0.6); ad.set_yticks(range(5)); ad.set_yticklabels([CTS[c] for c in CT], fontsize=FS_S)
ad.tick_params(axis="y", length=0)
ad.set_xlim(-5.2, 7.6); ad.set_xticks([-4, -2, 0, 2, 4, 6]); ad.set_xticklabels([m(v, 0, sign=False) for v in [-4, -2, 0, 2, 4, 6]])
ad.set_xlabel("Control \u2212 DCM GC score (SD, 95% CI)", fontsize=FS_T, labelpad=1.5)
# key with per-study counts
cnt = gse[gse.cell_type == "Fibroblast"].groupby(["study", "case"]).size().unstack()
cntR = gse[gse.study == "Reichart 2022"].groupby("cell_type").case.apply(lambda s: int((s == 1).sum()))
reich = R5[R5.study == "Reichart 2022"]
kx = 0.0; ky = y0 + 6 + 1.0
ak = ax_mm(x0 + 61, ky, w - 61, h - 6 - 9.5); ak.axis("off"); ak.set_xlim(0, 1); ak.set_ylim(1, 0)
keys = [("Chaffin 2022", f"Chaffin {int(cnt.loc['Chaffin 2022', 0])}/{int(cnt.loc['Chaffin 2022', 1])}",
         "Ctrl > DCM in\nMye, Fib, Peri"),
        ("Koenig 2022", f"Koenig {int(cnt.loc['Koenig 2022', 0])}/{int(cnt.loc['Koenig 2022', 1])}", "No cell type differs"),
        ("Reichart 2022", f"Reichart {int(cnt.loc['Reichart 2022', 0])}/{cntR.min()}\u2013{cntR.max()}",
         f"DCM > Ctrl in all 5,\n{m(-reich.beta.max(), sign=False)}\u2013{m(-reich.beta.min(), sign=False)} SD,\nall P\u2264{np.ceil(reich.p.max() * 1000) / 1000:.3f}")]
KH = h - 6 - 9.5            # key axes height in mm
LINE = 3.2 / KH             # one 8-pt line in axes units
ak.text(0.0, 0.0, "Study (n ctrl/DCM)", fontsize=FS_G, va="top", ha="left", color="#555555")
yk = LINE * 1.25
for s, l1, l2 in keys:
    ak.plot(0.05, yk + LINE * 0.45, marker="o", ms=4.6, mfc=GSE_STUDY[SNAME[s]],
            mec="#1F4F4A" if s == "Reichart 2022" else GSE_STUDY[SNAME[s]], mew=0.4, ls="none", clip_on=False)
    ak.text(0.13, yk, l1, fontsize=FS_G, va="top", ha="left", fontweight="bold")
    ak.text(0.13, yk + LINE, l2, fontsize=FS_G, va="top", ha="left", linespacing=1.1)
    yk += LINE * (1 + l2.count("\n") + 1) + LINE * 0.35
ak.add_patch(Rectangle((0.0, yk + LINE * 0.15), 0.10, LINE * 0.7, fc="#E3E3E3", ec="none", clip_on=False))
ak.plot([0.05, 0.05], [yk + LINE * 0.15, yk + LINE * 0.85], color="#8C8C8C", lw=1.0, clip_on=False)
ak.text(0.13, yk, "Pooled, study-\nadjusted (as in A)", fontsize=FS_G, va="top", ha="left", linespacing=1.1)
for s in STUDIES:
    reg("D", f"{SNAME[s]} {int(cnt.loc[s, 0])} ctrl/{int(cnt.loc[s, 1]) if s != 'Reichart 2022' else str(cntR.min()) + '-' + str(cntR.max())} DCM",
        f"{int(cnt.loc[s, 0])}/{int(cnt.loc[s, 1]) if s != 'Reichart 2022' else str(cntR.min()) + '-' + str(cntR.max())}",
        "perdonor_scores_GSE290367_DCM.csv", "per-study donor counts (printed as n ctrl/DCM)", "N", "new")
reg("D", f"{m(-reich.beta.max(), sign=False)}-{m(-reich.beta.min(), sign=False)} SD", f"{-reich.beta.max():.4f}..{-reich.beta.min():.4f}",
    "rev_posthoc_models R5 Reichart", "range DCM - control", "T", "yes")
reg("D", f"all P<={np.ceil(reich.p.max() * 1000) / 1000:.3f}", reich.p.max(), "rev_posthoc_models R5 Reichart", "max P", "T", "yes")
ko = R5[R5.study == "Koenig 2022"]
assert (ko.p > 0.05).all()
SRC["3D"] = pd.concat([R5[["ct", "study", "n", "beta", "ci_lo", "ci_hi", "p"]].rename(columns={"beta": "ctrl_minus_DCM"}),
                       pd.DataFrame([dict(ct=ct, study="pooled_study_adjusted", n=None, ctrl_minus_DCM=ck("GSE290367_DCM", ct)["est"],
                                          ci_lo=ck("GSE290367_DCM", ct)["lo"], ci_hi=ck("GSE290367_DCM", ct)["hi"], p=ck("GSE290367_DCM", ct)["p"]) for ct in CT])])

# ================================================================== E  GSE290367 per donor by study
x0, y0, w, h = cx(6), ROW_Y[2], cw(6), ROW_H[2]
panel_frame("E", x0, y0, w, h, "GSE290367 per donor, by study", tag="post hoc")
fw = (w - 11 - 4 * 1.6) / 5
for i, ct in enumerate(CT):
    a = ax_mm(x0 + 11 + i * (fw + 1.6), y0 + 6 + 4.5, fw, h - 6 - 4.5 - 13.5)
    d = gse[gse.cell_type == ct]
    for j, s_ in enumerate(STUDIES):
        for g, col, dx in ((0, CONTROL, -0.19), (1, GSE_STUDY[SNAME[s_]], 0.19)):
            v = d.loc[(d.study == s_) & (d.case == g), "score_GC3"].values
            a.scatter(j + dx + rng.uniform(-0.09, 0.09, len(v)), v, s=4.5, marker="o", fc=col,
                      ec="#1F4F4A" if (g == 1 and s_ == "Reichart 2022") else "none", lw=0.25, alpha=0.9, zorder=2, clip_on=False)
            a.plot([j + dx - 0.15, j + dx + 0.15], [np.median(v)] * 2, color=TEXT, lw=1.0, zorder=3)
    zero_line(a, "h")
    a.set_xlim(-0.55, 2.55); a.set_ylim(-10, 6)
    a.set_xticks([0, 1, 2]); a.set_xticklabels(["C", "K", "R"], fontsize=FS_S)
    a.set_title(CTS[ct], fontsize=FS_S, pad=1.5)
    a.set_yticks([-8, -4, 0, 4])
    if i == 0:
        a.set_yticklabels(["\u22128", "\u22124", "0", "4"]); a.set_ylabel("GC score (control SD)", fontsize=FS_T, labelpad=1)
    else:
        a.set_yticklabels([])
fig.text((x0 + w) / FIG_W, 1 - (y0 + h - 3.6) / FIG_H, "C, Chaffin; K, Koenig; R, Reichart; bar, median", fontsize=FS_G,
         ha="right", va="bottom", color="#555555")
fig.text((x0 + w) / FIG_W, 1 - (y0 + h - 0.3) / FIG_H, "left of each pair, ctrl (gray); right, DCM (study shade, as in D)",
         fontsize=FS_G, ha="right", va="bottom", color="#555555")
SRC["3E"] = gse[["cell_type", "unit_id", "study", "case", "age", "sex", "score_GC3"]]

# ================================================================== F  GC target genes heatmap
x0, y0, w, h = cx(0), ROW_Y[3], cw(7), ROW_H[3]
panel_frame("F", x0, y0, w, h, "GC target genes, per donor")
GENES = ["FKBP5", "TSC22D3", "ZBTB16", "PER1", "KLF9", "DDIT4", "KLF15", "MT2A"]
blocks = [("SCP3342", "SCP3342_HFpEF", "Fibroblast", scp, "SCP3342 Fib", "HFpEF"),
          ("Hahn", "Hahn2021_RVS_HFpEF", "bulk", hahn, "Hahn", "HFpEF"),
          ("MAGNet", "MAGNet_exclOverlap_DCM", "bulk", mag, "MAGNet", "DCM"),
          ("GSE290367", "GSE290367_DCM", "Myeloid", gse, "GSE Mye", "DCM")]
gl, gr = 21.0, 12.5
hw_tot = w - gl - gr; gap = 1.8; bw = (hw_tot - 3 * gap) / 4
ytr = y0 + 6 + 8.0
srcF = []
for k, (coh, ds, ctp, sdf, lab, cl) in enumerate(blocks):
    bx = x0 + gl + k * (bw + gap)
    sd = sdf[sdf.cell_type == ctp][["unit_id", "case", "score_GC3"] + (["study"] if coh == "GSE290367" else []) + (["pool"] if coh == "SCP3342" else [])]
    if coh == "GSE290367":   # ordered by study, then group, then score (display only)
        sd = sd.assign(_s=sd.study.map({s_: i for i, s_ in enumerate(STUDIES)})).sort_values(["_s", "case", "score_GC3"],
                                                                                            ascending=[True, True, False]).drop(columns="_s")
    else:
        sd = sd.sort_values(["case", "score_GC3"], ascending=[True, False])
    sub = ge[(ge.dataset == ds) & (ge.cell_type == ctp)].pivot(index="gene", columns="unit_id", values="z_ref")
    sub.columns = sub.columns.astype(str); sd["unit_id"] = sd.unit_id.astype(str)
    Z = sub.loc[GENES, sd.unit_id].values
    ah = ax_mm(bx, ytr + 5.4, bw, h - 6 - 8.0 - 5.4 - 3.0)
    im2 = ah.imshow(np.clip(Z, -CAP, CAP), cmap=GCDIV, vmin=-CAP, vmax=CAP, aspect="auto", interpolation="nearest")
    nanrows = [i for i in range(len(GENES)) if np.isnan(Z[i]).all()]
    for i in range(len(GENES)):
        if np.isnan(Z[i]).any():
            js = np.where(np.isnan(Z[i]))[0]
            # contiguous NaN cells -> hatched rectangles
            st = js[0]; pv = js[0]
            for jj in list(js[1:]) + [None]:
                if jj is None or jj != pv + 1:
                    ah.add_patch(Rectangle((st - .5, i - .5), pv - st + 1, 1, fc=GC_CLASS["not_measured"], ec="#9A9A9A",
                                           hatch="//////", lw=0))
                    if jj is not None: st = jj
                if jj is not None: pv = jj
    ah.axhline(2.5, color="white", lw=1.6)
    ah.set_xticks([]); [s.set_visible(False) for s in ah.spines.values()]
    if k == 0:
        ah.set_yticks(range(len(GENES))); ah.set_yticklabels(GENES, fontsize=FS_G, fontstyle="italic"); ah.tick_params(axis="y", length=0, pad=1.5)
    else:
        ah.set_yticks([])
    nc = int((sd.case == 0).sum())
    if coh == "GSE290367":
        cuts = list(np.where(sd.study.values[1:] != sd.study.values[:-1])[0] + 0.5)
    else:
        cuts = [nc - 0.5]
    for cut in cuts: ah.axvline(cut, color="white", lw=1.2)
    # tracks
    ag = ax_mm(bx, ytr, bw, 2.2)
    ag.imshow([[to_rgb(CONTROL) if c == 0 else to_rgb(COHORT[coh]) for c in sd.case]], aspect="auto", interpolation="nearest")
    for cut in cuts: ag.axvline(cut, color="white", lw=1.2)
    at = ax_mm(bx, ytr + 2.7, bw, 2.2)
    if coh == "GSE290367":
        at.imshow([[to_rgb(GSE_STUDY[SNAME[s]]) for s in sd.study]], aspect="auto", interpolation="nearest")
    elif coh == "SCP3342":
        pn = sd.pool.str.replace("pool", "").astype(int)
        at.imshow([[to_rgb(pshade[p]) for p in pn]], aspect="auto", interpolation="nearest")
    else:
        at.axis("off")
    for cut in cuts: at.axvline(cut, color="white", lw=1.2)
    for a in (ag, at):
        a.set_xticks([]); a.set_yticks([]); [s.set_visible(False) for s in a.spines.values()]
    if k == 0:
        for a, l in ((ag, "Group"), (at, "Pool/study")):
            a.set_yticks([0]); a.set_yticklabels([l], fontsize=FS_G); a.tick_params(axis="y", length=0, pad=1.5)
    fig.text((bx + bw / 2) / FIG_W, 1 - (ytr - 0.6) / FIG_H, lab.replace("\n", " "), fontsize=FS_S,
             ha="center", va="bottom", fontweight="bold", color=COHORT[coh] if coh != "Hahn" else TEXT)
    fig.text((bx + bw / 2) / FIG_W, 1 - (y0 + h - 0.3) / FIG_H, f"{nc}/{len(sd) - nc}",
             fontsize=FS_G, ha="center", va="bottom")
    reg("F", f"{nc}/{len(sd) - nc}", f"{nc}/{len(sd) - nc}", ds, "n per group in heatmap block", "T" if coh != "GSE290367" else "S",
        "yes" if coh != "GSE290367" else "new")
    long = sub.loc[GENES, sd.unit_id].reset_index().melt(id_vars="gene", var_name="unit_id", value_name="z_ref")
    long = long.merge(sd, on="unit_id"); long.insert(0, "block", f"{coh}_{ctp}")
    long["column_order"] = long.unit_id.map({u: i for i, u in enumerate(sd.unit_id)}); srcF.append(long)
fig.text((x0 + gl - 1.5) / FIG_W, 1 - (y0 + h - 0.3) / FIG_H, "n ctrl/case", fontsize=FS_G, ha="right", va="bottom", style="italic")
# bracket for score genes (left)
ahs = ax_mm(x0, ytr + 5.4, 3.0, h - 6 - 8.0 - 5.4 - 3.0); ahs.axis("off"); ahs.set_ylim(7.5, -0.5); ahs.set_xlim(0, 1)
for (a0, a1, lab) in ((0, 2, "GC score"), (3, 7, "8-gene adds")):
    ahs.plot([0.95, 0.75, 0.75, 0.95], [a0 - .35, a0 - .35, a1 + .35, a1 + .35], color=TEXT, lw=0.6, clip_on=False)
    ahs.text(0.35, (a0 + a1) / 2, lab, rotation=90, fontsize=FS_G, ha="center", va="center")
cax = ax_mm(x0 + w - gr + 2.0, ytr + 5.4, 2.2, h - 6 - 8.0 - 5.4 - 3.0)
cb = fig.colorbar(im2, cax=cax, ticks=[-3, 0, 3]); cb.outline.set_linewidth(0.4)
cb.ax.set_yticklabels(["\u2264\u22123", "0", "\u2265+3"]); cb.ax.tick_params(labelsize=FS_G, length=2)
cb.ax.set_title("z, SD", fontsize=FS_G, pad=6, loc="left")
_gm = gse[(gse.cell_type == "Myeloid") & (gse.case == 1)]
reg("F (legend)", f"{int((_gm.study == 'Reichart 2022').sum())} of {len(_gm)} DCM columns are Reichart", f"{int((_gm.study == 'Reichart 2022').sum())}/{len(_gm)}",
    "perdonor_scores_GSE290367_DCM.csv", "count of Reichart DCM donors in the GSE290367 myeloid block", "N", "new")
SRC["3F"] = pd.concat(srcF)

# ================================================================== G  age within controls
GH_ = 30.0
x0, y0, w, h = cx(7), ROW_Y[3], cw(5), GH_
panel_frame("G", x0, y0, w, h, "No detectable age association")
fig.text((x0 + 9.0) / FIG_W, 1 - (y0 + 6.3) / FIG_H, "controls only; post hoc", fontsize=FS_G, ha="left", va="top", **POSTHOC_TAG)
ag_ = ax_mm(x0 + 11.0, y0 + 6 + 4.6, 27.0, h - 6 - 4.6 - 8.0)
gG = [("GSE290367", gse[(gse.cell_type == "Myeloid") & (gse.case == 0)], "^", "GSE290367_DCM", "Myeloid", "GSE290367 Mye"),
      ("MAGNet", mag[mag.case == 0], "s", "MAGNet_exclOverlap_DCM", "bulk", "MAGNet"),
      ("Hahn", hahn[hahn.case == 0], "s", "Hahn2021_RVS_HFpEF", "bulk", "Hahn biopsy")]
TH = h - 6 - 0.2
txt = ax_mm(x0 + 40.5, y0 + 6 + 0.2, w - 40.5, TH); txt.axis("off"); txt.set_ylim(1, 0); txt.set_xlim(0, 1)
L1 = 2.95 / TH
srcG = []; yl = 0.0
R3 = s14[s14.test == "R3_GC3_age_in_reference"]
def per_dec(v):
    t = m(v, 2)
    return "0.00" if t in ("\u22120.00", "+0.00") else t
for coh, d, mk, ds, ctp, lab in gG:
    ag_.scatter(d.age, d.score_GC3, s=5, marker=mk, fc=COHORT[coh], ec="#555555" if coh == "Hahn" else "none", lw=0.3, alpha=0.85, zorder=2)
    txt.text(0.0, yl, f"{lab} ({len(d)})", fontsize=FS_G, fontweight="bold", va="top", color=COHORT[coh] if coh != "Hahn" else TEXT)
    reg("G", f"({len(d)})", len(d), ds, "n controls with age (plotted)", "S", "new")
    if coh == "GSE290367":
        rr = R3[R3.dataset == ds]
        t1 = f"{per_dec(rr.beta.min())} to {per_dec(rr.beta.max())} /decade"
        txt.text(0.0, yl + L1, t1, fontsize=FS_G, va="top")
        txt.text(0.0, yl + 2 * L1, "(range, 5 cell types)", fontsize=FS_G, va="top", color="#555555")
        reg("G", t1, f"{rr.beta.min():.4f}..{rr.beta.max():.4f}", "XL:S14a R3_GC3_age_in_reference GSE290367 (5 cell types)",
            "range of age slopes, SD per decade (S-R7: -0.10 to +0.11 SD per decade)", "S", "yes")
        yl += 3 * L1 + 0.02
    else:
        r = R3[R3.dataset == ds].iloc[0]
        t1 = f"{per_dec(r.beta * 10)} /decade"
        txt.text(0.0, yl + L1, t1, fontsize=FS_G, va="top")
        reg("G", t1, round(r.beta * 10, 5), f"XL:S14a R3_GC3_age_in_reference {ds}",
            "age slope x 10 (SD per decade; S-R7 states it per year: " + ("+0.002/y" if coh == "MAGNet" else "0.000/y") + ")", "N", "new")
        yl += 2 * L1 + 0.02
    srcG.append(d[["unit_id", "age", "score_GC3"]].assign(cohort=coh))
allp = R3.p
pm = f"{allp.min():.2f}"
txt.text(0.0, yl + 0.01, f"min P={pm}", fontsize=FS_G, va="top", fontweight="bold")
reg("G", f"min P={pm}", allp.min(), "XL:S14a R3_GC3_age_in_reference", "min P over 7 rows (text: all P>=0.32)", "T", "yes")
zero_line(ag_, "h")
ag_.set_xlabel("Age (y)", fontsize=FS_T, labelpad=0.5); ag_.set_ylabel("GC score", fontsize=FS_T, labelpad=0.5)
ag_.set_xlim(10, 88); ag_.set_xticks([20, 40, 60, 80]); ag_.set_ylim(-4.5, 4.5); ag_.set_yticks([-3, 0, 3]); ag_.set_yticklabels(["\u22123", "0", "3"])
SRC["3G"] = pd.concat(srcG)

# ================================================================== H  medications within HFpEF
x0, y0, w, h = cx(7), ROW_Y[3] + GH_ + 0.5, cw(5), ROW_H[3] - GH_ - 0.5
panel_frame("H", x0, y0, w, h, "Medications within HFpEF", tag="post hoc")
fig.text((x0) / FIG_W, 1 - (y0 + 6.3) / FIG_H, "Hahn biopsy (41); MRA and steroid use not recorded", fontsize=FS_G, ha="left", va="top",
         **POSTHOC_TAG)
hm = hahn.merge(clin[["SequenceRunId", "ACEi/ARB", "Beta Blocker"]], left_on="unit_id", right_on="SequenceRunId", how="left")
hm = hm[hm.case == 1]
reg("H", "Hahn biopsy (41)", len(hm), "perdonor_scores_Hahn2021_RVS_HFpEF.csv", "n HFpEF", "T", "yes")
ah_ = ax_mm(x0 + 19.0, y0 + 6 + 4.4, 24.0, h - 6 - 4.4 - 8.5)
labs, srcH = [], []
for i, (col_, nm, adj) in enumerate([("ACEi/ARB", "ACEi/ARB", "ACEi/ARB"), ("Beta Blocker", "\u03b2-blocker", "Beta Blocker")]):
    r = s14[(s14.test == "R4_med_within_HFpEF") & (s14.adjustment == adj)].iloc[0]
    for k, lv in enumerate(["No", "Yes"]):
        v = hm.loc[hm[col_] == lv, "score_GC3"].values; yy = i * 2.5 + k
        ah_.scatter(v, yy + rng.uniform(-0.2, 0.2, len(v)), s=5, marker="s", fc=COHORT["Hahn"], ec="#555555", lw=0.3, zorder=2)
        ah_.plot([np.median(v)] * 2, [yy - 0.36, yy + 0.36], color=TEXT, lw=1.0, zorder=3)
        labs.append((yy, f"{nm} no ({len(v)})" if k == 0 else f"yes ({len(v)})"))
        reg("H", f"{lv.lower()} ({len(v)})", len(v), "ClinicalCharacteristicsRNAseqSamples.csv", f"{nm} {lv} count, HFpEF", "N", "new")
        srcH.append(pd.DataFrame(dict(medication=nm, use=lv, unit_id=hm.loc[hm[col_] == lv, "unit_id"].values, score_GC3=v)))
    ah_.text(1.04, i * 2.5 + 0.5, f"yes\u2212no {m(r.beta)}\n({m(r.ci_lo, sign=False)} to {m(r.ci_hi, sign=False)})", transform=ah_.get_yaxis_transform(),
             fontsize=FS_G, va="center", ha="left", linespacing=1.05)
    reg("H", f"{m(r.beta)} ({m(r.ci_lo, sign=False)} to {m(r.ci_hi, sign=False)})", round(r.beta, 4), "XL:S14a R4_med_within_HFpEF",
        "yes - no, GC3 (S14a model)", "T", "yes")
ah_.set_yticks([l[0] for l in labs]); ah_.set_yticklabels([l[1] for l in labs], fontsize=FS_G); ah_.tick_params(axis="y", length=0, pad=1)
ah_.set_ylim(3.75, -0.45); ah_.set_xlim(-8.8, -1.2); ah_.set_xticks([-8, -6, -4, -2])
ah_.set_xticklabels(["\u22128", "\u22126", "\u22124", "\u22122"])
ah_.set_xlabel("GC score (control SD)", fontsize=FS_T, labelpad=0.5)
ah_.spines["left"].set_visible(False)
SRC["3H"] = pd.concat(srcH)

# ------------------------------------------------------------------ save
save_exact(fig, "Figure_3")
# JACC: RGB (not RGBA) TIFF, 600 dpi, LZW
from PIL import Image as _Im
_t = _Im.open("Figure_3.tif"); _t.convert("RGB").save("Figure_3.tif", compression="tiff_lzw", dpi=(600, 600))

pd.DataFrame(REG).to_csv("Figure_3_number_registry.csv", index=False)
with pd.ExcelWriter("Figure_3_source_data.xlsx") as xw:
    for k, d in SRC.items():
        d.to_excel(xw, sheet_name=f"Fig{k}", index=False)
