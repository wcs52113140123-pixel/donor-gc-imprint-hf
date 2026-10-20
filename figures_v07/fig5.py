"""Figure 5 (v0.7, reader-facing labels). Identifiability determines whether the imprint can be separated from disease effects.
v0.7: panel H badge reworded in plain language (no internal test code); panel H heatmap rows 4.3 -> 4.0 mm and
heatmap shifted 2 mm down to make room for the two-line badge; all data, statistics and printed numbers unchanged.
Standalone: rebuilds Figure_5.{pdf,png,tif}, Figure_5_source_data.xlsx, Figure_5_number_registry.csv
from the project input artifacts (paths below). Display-level derivations only (KDE, counts, ranks); no new tests.
"""
import os, tarfile, re
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch, Polygon
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

import os  # public package: inputs resolved relative to this script (override with JBTS_INPUTS)
HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.environ.get("JBTS_INPUTS", os.path.join(HERE, "..", "figures_v05", "inputs"))
STYLE = os.path.join(HERE, "jbts_style.py")  # v0.6: >=10 pt lettering
TAR_FI = os.path.join(INP, "figure_inputs_v05.tar.gz")
TAR_PD = os.path.join(INP, "figure_inputs_v05_perdonor.tar.gz")
XLP = os.path.join(INP, "Supplemental_Tables_S1-S14.xlsx")
WORK = os.environ.get("FIG_DATA_DIR", "data")
for t in (TAR_FI, TAR_PD):
    root = "figure_inputs_v05_perdonor" if "perdonor" in t else "figure_inputs_v05"
    if not os.path.isdir(os.path.join(WORK, root)):
        os.makedirs(WORK, exist_ok=True)
        with tarfile.open(t) as tf:
            tf.extractall(WORK)
FI = f"{WORK}/figure_inputs_v05/04_results/"
PD = f"{WORK}/figure_inputs_v05_perdonor/v13_rerun/"

exec(open(STYLE).read())
FONT = apply_jbts_style()
mpl.rcParams["savefig.pad_inches"] = 0
FS_T, FS_AX, FS_TK, FS_LG, FS_SM = SIZES["panel_title"], SIZES["axis_label"], SIZES["tick"], SIZES["legend"], SIZES["min_gene_label"]

W, H = 183.0, 240.0
fig = plt.figure(figsize=(W * MM, H * MM))
fig.patches.append(Rectangle((0, 0), 1, 1, transform=fig.transFigure, fc="white", ec="none", zorder=-10))


def axmm(x, y, w, h, **kw):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **kw)


def header(letter, x, y, title, posthoc=False, right=None, tag_line=1):
    fig.text(x / W, 1 - y / H, letter, fontsize=SIZES["panel_letter"], fontweight="bold", va="top", ha="left")
    fig.text((x + 5.6) / W, 1 - (y + 0.7) / H, title, fontsize=FS_T, fontweight="bold", va="top", ha="left",
             linespacing=1.02)
    if posthoc:
        fig.text(right / W, 1 - (y + 0.9 + (tag_line - 1) * 4.6) / H, "post hoc", fontsize=FS_LG, va="top", ha="right", **POSTHOC_TAG)


REG, SRC = [], {}


def reg(panel, text, value, source, computation, status, lock):
    REG.append(dict(figure="Figure 5", panel=panel, printed_text=text, value=value, source_file=source,
                    computation=computation, matches_locked_text=status, lock_status=lock))


def coh_of(ds):
    return "SCP3342" if ds.startswith("SCP") else "Hahn" if ds.startswith("Hahn") else "MAGNet" if ds.startswith("MAG") else "GSE290367"


CT_SHORT = {"Cardiomyocyte": "CM", "Fibroblast": "Fib", "Endothelial": "EC", "Pericyte": "Peri", "Myeloid": "Mye",
            "VascStromal": "VSC", "bulk": "bulk"}
CT_LONG = {"Cardiomyocyte": "Cardiomyocyte", "Fibroblast": "Fibroblast", "Endothelial": "Endothelial",
           "Pericyte": "Pericyte", "Myeloid": "Myeloid", "VascStromal": "Vasc. stromal"}

# ---------------- data ----------------
S9a = pd.read_excel(XLP, "S9a_impact_audit")
S9b = pd.read_excel(XLP, "S9b_identifiability")
S10a = pd.read_excel(XLP, "S10a_positive_controls")
gw = pd.read_csv(FI + "validation_v1.3/v13_genewise_before_after.csv.gz")
meta = pd.read_csv(FI + "validation_v1.2/steroid_reference/meta_all_genes.csv.gz")
CLS = dict(zip(meta.gene, meta["class"]))
nulls = pd.read_csv(PD + "matched_null_draws.csv.gz")
anull = pd.read_csv(PD + "rev_posthoc_audit_null.csv")
a3 = pd.read_csv(FI + "validation_v1.3/v13_A3_positive_controls_descriptive.csv")
scp = pd.read_csv(PD + "perdonor_scores_SCP3342_HFpEF.csv")
hahn = pd.read_csv(PD + "perdonor_scores_Hahn2021_RVS_HFpEF.csv")
mag = pd.read_csv(PD + "perdonor_scores_MAGNet_exclOverlap_DCM.csv")
gse = pd.read_csv(PD + "perdonor_scores_GSE290367_DCM.csv")

# layout (mm): band 1 y 0-68 (A, B; C spans bands 1-2), band 2 y 68-98 (D), band 3 y 98-166 (E, F), band 4 y 166-240 (G, H, I)
CW = W / 12
LT = dict(fontsize=FS_SM)

# =====================================================================
# A  schematic: identifiable vs aliased
# =====================================================================
xA, yA, wA, hA = 0, 0, 46.0, 68.0
header("A", xA, yA, "Identifiable\nvs aliased", posthoc=True, right=xA + wA - 1.0, tag_line=2)
axA = axmm(xA + 0.6, yA + 11.0, wA - 1.2, hA - 11.0)
axA.set_xlim(0, wA - 1.2); axA.set_ylim(hA - 11.0, 0); axA.axis("off")


def toy_lane(ax, y0, title, ccase, sep, kept):
    ax.text(0.3, y0, title, fontsize=FS_SM, fontweight="bold", va="top")
    xs = np.linspace(0, 19, 200)
    base = y0 + 11.5
    for mu, col in [(9.5 - sep / 2, CONTROL), (9.5 + sep / 2, ccase)]:   # toy curves (illustrative; no data)
        yy = base - 6.6 * np.exp(-0.5 * ((xs - mu) / 2.0) ** 2)
        ax.fill_between(xs + 0.3, yy, base, color=col, alpha=0.45, lw=0)
        ax.plot(xs + 0.3, yy, color=col, lw=0.9)
    ax.plot([0.3, 19.3], [base, base], color=TEXT, lw=0.6)
    ax.text(9.8, base + 0.4, "GC score", ha="center", va="top", **LT)
    ax.add_patch(FancyArrowPatch((20.2, base - 2.6), (26.4, base - 2.6), arrowstyle="-|>", mutation_scale=7, lw=0.8, color=TEXT))
    ax.text(23.3, base - 3.4, "+GC", ha="center", va="bottom", **LT)
    bx1, bx2, bw = 28.4, 37.4, 3.4
    hb, ha_ = 6.8, (6.0 if kept else 0.5)
    ax.add_patch(Rectangle((bx1, base - hb), bw, hb, fc=ccase, ec="none"))
    ax.add_patch(Rectangle((bx2, base - ha_), bw, ha_, fc="white", ec=ccase, lw=0.9))
    ax.plot([27.2, 42.0], [base, base], color=TEXT, lw=0.6)
    ax.text(bx1 + bw / 2, base + 0.4, "base", ha="center", va="top", **LT)
    ax.text(bx2 + bw / 2, base + 0.4, "+GC", ha="center", va="top", **LT)
    ax.text(0.3, base + 4.6, "Disease effect retained" if kept else "Effect removed (VIF ↑)", va="top", style="italic", **LT)


toy_lane(axA, 0.0, "Identifiable", COHORT["MAGNet"], 2.4, True)
toy_lane(axA, 21.6, "Aliased", COHORT["SCP3342"], 11.0, False)
axA.plot([0.3, wA - 2.0], [20.6, 20.6], color="#CCCCCC", lw=0.5)
axA.text(0.3, 43.4, "Non-identifiable if\nAUC ≥0.9 or ≤0.1,\nor VIF ≥5 (illustrative)", va="top", linespacing=1.02, **LT)
reg("A", "AUC ≥0.9 or ≤0.1", 0.9, "Methods (Identifiability)", "threshold, stated", "yes", "T")
reg("A", "VIF ≥5", 5, "Methods (Identifiability)", "threshold, stated", "yes", "T")
SRC["5A"] = pd.DataFrame({"note": ["Schematic, no data. Toy Gaussian curves and bars are illustrative only."],
                          "thresholds": ["AUC >=0.9 or <=0.1, or VIF >=5 (illustrative, proposed post hoc)"]})

# =====================================================================
# B  identifiability map
# =====================================================================
xB, yB, wB, hB = 46.0, 0, 76.0, 68.0
header("B", xB + 1.0, yB, "Identifiability map: HFpEF\ncohorts are aliased", posthoc=True, right=xB + wB - 1.0)
axB = axmm(xB + 15.5, yB + 11.5, 27.0, 42.0)
axB.set_xlim(-0.03, 1.04); axB.set_yscale("log"); axB.set_ylim(0.88, 27)
axB.axvspan(0.9, 1.04, color="#EDEDED", lw=0, zorder=0)
axB.axvspan(-0.03, 0.1, color="#EDEDED", lw=0, zorder=0)
axB.axhspan(5, 27, color="#EDEDED", lw=0, zorder=0)
for v in (0.1, 0.9):
    axB.axvline(v, color="#7F7F7F", lw=0.6, ls="--", zorder=1)
axB.axhline(5, color="#7F7F7F", lw=0.6, ls="--", zorder=1)
axB.set_yticks([1, 2, 5, 10, 20]); axB.set_yticklabels(["1", "2", "5", "10", "20"])
axB.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
axB.set_xticks([0, 0.5, 1.0]); axB.set_xticklabels(["0", "0.5", "1.0"])
axB.set_xlabel("AUC (control > case)", fontsize=FS_AX, labelpad=2)
axB.set_ylabel("Variance inflation,\ndisease coefficient", fontsize=FS_AX, labelpad=2, linespacing=1.0)
srcB = S9b.copy(); srcB["cohort"] = srcB.dataset.map(coh_of)


def coh_marker(c, ax, x, y, ms=5):
    if c == "GSE290367":
        ax.plot(x, y, "^", ms=ms, mfc="white", mec=COHORT[c], mew=1.0, zorder=4, clip_on=False)
    else:
        ax.plot(x, y, "o" if c == "SCP3342" else "s", ms=ms, mfc=COHORT[c], mec=TEXT if c == "Hahn" else COHORT[c],
                mew=0.5, zorder=4, clip_on=False)


for _, r in srcB.iterrows():
    coh_marker(r.cohort, axB, r.AUC_GC3_ctrl_gt_case, r.VIF_case_coef)
scpB = srcB[srcB.cohort == "SCP3342"]
y0b, y1b = scpB.VIF_case_coef.min(), scpB.VIF_case_coef.max()
axB.plot([0.925, 0.925], [y0b, y1b], color=COHORT["SCP3342"], lw=0.9)
for yb_ in (y0b, y1b):
    axB.plot([0.925, 0.95], [yb_, yb_], color=COHORT["SCP3342"], lw=0.9)
axB.text(0.89, (y0b * y1b) ** 0.5, "SCP3342 (6)", ha="right", va="center", color=COHORT["SCP3342"],
         linespacing=1.0, **LT)
colx = xB + 45.0
entries = [("SCP3342", "SCP3342\nAUC 0.974–0.998\nVIF 5.1–12.4", 3), ("Hahn", "Hahn biopsy\nAUC 1.000", 2),
           ("MAGNet", "MAGNet\nAUC 0.80\nVIF 1.33", 3),
           ("GSE290367", "GSE290367\nVIF 1.00–1.09\npooled AUC\n0.20–0.37 (not\ninterpretable)", 5)]
axBk = axmm(colx, yB + 9.6, wB - 45.0, hB - 9.6); axBk.axis("off"); axBk.set_xlim(0, wB - 45.0); axBk.set_ylim(hB - 9.6, 0)
yy = 0.4
for c, txt, nl in entries:
    coh_marker(c, axBk, 1.1, yy + 1.9, ms=4.5)
    axBk.text(2.6, yy, txt, va="top", ha="left", linespacing=1.0, color=TEXT, **LT)
    yy += nl * 4.07 + 1.0
for t, v, comp, st in [("AUC 0.974–0.998", "0.9737;0.9977", "min/max AUC_GC3_ctrl_gt_case, SCP3342 6 compartments", "yes"),
                       ("VIF 5.1–12.4", "5.075;12.364", "min/max VIF_case_coef, SCP3342", "yes"),
                       ("Hahn biopsy AUC 1.000", 1.0, "AUC_GC3_ctrl_gt_case Hahn", "yes"),
                       ("MAGNet AUC 0.80", 0.8018, "AUC MAGNet", "yes"), ("MAGNet VIF 1.33", 1.3349, "VIF MAGNet", "yes"),
                       ("GSE290367 VIF 1.00–1.09", "1.0026;1.0932", "min/max VIF GSE290367 5 cell types", "yes"),
                       ("pooled AUC 0.20–0.37", "0.1969;0.3705", "min/max AUC GSE290367", "yes")]:
    reg("B", t, v, "Supplemental_Tables S9b_identifiability", comp, st, "T")
reg("B", "(point position only) Hahn VIF", 4.424, "S9b", "VIF_case_coef Hahn; not printed", "new", "S")
reg("B", "SCP3342 (6) [6 compartments]", 6, "S9b", "row count", "new", "S")
SRC["5B"] = srcB[["dataset", "cell_type", "cohort", "n", "n_case", "AUC_GC3_ctrl_gt_case", "VIF_case_coef"]]

# =====================================================================
# C  ridgelines of per-donor GC score
# =====================================================================
xC, yC, wC, hC = 122.0, 0, 61.0, 98.0
header("C", xC + 1.0, yC, "Donor GC scores:\nseparate vs overlap", posthoc=True, right=xC + wC - 0.5)
axC = axmm(xC + 24.0, yC + 25.0, wC - 25.5, 61.0)
ridges = [("SCP3342 Fib", scp[scp.cell_type == "Fibroblast"], COHORT["SCP3342"]),
          ("SCP3342 Mye", scp[scp.cell_type == "Myeloid"], COHORT["SCP3342"]),
          ("Hahn biopsy", hahn, COHORT["Hahn"]),
          ("MAGNet", mag, COHORT["MAGNet"]),
          ("GSE Chaffin", gse[(gse.cell_type == "Myeloid") & (gse.study == "Chaffin 2022")], GSE_STUDY["Chaffin"]),
          ("GSE Koenig", gse[(gse.cell_type == "Myeloid") & (gse.study == "Koenig 2022")], GSE_STUDY["Koenig"]),
          ("GSE Reichart", gse[(gse.cell_type == "Myeloid") & (gse.study == "Reichart 2022")], GSE_STUDY["Reichart"])]
XR = (-10.5, 6.5)
xs = np.linspace(*XR, 400)
srcC = []
for i, (lab, d, col) in enumerate(ridges):
    base = -i
    for case, c in [(0, CONTROL), (1, col)]:
        v = d.loc[d.case == case, "score_GC3"].dropna().values
        yy = gaussian_kde(v)(xs); yy = yy / yy.max() * 0.80
        axC.fill_between(xs, base, base + yy, color=c, alpha=0.55, lw=0, zorder=10 - i)
        axC.plot(xs, base + yy, color=c if case else "#7F7F7F", lw=0.8, zorder=10 - i)
        off = -0.07 if case == 0 else -0.19
        axC.vlines(v, base + off - 0.09, base + off, color=c if case else "#7F7F7F", lw=0.5, zorder=11)
        srcC.append(pd.DataFrame({"ridge": lab, "unit_id": d.loc[d.case == case, "unit_id"].values, "case": case,
                                  "score_GC3": d.loc[d.case == case, "score_GC3"].values}))
    n0, n1 = int((d.case == 0).sum()), int((d.case == 1).sum())
    axC.text(XR[0] - 0.4, base + 0.30, lab, ha="right", va="bottom", color=TEXT, clip_on=False, **LT)
    axC.text(XR[0] - 0.4, base + 0.27, f"{n0}/{n1}", ha="right", va="top", color="#555555", clip_on=False, **LT)
axC.text(XR[0] - 0.4, 0.98, "Ctrl/case n", ha="right", va="bottom", color="#555555", style="italic", clip_on=False, **LT)
axC.axvline(0, color="#7F7F7F", lw=0.6, ls="--", zorder=0)
axC.set_xlim(*XR); axC.set_ylim(-len(ridges) + 1 - 0.35, 0.92)
axC.set_yticks([]); axC.spines["left"].set_visible(False)
axC.set_xticks([-10, -5, 0, 5])
axC.set_xlabel("GC score (control SD)", fontsize=FS_AX, labelpad=2)
axCk = axmm(xC + 2.0, yC + 11.0, wC - 3.0, 9.0); axCk.axis("off")
axCk.legend(handles=[Rectangle((0, 0), 1, 1, fc=CONTROL, alpha=0.6, ec="none", label="Control"),
                     Rectangle((0, 0), 1, 1, fc=COHORT["SCP3342"], alpha=0.55, ec="none", label="Case (cohort color)")],
            loc="upper left", ncol=1, fontsize=FS_SM, handlelength=1.0, handletextpad=0.4, borderaxespad=0, labelspacing=0.2)
SRC["5C"] = pd.concat(srcC)
reg("C", "per-ridge n (Ctrl/case)", "24/18;24/15;24/41;102/128;16/11;23/10;6/45", "perdonor_scores_*.csv",
    "counts of donors per group", "new", "N")

# =====================================================================
# D  SCP3342 fibroblast GC score per pool
# =====================================================================
xD, yD, wD, hD = 0, 68.0, 122.0, 30.0
header("D", xD, yD + 0.6, "SCP3342 fibroblasts: separation in every pool", posthoc=True, right=xD + wD - 1.0)
axD = axmm(xD + 15.5, yD + 7.8, 80.0, 13.0)
fib = scp[scp.cell_type == "Fibroblast"].copy()
fib["pool_n"] = fib.pool.str.replace("pool", "").astype(int)
pools = sorted(fib.pool_n.unique())
rng = np.random.default_rng(1)
for j, p in enumerate(pools):
    for case, c, off in [(0, CONTROL, -0.17), (1, COHORT["SCP3342"], 0.17)]:
        v = fib[(fib.pool_n == p) & (fib.case == case)].score_GC3.values
        jit = rng.uniform(-0.07, 0.07, len(v))
        axD.scatter(j + off + jit, v, s=14, c=c, edgecolors="white" if case else "#6F6F6F", linewidths=0.4, zorder=3)
axD.axhline(0, color="#7F7F7F", lw=0.6, ls="--", zorder=0)
axD.set_xticks(range(len(pools))); axD.set_xticklabels([str(p) for p in pools])
axD.set_xlim(-0.5, len(pools) - 0.5); axD.set_ylim(-9.5, 2.5); axD.set_yticks([-8, -4, 0])
axD.set_xlabel("Sequencing pool", fontsize=FS_AX, labelpad=1)
axD.set_ylabel("GC score\n(SD)", fontsize=FS_AX, labelpad=2, linespacing=1.0)
axDk = axmm(xD + 97.5, yD + 9.0, 19.0, 10.0); axDk.axis("off")
axDk.legend(handles=[Line2D([], [], marker="o", ls="none", mfc=CONTROL, mec="#6F6F6F", mew=0.4, ms=4.5, label="Control"),
                     Line2D([], [], marker="o", ls="none", mfc=COHORT["SCP3342"], mec="white", mew=0.4, ms=4.5, label="HFpEF")],
            loc="center left", fontsize=FS_SM, handletextpad=0.2, borderaxespad=0, labelspacing=0.3)
SRC["5D"] = fib[["unit_id", "pool", "case", "disease", "score_GC3"]].sort_values(["pool", "case"])

# =====================================================================
# E  volcano SCP3342 myeloid
# =====================================================================
xE, yE, wE, hE = 0, 98.0, 61.0, 68.0
header("E", xE, yE + 0.5, "SCP3342 Mye: GC-induced\ngenes lead the down list")
fig.text((xE + 5.6) / W, 1 - (yE + 10.3) / H, "Top 100 down: 60% GC-induced\n(3% background); OR 14.2",
         va="top", ha="left", linespacing=1.02, **LT)
axE = axmm(xE + 12.5, yE + 19.6, wE - 14.0, 25.4)
g = gw.query("dataset=='SCP3342_HFpEF' and cell_type=='Myeloid' and model=='base'").copy()
g["cls"] = g.gene.map(CLS).fillna("not_in_reference")
g["nlp"] = -np.log10(g.p)
top = g[g.gene.isin(set(meta.gene)) & (g.beta < 0)].sort_values("p").head(100)
g["top100_down"] = g.gene.isin(top.gene)
order = [("not_in_reference", "#E3E3E3", "#E3E3E3", 3), ("GC_insensitive", GC_CLASS["insensitive"], GC_CLASS["insensitive"], 3),
         ("ambiguous", "white", AMBIGUOUS_EDGE, 3), ("GC_repressed", GC_CLASS["repressed"], GC_CLASS["repressed"], 5),
         ("GC_induced", GC_CLASS["induced"], GC_CLASS["induced"], 6)]
for k, fc, ec, s in order:
    d = g[g.cls == k]
    axE.scatter(d.beta, d.nlp, s=s, c=fc, edgecolors=ec, linewidths=0.3, rasterized=True, zorder=2 if k != "GC_induced" else 3)
tt = g[g.top100_down]
axE.scatter(tt.beta, tt.nlp, s=16, facecolors="none", edgecolors=TEXT, linewidths=0.45, zorder=4)
lab_genes = ["CD163", "DPYD", "FMN1", "SLA", "MERTK", "VSIG4"]
assert all(CLS.get(x) == "GC_induced" and x in set(top.gene) for x in lab_genes)
lab_df = g.set_index("gene").loc[lab_genes].sort_values("nlp", ascending=False)
ly = np.linspace(20.6, 2.2, len(lab_genes))
for (gn, r), yl in zip(lab_df.iterrows(), ly):
    axE.annotate(gn, xy=(r.beta, r.nlp), xytext=(-14.6, yl), style="italic", ha="left", va="center",
                 arrowprops=dict(arrowstyle="-", lw=0.4, color="#888888", shrinkA=0.5, shrinkB=1.5), **LT)
axE.set_xlim(-15, 5.5); axE.set_ylim(-0.3, 21.6)
axE.set_xticks([-10, -5, 0, 5]); axE.set_yticks([0, 10, 20])
axE.axvline(0, color="#7F7F7F", lw=0.6, ls="--", zorder=0)
axE.set_xlabel("HFpEF − Ctrl (SD)", fontsize=FS_AX, labelpad=1)
axE.set_ylabel("−log10 P", fontsize=FS_AX, labelpad=1)
keyE = [Line2D([], [], marker="o", ls="none", mfc=GC_CLASS["induced"], mec=GC_CLASS["induced"], ms=4, label="GC-induced"),
        Line2D([], [], marker="o", ls="none", mfc=GC_CLASS["repressed"], mec=GC_CLASS["repressed"], ms=4, label="GC-repressed"),
        Line2D([], [], marker="o", ls="none", mfc=GC_CLASS["insensitive"], mec=GC_CLASS["insensitive"], ms=4, label="insensitive"),
        Line2D([], [], marker="o", ls="none", mfc="white", mec=AMBIGUOUS_EDGE, ms=4, label="ambiguous"),
        Line2D([], [], marker="o", ls="none", mfc="#E3E3E3", mec="#E3E3E3", ms=4, label="not classified"),
        Line2D([], [], marker="o", ls="none", mfc="none", mec=TEXT, ms=5, mew=0.6, label="top 100 down")]
axEk = axmm(xE + 1.0, yE + 54.6, wE - 1.5, 13.0); axEk.axis("off")
axEk.legend(handles=keyE, loc="upper left", ncol=2, fontsize=FS_SM, handletextpad=0.1, columnspacing=0.6,
            borderaxespad=0, handlelength=0.9, labelspacing=0.15)
rowE = S9a.query("dataset=='SCP3342_HFpEF' and cell_type=='Myeloid' and model=='base'").iloc[0]
assert abs(top.gene.map(CLS).eq("GC_induced").mean() - rowE.top100down_frac_GCinduced) < 1e-9
reg("E", "60% GC-induced", rowE.top100down_frac_GCinduced, "S9a; recomputed from v13_genewise_before_after + meta_all_genes",
    "top-100 rule: in reference, beta<0, sort p, first 100; share GC_induced", "yes", "T")
reg("E", "(3% background)", round(rowE.bg_frac_GCinduced, 4), "S9a", "bg_frac_GCinduced", "yes", "T")
reg("E", "OR 14.2", round(rowE.OR_down_GCinduced, 3), "S9a", "OR_down_GCinduced", "yes", "T")
SRC["5E"] = g[["gene", "beta", "t", "p", "fdr", "cls", "top100_down"]].rename(columns={"cls": "GC_class"})

# =====================================================================
# F  top-100 share vs expression-matched null
# =====================================================================
xF, yF, wF, hF = 61.0, 98.0, 122.0, 68.0
header("F", xF + 1.0, yF + 0.5, "Top-100 GC-induced share vs a matched null", posthoc=True, right=xF + wF - 1.0)
order_F = ([("SCP3342_HFpEF", c) for c in ["Cardiomyocyte", "Fibroblast", "Endothelial", "Pericyte", "Myeloid", "VascStromal"]]
           + [("Hahn2021_RVS_HFpEF", "bulk"), ("MAGNet_exclOverlap_DCM", "bulk")]
           + [("GSE290367_DCM", c) for c in ["Cardiomyocyte", "Fibroblast", "Endothelial", "Pericyte", "Myeloid"]])
xpos, xp, prev = [], 0.0, None
for ds, ct in order_F:
    if prev is not None and coh_of(ds) != prev:
        xp += 0.9
    xpos.append(xp); xp += 1.0; prev = coh_of(ds)
XPOS = {k: x for k, x in zip(order_F, xpos)}
axFu = axmm(xF + 15.5, yF + 20.2, wF - 17.0, 6.0)
axFl = axmm(xF + 15.5, yF + 28.0, wF - 17.0, 19.5)
srcF, n_exceed = [], 0
for (ds, ct), x in zip(order_F, xpos):
    c = COHORT[coh_of(ds)]
    dr = nulls.query("dataset==@ds and cell_type==@ct and model=='base'").top100_frac_GCinduced.values * 100
    vp = axFl.violinplot([dr], positions=[x], widths=0.78, showextrema=False, bw_method=0.35)
    for b in vp["bodies"]:
        b.set_facecolor("#D9D9D9"); b.set_edgecolor("#A0A0A0"); b.set_linewidth(0.4); b.set_alpha(1)
    an = anull.query("dataset==@ds and cell_type==@ct and model=='base'").iloc[0]
    q = an.matched_null_q975 * 100
    sa = S9a.query("dataset==@ds and cell_type==@ct")
    ob = sa[sa.model == "base"].top100down_frac_GCinduced.iloc[0] * 100
    oa = sa[sa.model == "GC3_adjusted"].top100down_frac_GCinduced.iloc[0] * 100
    bg = sa[sa.model == "base"].bg_frac_GCinduced.iloc[0] * 100
    assert abs(ob - an.top100_frac_GCinduced * 100) < 1e-9
    for ax in (axFl, axFu):
        ax.plot([x - 0.36, x + 0.36], [bg, bg], color="#7F7F7F", lw=0.8, ls=(0, (2, 1.5)), zorder=2)
        ax.plot([x - 0.22, x + 0.22], [q, q], color=TEXT, lw=1.1, zorder=3)
        ax.plot(x + 0.14, oa, "o", ms=4.5, mfc="white", mec=c, mew=1.0, zorder=4)
        ax.plot(x - 0.14, ob, "o", ms=4.5, mfc=c, mec=TEXT if coh_of(ds) == "Hahn" else c, mew=0.5, zorder=5)
    exceed = ob > q + 1e-9
    n_exceed += exceed
    if exceed:
        (axFu if ob > 50 else axFl).text(x - 0.14, ob + 0.9 if ob < 50 else ob + 0.6, "†", ha="center", va="bottom", **LT)
    srcF.append(dict(dataset=ds, cell_type=ct, observed_base=ob / 100, observed_GC3_adjusted=oa / 100, bg=bg / 100,
                     matched_null_q975=q / 100, matched_null_mean=an.matched_null_mean, exceeds_q975=bool(exceed),
                     n_null_draws=len(dr)))
assert n_exceed == 9
axFl.set_ylim(0, 27); axFu.set_ylim(56, 64)
for ax in (axFl, axFu):
    ax.set_xlim(xpos[0] - 0.6, xpos[-1] + 0.6)
axFu.spines["bottom"].set_visible(False); axFu.set_xticks([]); axFu.set_yticks([60])
axFl.set_yticks([0, 10, 20])
for ax, yv in [(axFl, 27), (axFu, 56)]:
    ax.plot([xpos[0] - 0.6], [yv], marker=[(-1, -0.5), (1, 0.5)], ms=6, ls="none", color=TEXT, mew=0.6, clip_on=False)
axFl.set_xticks(xpos); axFl.set_xticklabels(["" if ct == "bulk" else CT_SHORT[ct] for _, ct in order_F], rotation=90)
axFl.set_ylabel("GC-induced share,\ntop 100 down (%)", fontsize=FS_AX, labelpad=2, linespacing=1.0)
axFl.yaxis.set_label_coords(-0.075, 0.62)
for cname, lab in [("SCP3342", "SCP3342 (HFpEF)"), ("Hahn", "Hahn"), ("MAGNet", "MAGNet"), ("GSE290367", "GSE290367 (DCM)")]:
    xs_ = [x for (ds, ct), x in zip(order_F, xpos) if coh_of(ds) == cname]
    yb = -13.6
    axFl.plot([min(xs_) - 0.35, max(xs_) + 0.35], [yb, yb], color=COHORT[cname], lw=1.6, clip_on=False, solid_capstyle="butt")
    axFl.text(np.mean(xs_), yb - 1.3, lab, ha="center", va="top", color=TEXT, **LT)
xm = XPOS[("SCP3342_HFpEF", "Myeloid")]
axFu.text(xm + 0.45, 60, "SCP3342 Mye: 60% vs 14% (97.5th pct)", va="center", ha="left", **LT)
xh = XPOS[("Hahn2021_RVS_HFpEF", "bulk")]; xmg = XPOS[("MAGNet_exclOverlap_DCM", "bulk")]
xgm = XPOS[("GSE290367_DCM", "Myeloid")]
axFl.text(xh - 0.32, 21, "21%", ha="right", va="center", **LT)
axFl.text(xmg + 0.32, 18.2, "18% (OR 2.70)", ha="left", va="center", **LT)
axFl.text(xgm + 0.6, 24.0, "10% (OR 1.89)", ha="right", va="center", **LT)
keyF = [Line2D([], [], marker="o", ls="none", mfc="#555555", mec="#555555", ms=4, label="standard model"),
        Line2D([], [], marker="o", ls="none", mfc="white", mec="#555555", ms=4, label="+ GC score"),
        Line2D([], [], color="#7F7F7F", lw=0.8, ls=(0, (2, 1.5)), label="background"),
        Rectangle((0, 0), 1, 1, fc="#D9D9D9", ec="#A0A0A0", lw=0.4, label="matched null (2,000 draws)"),
        Line2D([], [], color=TEXT, lw=1.1, label="97.5th percentile"),
        Line2D([], [], marker="$†$", ls="none", color=TEXT, ms=5, label="above 97.5th pct (9 of 13)")]
axFk = axmm(xF + 15.5, yF + 6.2, wF - 17.0, 13.0); axFk.axis("off")
axFk.legend(handles=keyF, loc="upper left", ncol=2, fontsize=FS_SM, handlelength=1.1, handletextpad=0.3, columnspacing=1.5,
            borderaxespad=0, labelspacing=0.15)
SRC["5F_summary"] = pd.DataFrame(srcF)
SRC["5F_null_draws"] = nulls[nulls.model == "base"]
reg("F", "above 97.5th pct (9 of 13); † marks", 9, "rev_posthoc_audit_null.csv + S9a", "count observed > matched_null_q975, base, 13 contrasts", "yes", "T")
reg("F", "SCP3342 Mye: 60% vs 14% (97.5th pct)", "0.60;0.14", "rev_posthoc_audit_null.csv", "SCP3342 myeloid top100 frac; q975", "yes", "T")
reg("F", "21%", 0.21, "S9a", "Hahn base top100 share", "yes", "T")
reg("F", "18% (OR 2.70)", "0.18;2.703", "S9a", "MAGNet base share; OR_down_GCinduced", "yes", "T")
reg("F", "10% (OR 1.89)", "0.10;1.891", "S9a", "GSE290367 myeloid share; OR", "yes", "T")
reg("F", "matched null (2,000 draws)", 2000, "matched_null_draws.csv.gz", "draws per contrast", "yes", "T")
reg("F", "other contrast shares (plotted, unlabeled) 2-18%", "0.02-0.18", "S9a", "plotted positions", "new", "S")

# =====================================================================
# G  slope chart n FDR<0.05
# =====================================================================
xG, yG, wG, hG = 0, 166.0, 55.0, 74.0
header("G", xG, yG + 0.5, "Significant genes\ncollapse in SCP3342")
axG = axmm(xG + 14.0, yG + 12.5, 13.5, 47.0)
srcG = []
HL = {("MAGNet_exclOverlap_DCM", "bulk"), ("SCP3342_HFpEF", "Fibroblast"), ("SCP3342_HFpEF", "Myeloid")}
for (ds, ct), d in S9a.groupby(["dataset", "cell_type"]):
    b = d[d.model == "base"].n_fdr05.iloc[0]; a = d[d.model == "GC3_adjusted"].n_fdr05.iloc[0]
    c = COHORT[coh_of(ds)]
    hl = (ds, ct) in HL
    axG.plot([0, 1], [b, a], color=c, lw=1.4 if hl else 0.7, alpha=1 if hl else 0.7, zorder=3 if hl else 2)
    axG.plot(0, b, "o", ms=3.5, mfc=c, mec=TEXT if coh_of(ds) == "Hahn" else c, mew=0.4, zorder=4)
    axG.plot(1, a, "o", ms=3.5, mfc="white", mec=c, mew=0.9, zorder=4)
    srcG.append(dict(dataset=ds, cell_type=ct, n_fdr05_base=b, n_fdr05_GC3_adjusted=a))
axG.set_yscale("log"); axG.set_ylim(0.7, 2.2e4); axG.set_xlim(-0.2, 1.2)
axG.set_yticks([1, 10, 100, 1000, 10000]); axG.set_yticklabels(["1", "10", "100", "1k", "10k"])
axG.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
axG.set_xticks([0, 1]); axG.set_xticklabels(["base", "+GC"])
axG.set_ylabel("Genes with FDR <0.05", fontsize=FS_AX, labelpad=2)
for txt, yl, ya, col in [("MAGNet\n9,112 → 6,226", 7000, 6226, COHORT["MAGNet"]),
                         ("SCP3342 Mye\n1,700 → 21", 50, 21, COHORT["SCP3342"]),
                         ("SCP3342 Fib\n8,404 → 5", 3.0, 5, COHORT["SCP3342"])]:
    axG.annotate(txt, xy=(1.0, ya), xytext=(1.45, yl), va="center", ha="left", color=col,
                 linespacing=1.0, annotation_clip=False,
                 arrowprops=dict(arrowstyle="-", lw=0.4, color="#888888", shrinkA=0.5, shrinkB=2.5), **LT)
SRC["5G"] = pd.DataFrame(srcG)
reg("G", "MAGNet 9,112 → 6,226", "9112;6226", "S9a", "n_fdr05 base / GC3_adjusted", "yes", "T")
reg("G", "SCP3342 Fib 8,404 → 5", "8404;5", "S9a", "n_fdr05", "yes", "T")
reg("G", "SCP3342 Mye 1,700 → 21", "1700;21", "S9a", "n_fdr05", "yes", "T")
reg("G", "other lines (unlabeled, e.g. Hahn 8,942 -> 3,835)", "8942;3835", "S9a", "plotted positions", "new", "S")

# =====================================================================
# H  canonical HF genes heatmap
# =====================================================================
xH, yH, wH, hH = 55.0, 166.0, 84.0, 74.0
header("H", xH + 0.5, yH + 0.5, "Canonical HF genes: preserved where\nidentifiable, lost in SCP3342")
genesH = ["NPPA", "NPPB", "MYH6", "POSTN", "ASPN", "SFRP4", "FRZB"]
EXPECTED = {"NPPA": 1, "NPPB": 1, "POSTN": 1, "ASPN": 1, "SFRP4": 1, "FRZB": 1, "MYH6": -1}  # planned directions (adjustment positive control)
colsH = [("MAGNet_exclOverlap_DCM", "bulk", "MAGNet\nbulk"), ("GSE290367_DCM", "Cardiomyocyte", "GSE\nCM"),
         ("GSE290367_DCM", "Fibroblast", "GSE\nFib"), ("SCP3342_HFpEF", "Fibroblast", "SCP\nFib"),
         ("SCP3342_HFpEF", "Cardiomyocyte", "SCP\nCM"), ("Hahn2021_RVS_HFpEF", "bulk", "Hahn\nbiopsy")]
mag_a3 = S10a[S10a.dataset == "MAGNet_exclOverlap_DCM"][["gene", "beta", "fdr", "dataset", "cell_type", "model"]]
hm_src = pd.concat([mag_a3, a3[["gene", "beta", "fdr", "dataset", "cell_type", "model"]]])
PAIR = 2.5
NXU = 6 * PAIR - 0.5
HW = 57.0
RHc = 4.0   # v0.6: 4.3 mm
X0H = xH + 12.5
TOPH = yH + 27.6   # v0.6: +25.6
MMU = lambda mm_: mm_ / RHc   # mm -> heatmap row units (keeps v0.6 mm offsets)
axH = axmm(X0H, TOPH, HW, 7 * RHc)
cmapH = diverging_cmap()
normH = mpl.colors.TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
srcH, counts = [], []
for j, (ds, ct, lab) in enumerate(colsH):
    x0p = j * PAIR
    for k, model in enumerate(["base", "GC3_adjusted"]):
        x = x0p + k
        nsig, nmeas, nopp = 0, 0, 0
        for i, gn in enumerate(genesH):
            r = hm_src[(hm_src.dataset == ds) & (hm_src.cell_type == ct) & (hm_src.model == model) & (hm_src.gene == gn)]
            if len(r) == 0:
                axH.add_patch(Rectangle((x, i), 1, 1, fc=GC_CLASS["not_measured"], ec="#BDBDBD", lw=0.0, hatch="//////"))
                srcH.append(dict(dataset=ds, cell_type=ct, model=model, gene=gn, beta=np.nan, fdr=np.nan, measured=False,
                                 expected_sign=EXPECTED[gn], call="n.d."))
                continue
            r = r.iloc[0]; nmeas += 1
            axH.add_patch(Rectangle((x, i), 1, 1, fc=cmapH(normH(np.clip(r.beta, -3, 3))), ec="white", lw=0.5))
            call = "ns"
            gcol = "white" if abs(r.beta) > 1.9 else TEXT
            if r.fdr < 0.05 and np.sign(r.beta) == EXPECTED[gn]:
                nsig += 1; call = "sig_expected"
                axH.plot(x + 0.5, i + 0.5, "o", ms=2.8, mfc=gcol, mec="none")
            elif r.fdr < 0.05:
                nopp += 1; call = "sig_opposite"
                axH.plot(x + 0.5, i + 0.5, "x", ms=3.4, mew=0.9, color=gcol)
            srcH.append(dict(dataset=ds, cell_type=ct, model=model, gene=gn, beta=r.beta, fdr=r.fdr, measured=True,
                             expected_sign=EXPECTED[gn], call=call))
        counts.append((j, k, nsig, nmeas, nopp))
        axH.text(x + 0.5, 7 + MMU(0.645), "base" if k == 0 else "+GC", ha="center", va="top", rotation=90, **LT)
        axH.text(x + 0.5, 7 + MMU(9.03), f"{nsig}/{nmeas}", ha="center", va="top", rotation=90, **LT)
    axH.text(x0p + 1 - (0.5 if j == 0 else 0), -MMU(1.634), lab, ha="center", va="bottom", linespacing=0.95, color=TEXT, **LT)
    axH.plot([x0p + 0.06, x0p + 1.94], [-MMU(0.688), -MMU(0.688)], color=COHORT[coh_of(ds)], lw=2.0, solid_capstyle="butt", clip_on=False)
axH.set_xlim(0, NXU); axH.set_ylim(7, 0)
axH.set_yticks(np.arange(7) + 0.5); axH.set_yticklabels(genesH, style="italic", fontsize=FS_SM)
axH.set_xticks([]); axH.tick_params(length=0, pad=1.5)
for s in axH.spines.values():
    s.set_visible(False)
xb_mm = X0H + (1.0 / NXU) * HW
BADGE_H = "Planned check: 7 of 7 canonical\nHF genes retained in MAGNet"   # v0.6: "A3 passed, 7/7"
fig.text((xH + 0.9) / W, 1 - (yH + 10.6) / H, BADGE_H, fontsize=FS_SM, color=TEXT, va="top", ha="left",
         linespacing=1.0, bbox=dict(boxstyle="round,pad=0.25,rounding_size=0.8", fc="#EDEDED", ec=TEXT, lw=0.6))
axHk = axmm(xH + 54.6, yH + 9.1, 32.0, 8.6); axHk.axis("off")
axHk.legend(handles=[Line2D([], [], marker="o", ls="none", mfc=TEXT, mec="none", ms=3, label="FDR<0.05"),
                     Line2D([], [], marker="x", ls="none", color=TEXT, mew=0.9, ms=3.6, label="opposite"),
                     Rectangle((0, 0), 1, 1, fc=GC_CLASS["not_measured"], ec="#BDBDBD", lw=0.4, hatch="//////", label="n.d.")],
            loc="upper left", ncol=2, fontsize=FS_SM, handlelength=0.8, handletextpad=0.25, borderaxespad=0, columnspacing=0.6,
            labelspacing=0.2)
caxH = axmm(X0H + HW + 1.4, TOPH, 2.0, 7 * RHc)
cb = mpl.colorbar.ColorbarBase(caxH, cmap=cmapH, norm=normH, extend="both", ticks=[-3, 0, 3])
cb.outline.set_linewidth(0.4); caxH.tick_params(labelsize=FS_SM, length=2, pad=1)
cb.set_label("Case − Ctrl (SD)", fontsize=FS_SM, labelpad=1.5)   # v0.6: labelpad=0
SRC["5H"] = pd.DataFrame(srcH)
cnt = {(colsH[j][2].replace("\n", " "), "base" if k == 0 else "+GC"): (s, m, o) for j, k, s, m, o in counts}
assert cnt[("MAGNet bulk", "base")][:2] == (7, 7) and cnt[("MAGNet bulk", "+GC")][:2] == (7, 7)
assert cnt[("SCP Fib", "base")][:2] == (4, 6) and cnt[("SCP Fib", "+GC")][:2] == (0, 6)
assert cnt[("SCP CM", "base")][:2] == (3, 5) and cnt[("SCP CM", "+GC")][:2] == (0, 5)
assert cnt[("Hahn biopsy", "base")][:2] == (5, 7) and cnt[("Hahn biopsy", "+GC")][:2] == (3, 7)
assert cnt[("GSE CM", "base")] == (3, 7, 1) and cnt[("GSE CM", "+GC")] == (3, 7, 1)
assert cnt[("GSE Fib", "base")][:2] == (7, 7) and cnt[("GSE Fib", "+GC")][:2] == (7, 7)
okc = a3.groupby(["dataset", "cell_type", "model"]).ok.sum().to_dict()
assert okc[("GSE290367_DCM", "Cardiomyocyte", "base")] == 3 and okc[("SCP3342_HFpEF", "Fibroblast", "base")] == 4
reg("H", "Planned check: 7 of 7 canonical HF genes retained in MAGNet; MAGNet 7/7 -> 7/7", "7;7", "S10a_positive_controls",
    "count expected-direction FDR<0.05 of 7 measured, MAGNet base and +GC (adjustment positive control met)", "yes", "T")
reg("H", "SCP Fib 4/6 -> 0/6", "4;0", "v13_A3_positive_controls_descriptive.csv", "count expected-direction FDR<0.05 of measured", "yes", "T")
reg("H", "Hahn NPPA lost (dot absent after adjustment)", "fdr 0.029 -> 0.977", "v13_A3_positive_controls_descriptive.csv", "NPPA fdr base/adj", "yes", "T")
reg("H", "GSE CM 3/7 -> 3/7 (NPPA, NPPB, MYH6); ASPN x = FDR<0.05 opposite direction", "3;3", "v13_A3_positive_controls_descriptive.csv (ok flag)", "count expected-direction FDR<0.05", "yes", "T")
reg("H", "GSE Fib 7/7 -> 7/7", "7;7", "v13_A3_positive_controls_descriptive.csv", "count", "new", "N")
reg("H", "SCP CM 3/5 -> 0/5", "3;0", "v13_A3_positive_controls_descriptive.csv", "count", "new", "S")
reg("H", "Hahn biopsy 5/7 -> 3/7", "5;3", "v13_A3_positive_controls_descriptive.csv", "count", "new", "S")

# =====================================================================
# I  gene-wise t before vs after
# =====================================================================
xI, yI, wI, hI = 139.0, 166.0, 44.0, 74.0
header("I", xI + 0.5, yI + 0.5, "Gene-wise t:\npreserved vs\ncollapsed")
srcI = []
for k, (ds, ct, lab, col) in enumerate([("MAGNet_exclOverlap_DCM", "bulk", "MAGNet", COHORT["MAGNet"]),
                                        ("SCP3342_HFpEF", "Fibroblast", "SCP3342 Fib", COHORT["SCP3342"])]):
    ax = axmm(xI + 13.0, yI + 17.5 + k * 25.5, wI - 14.0, 20.0)
    d = gw.query("dataset==@ds and cell_type==@ct").pivot(index="gene", columns="model", values="t").dropna()
    d["cls"] = d.index.map(CLS)
    ax.hexbin(d.base, d.GC3_adjusted, gridsize=34, bins="log", cmap="Greys", mincnt=1, linewidths=0, extent=(-24, 24, -24, 24),
              rasterized=True, vmax=2000)
    gi = d[d.cls == "GC_induced"]
    ax.scatter(gi.base, gi.GC3_adjusted, s=2.5, c=GC_CLASS["induced"], lw=0, zorder=3, rasterized=True)
    ax.plot([-24, 24], [-24, 24], color="#7F7F7F", lw=0.6, ls="--", zorder=4)
    ax.axhline(0, color="#BBBBBB", lw=0.4, zorder=1)
    ax.set_xlim(-24, 24); ax.set_ylim(-24, 24); ax.set_xticks([-20, 0, 20]); ax.set_yticks([-20, 0, 20])
    ax.plot(-21.0, 19.5, "s", ms=4, mfc=col, mec="none", clip_on=False)
    ax.text(-18.0, 19.5, lab, va="center", ha="left", color=TEXT, fontweight="bold", **LT)
    if k == 0:
        ax.text(23, -23, "residual\nOR 2.07", ha="right", va="bottom", linespacing=1.0, **LT)
        ax.set_xticklabels([])
    else:
        ax.set_xlabel("t, standard model", fontsize=FS_AX, labelpad=1)
    ax.set_ylabel("t, + GC", fontsize=FS_AX, labelpad=0)
    d2 = d.reset_index(); d2.insert(0, "contrast", f"{ds}|{ct}"); srcI.append(d2)
SRC["5I"] = pd.concat(srcI)
orI = S9a.query("dataset=='MAGNet_exclOverlap_DCM' and model=='GC3_adjusted'").OR_down_GCinduced.iloc[0]
assert f"{orI:.2f}" == "2.07"
reg("I", "residual OR 2.07", round(orI, 4), "S9a", "OR_down_GCinduced MAGNet GC3_adjusted", "yes", "T")

for p, t in [("B", "x 0, 0.5, 1.0; y 1, 2, 5, 10, 20"), ("C", "x -10, -5, 0, 5"), ("D", "y -8, -4, 0; pools 2-10 (identifiers)"),
             ("E", "x -10, -5, 0, 5; y 0, 10, 20"), ("F", "y 0, 10, 20, 60"), ("G", "y 1, 10, 100, 1k, 10k"),
             ("H", "color bar -3, 0, 3"), ("I", "x/y -20, 0, 20")]:
    reg(p, "axis tick labels: " + t, "", "axis", "axis scale only", "new", "N")
for t, v, src_, st, lk in [("13 contrasts", 13, "S9b rows", "yes", "T"), ("Hahn biopsy 24/41 (Ctrl/case)", "24;41", "S9b n/n_case", "yes", "T"),
                          ("MAGNet 102/128 (Ctrl/case)", "102;128", "S9b n/n_case", "yes", "T"),
                          ("SCP3342 6 compartments; GSE290367 5 cell types", "6;5", "S9b", "new", "S"),
                          ("2,000 draws", 2000, "matched_null_draws", "yes", "T"), ("top 100", 100, "Methods", "yes", "T"),
                          ("3 score genes", 3, "Methods", "yes", "T"), ("8 pools", 8, "perdonor_scores_SCP3342 pool / data_inventory", "new", "S"),
                          ("thresholds AUC >=0.9 or <=0.1, VIF >=5 (B gray zone)", "0.9;0.1;5", "Methods", "yes", "T")]:
    reg("legend", t, v, src_, "legend text", st, lk)


# ---------------- QA ----------------
def qa(fig):
    fig.canvas.draw(); r = fig.canvas.get_renderer(); fb = fig.bbox
    T = [t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text().strip()]
    bb = [(t, t.get_window_extent(r)) for t in T]
    out = [t.get_text()[:30] for t, b in bb if b.x0 < fb.x0 - 1 or b.x1 > fb.x1 + 1 or b.y0 < fb.y0 - 1 or b.y1 > fb.y1 + 1]
    ov = []
    for i, (a, ba) in enumerate(bb):
        if isinstance(a, mpl.text.Annotation):
            continue
        for b, bbb in bb[i + 1:]:
            if isinstance(b, mpl.text.Annotation):
                continue
            ix = min(ba.x1, bbb.x1) - max(ba.x0, bbb.x0); iy = min(ba.y1, bbb.y1) - max(ba.y0, bbb.y0)
            if ix > 1 and iy > 1:
                ov.append((a.get_text()[:25], b.get_text()[:25]))
    return out, ov


BANNED = re.compile(r"\b(G[1-5]|G4[ab]|C[1-3]|K[12]|A[1-4])\b|\bv1\.[0-3]\b|version 1\.[0-3]|SHA-?256|\b(frozen|freez\w*|lock\w*)\b"
                    r"|pipeline version|\b(passed|failed)\b", re.I)


def fig_texts(fig):
    return [t.get_text() for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text().strip()]


def number_trace(fig, reg_rows):
    """Every numeric token printed in the figure must occur in a registry printed_text entry."""
    norm = lambda z: z.replace("\u2212", "-").replace("\u2013", "-").replace(",", "")
    pool = norm(" | ".join(f"{r['printed_text']} | {r['value']}" for r in reg_rows))
    pool_tok = set(re.findall(r"\d+(?:\.\d+)?", pool))
    miss = []
    for t in fig_texts(fig):
        for tok in re.findall(r"(?<![A-Za-z\d.])\d+(?:\.\d+)?", norm(t)):   # skip gene symbols (CD163)
            if tok not in pool_tok:
                miss.append((t[:40], tok))
    return miss


_o, _v = qa(fig)
_m = min_font_audit(fig, 10)
_b = [t for t in fig_texts(fig) if BANNED.search(t)]
_n = number_trace(fig, REG)
print("outside:", _o); print("overlaps:", _v); print("min_font_audit:", _m)
print("banned terms in figure:", _b); print("untraced printed numbers:", _n)
assert _m == [] and _o == [] and _v == [] and _b == [] and _n == [], (_m, _o, _v, _b, _n)
save_figure(fig, "Figure_5")
pd.DataFrame(REG).to_csv("Figure_5_number_registry.csv", index=False)
with pd.ExcelWriter("Figure_5_source_data.xlsx") as xw:
    for k, v in SRC.items():
        v.to_excel(xw, sheet_name=k[:31], index=False)
print("saved Figure_5", FONT)
