"""Figure 6 (R5). Genes interpreted as HF or LVAD-recovery biology are glucocorticoid targets in vitro.
Standalone: rebuilds Figure_6.{pdf,png,tif}, Figure_6_source_data.xlsx, Figure_6_number_registry.csv
from the project input artifacts. Display-level derivations only (ranks, counts); no new tests.
"""
import os, tarfile
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch, Polygon
from matplotlib.lines import Line2D

import os  # public package: inputs resolved relative to this script (override with JBTS_INPUTS)
HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.environ.get("JBTS_INPUTS", os.path.join(HERE, "..", "figures_v05", "inputs"))
STYLE = os.path.join(HERE, "jbts_style.py")  # v0.6: >=10 pt lettering
TAR_FI = os.path.join(INP, "figure_inputs_v05.tar.gz")
TAR_PD = os.path.join(INP, "figure_inputs_v05_perdonor.tar.gz")
XLP = os.path.join(INP, "Supplemental_Tables_S1-S14.xlsx")
# Per-donor z_ref for PDK4/TXNIP/FOXO3/SGK1: frozen v1.3 dump re-run on the server with DG extended by these 4 genes
# (server: 04_results/figure_inputs_v05/v13_rerun_lvad4/; all frozen outputs reproduced with max |diff| = 0, see verify_lvad4.csv)
GE_LVAD4 = os.path.join(INP, "perdonor_gene_expression_GC8_clock10_lvad4.csv.gz")
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
mpl.rcParams["hatch.linewidth"] = 0.5
FS_T, FS_AX, FS_TK, FS_LG, FS_SM = SIZES["panel_title"], SIZES["axis_label"], SIZES["tick"], SIZES["legend"], SIZES["min_gene_label"]

W, H = 183.0, 240.0
fig = plt.figure(figsize=(W * MM, H * MM))
fig.patches.append(Rectangle((0, 0), 1, 1, transform=fig.transFigure, fc="white", ec="none", zorder=-10))
CW, RH = W / 12, H / 7


def axmm(x, y, w, h, **kw):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **kw)


def header(letter, x, y, title, right=None, tag=None):
    fig.text(x / W, 1 - y / H, letter, fontsize=SIZES["panel_letter"], fontweight="bold", va="top", ha="left")
    fig.text((x + 5.6) / W, 1 - (y + 0.7) / H, title, fontsize=FS_T, fontweight="bold", va="top", ha="left", linespacing=1.02)
    if tag:
        fig.text(right / W, 1 - (y + 0.9) / H, tag, fontsize=FS_LG, va="top", ha="right", **POSTHOC_TAG)


REG, SRC = [], {}


def reg(panel, text, value, source, computation, status, lock):
    REG.append(dict(figure="Figure 6", panel=panel, printed_text=text, value=value, source_file=source,
                    computation=computation, matches_locked_text=status, lock_status=lock))


CLS_COL = {"GC_induced": GC_CLASS["induced"], "GC_repressed": GC_CLASS["repressed"], "ambiguous": "#FFFFFF",
           "GC_insensitive": GC_CLASS["insensitive"], "not measured in >=2": GC_CLASS["not_measured"]}
CLS_LAB = {"GC_induced": "GC-induced", "GC_repressed": "GC-repressed", "ambiguous": "Ambiguous",
           "not measured in >=2": "Not measured"}
POS_CTRL = {"FKBP5", "TSC22D3"}


def cls_marker(ax, x, y, cls, gene, ms=4.8, z=5):
    fc = CLS_COL[cls]; ec = AMBIGUOUS_EDGE if cls == "ambiguous" else (TEXT if gene in POS_CTRL else fc)
    if gene in POS_CTRL:  # positive control: thick black ring (hatching is reserved for missing data)
        ax.plot(x, y, "o", ms=ms * 1.1, mfc=fc, mec=TEXT, mew=1.4, zorder=z)
    else:
        ax.plot(x, y, "o", ms=ms, mfc=fc, mec=ec, mew=0.7 if cls == "ambiguous" else 0.4, zorder=z)


# ---------------- data ----------------
S10b = pd.read_excel(XLP, "S10b_literature_genes")
A4 = pd.read_csv(FI + "validation_v1.3/A4_literature_gene_GC_classification.csv")
assert (S10b[["source", "gene", "GC_class"]].values == A4[["source", "gene", "GC_class"]].values).all()
S5 = pd.read_excel(XLP, "S5_GC_classification")
meta = pd.read_csv(FI + "validation_v1.2/steroid_reference/meta_all_genes.csv.gz")
assert len(S5) == len(meta) == 13368
pdl = pd.read_csv(FI + "validation_v1.2/steroid_reference/per_dataset_log2fc.csv.gz")
S3 = pd.read_excel(XLP, "S3_gene_sets")
ge = pd.read_csv(GE_LVAD4)
_ge_old = pd.read_csv(PD + "perdonor_gene_expression_GC8_clock10.csv.gz")
_chk = ge[ge.gene.isin(_ge_old.gene.unique())].reset_index(drop=True)
assert np.allclose(_chk.z_ref.values, _ge_old.z_ref.values, equal_nan=True, rtol=0, atol=0)
scp = pd.read_csv(PD + "perdonor_scores_SCP3342_HFpEF.csv")
hahn = pd.read_csv(PD + "perdonor_scores_Hahn2021_RVS_HFpEF.csv")
mag = pd.read_csv(PD + "perdonor_scores_MAGNet_exclOverlap_DCM.csv")
gse = pd.read_csv(PD + "perdonor_scores_GSE290367_DCM.csv")


def src_key(s):
    return ("CircRes" if s.startswith("Circ Res") else "NatCommun" if s.startswith("Nat Commun")
            else "LVAD" if s.startswith("JACC") else "Amrute")


S10b["src"] = S10b.source.map(src_key)
GCLS = dict(zip(S10b.gene, S10b.GC_class))

# =====================================================================
# A  alluvial: source -> gene -> in vitro class
# =====================================================================
xA, yA, wA, hA = 0, 0, W, 62.0
header("A", xA, yA, "Genes read as HF or recovery biology map onto the in vitro GC classes")
axA = axmm(0, 7.0, W, hA - 7.0); axA.set_xlim(0, W); axA.set_ylim(hA - 7.0, 0); axA.axis("off")
SOURCES = [("CircRes", "HFpEF snRNA-seq (ref 3)", "FKBP5 down; PGR, MYOF\nup in CM", 46.0),
           ("NatCommun", "HF map (ref 2)", "FKBP5: stable\nEC marker", 29.0),
           ("LVAD", "LVAD (ref 22)", "'GC receptor signature'\nup: 8 of 9 GC-induced", 46.0),
           ("Amrute", "Recovery (ref 23)", "up in all 5 cell types:\n3 of 5 GC-induced", 43.0)]
GENES = ["PGR", "MYOF", "SGK1", "FKBP5", "TSC22D3", "ZBTB16", "PDK4", "DDIT4", "KLF9", "TXNIP", "FOXO3", "FGD4",
         "PTPN13", "ZMYND8"]
CLASSES = ["not measured in >=2", "GC_repressed", "GC_induced", "ambiguous"]
assert set(GENES) == set(S10b.gene)
links = S10b[["src", "gene", "GC_class"]].values.tolist()
UW = 2.3  # ribbon width (mm) per link
gap = (W - sum(s[3] for s in SOURCES) - 2.0) / 3
sx, x = {}, 1.0
for k, t1, t2, w in SOURCES:
    sx[k] = (x, w); x += w + gap
Y_SRC0, Y_SRC1 = 0.4, 13.8
Y_G0, Y_G1 = 25.5, 27.1
Y_C0, Y_C1 = 44.6, 54.2
slot = (W - 4.0) / len(GENES)
gx = {g: 2.0 + slot * (i + 0.5) for i, g in enumerate(GENES)}
gin = {g: sum(1 for l in links if l[1] == g) for g in GENES}
for k, t1, t2, w in SOURCES:
    x0, w0 = sx[k]
    axA.add_patch(FancyBboxPatch((x0, Y_SRC0), w0, Y_SRC1 - Y_SRC0, boxstyle="round,pad=0,rounding_size=1.2",
                                 fc="#F4F4F4", ec="#9A9A9A", lw=0.6))
    axA.text(x0 + 1.2, Y_SRC0 + 0.7, t1, fontsize=FS_SM, fontweight="bold", va="top", ha="left")
    axA.text(x0 + 1.2, Y_SRC0 + 5.0, t2, fontsize=FS_SM, va="top", ha="left", linespacing=1.0, style="italic")


def ribbon(ax, xt, yt, xb, yb, w, fc, ec, hatch=None, alpha=0.6):
    t = np.linspace(0, 1, 40); s = t * t * (3 - 2 * t)
    y = yt + (yb - yt) * t
    xl = xt + (xb - xt) * s
    P = np.vstack([np.c_[xl, y], np.c_[xl[::-1] + w, y[::-1]]])
    ax.add_patch(Polygon(P, closed=True, fc=fc, ec=ec, lw=0.3, alpha=alpha, hatch=hatch, zorder=2))


# source slots ordered by target gene x; gene slots ordered by source x
src_order = [s[0] for s in SOURCES]
ribs = []
gene_slots_used = {g: 0 for g in GENES}
for k in src_order:
    L = sorted([l for l in links if l[0] == k], key=lambda l: gx[l[1]])
    x0, w0 = sx[k]; tot = len(L) * UW
    start = x0 + (w0 - tot) / 2
    for i, l in enumerate(L):
        ribs.append([k, l[1], l[2], start + i * UW])
for g in GENES:
    R = sorted([r for r in ribs if r[1] == g], key=lambda r: src_order.index(r[0]))
    base = gx[g] - gin[g] * UW / 2
    for i, r in enumerate(R):
        r.append(base + i * UW)
for k, g, c, xt, xb in ribs:
    fc = CLS_COL[c]; ec = AMBIGUOUS_EDGE if c in ("ambiguous", "not measured in >=2") else "none"
    ribbon(axA, xt, Y_SRC1, xb, Y_G0, UW, fc if c != "ambiguous" else "#F2F2F2", ec,
           hatch="//////" if c == "not measured in >=2" else None, alpha=0.75 if c in ("GC_induced", "GC_repressed") else 1)
# gene nodes + labels
for g in GENES:
    c = GCLS[g]; wN = gin[g] * UW
    axA.add_patch(Rectangle((gx[g] - wN / 2, Y_G0), wN, Y_G1 - Y_G0, fc=CLS_COL[c] if c != "ambiguous" else "white",
                            ec=TEXT if c in ("ambiguous", "not measured in >=2") else CLS_COL[c], lw=0.5, zorder=3,
                            hatch="//////" if c == "not measured in >=2" else None))
    axA.text(gx[g], Y_G1 + 0.5 + 4.5 * (GENES.index(g) % 2), g, fontsize=FS_SM, style="italic", ha="center", va="top", zorder=4,
             bbox=dict(fc="white", ec="none", pad=0.3))  # alternate rows so 10-pt labels do not collide
# gene -> class ribbons and class nodes
cls_genes = {c: [g for g in GENES if GCLS[g] == c] for c in CLASSES}
for c in CLASSES:
    gs = cls_genes[c]
    x0 = min(gx[g] - gin[g] * UW / 2 for g in gs) - 1.2
    x1 = max(gx[g] + gin[g] * UW / 2 for g in gs) + 1.2
    if c == "not measured in >=2":
        x0, x1 = 0.3, 16.7
    if c == "GC_repressed":
        x0, x1 = 17.3, 40.6
    if c == "ambiguous":
        x0, x1 = x0 - 1.0, min(x1 + 1.0, W - 0.3)
    for g in gs:
        wN = gin[g] * UW
        ribbon(axA, gx[g] - wN / 2, Y_G1 + 1.0, gx[g] - wN / 2, Y_C0, wN, CLS_COL[c] if c != "ambiguous" else "#F2F2F2",
               AMBIGUOUS_EDGE if c in ("ambiguous", "not measured in >=2") else "none",
               hatch="//////" if c == "not measured in >=2" else None, alpha=0.45 if c in ("GC_induced", "GC_repressed") else 1)
    fc = CLS_COL[c] if c in ("GC_induced", "GC_repressed") else ("white" if c == "ambiguous" else GC_CLASS["not_measured"])
    axA.add_patch(FancyBboxPatch((x0, Y_C0), x1 - x0, Y_C1 - Y_C0, boxstyle="round,pad=0,rounding_size=1.0", fc=fc,
                                 ec=TEXT if c in ("ambiguous", "not measured in >=2") else fc, lw=0.6, zorder=3))
    n = len(gs)
    txt = f"{CLS_LAB[c]}: {n} gene{'s' if n > 1 else ''}"
    if c == "not measured in >=2":
        txt = "Not\nmeasured"
    elif c == "GC_repressed":
        txt = "GC-repressed\n(2 genes)"
    elif c == "ambiguous":
        txt = "Ambiguous\n(2 genes)"
    tc = "white" if c in ("GC_induced", "GC_repressed") else TEXT
    axA.text((x0 + x1) / 2, (Y_C0 + Y_C1) / 2, txt, fontsize=FS_SM, ha="center", va="center", color=tc,
             fontweight="bold" if c == "GC_induced" else "normal", zorder=4, linespacing=1.0)
SRC["6A_links"] = S10b[["source", "src", "gene", "k", "meta_log2FC", "ci", "GC_class"]]
reg("A", "8 of 9 GC-induced (LVAD 'GC receptor signature')", "8/9", "S10b_literature_genes", "count GC_induced, JACC BTS 2026 rows", "yes", "T")
reg("A", "3 of 5 GC-induced (recovery)", "3/5", "S10b_literature_genes", "count GC_induced, Amrute rows", "yes", "T")
reg("A", "FKBP5 -> GC-induced", "GC_induced", "S10b", "class", "yes", "T")
reg("A", "SGK1 -> GC-repressed (exception)", "GC_repressed", "S10b", "class", "yes", "T")
reg("A", "PGR not measured in >=2 datasets", "not measured", "S10b", "class", "yes", "T")
reg("A", "MYOF -> GC-repressed", "GC_repressed", "S10b", "class", "new", "S")
reg("A", "GC-induced: 9 genes; GC-repressed (2 genes); Ambiguous (2 genes); Not measured", "9;2;2;1", "S10b", "unique genes per class (14 unique of 18 links)", "new", "N")
reg("A", "FKBP5 (4), ZBTB16 (2) source counts", "4;2", "S10b", "number of sources reporting the gene", "new", "N")
reg("A", "'5 genes up in all 5 cell types'", "5;5", "S10b", "paraphrase of cited study (gene count, cell types)", "new", "N")
reg("A", "ref 3, ref 2, ref 22, ref 23 (citation numbers)", "3;2;22;23", "lead instruction: final typeset numbering", "citation identifiers, not data", "new", "N")

# =====================================================================
# B  forest of meta log2FC
# =====================================================================
xB, yB, wB, hB = 0, 62.0, 76.0, 100.0
header("B", xB, yB + 1.0, "In vitro GC effect")
mt = meta.set_index("gene")
groups = [("LVAD signature (9)", ["ZBTB16", "PDK4", "FKBP5", "DDIT4", "KLF9", "TSC22D3", "TXNIP", "FOXO3", "SGK1"]),
          ("Recovery (5; + FKBP5, ZBTB16)", ["FGD4", "ZMYND8", "PTPN13"]),
          ("HFpEF / HF map (+ FKBP5)", ["MYOF", "PGR"])]
rows = []  # (y, kind, label)
yy = 0
for gname, gl in groups:
    rows.append((yy, "hdr", gname)); yy += 1
    for g in gl:
        rows.append((yy, "gene", g)); yy += 1
NROW = yy
ROW_H = 4.2
TOPY = yB + 17.0
axB = axmm(xB + 20.5, TOPY, 43.0, NROW * ROW_H)
axB.set_ylim(NROW - 0.5, -0.5); axB.set_xlim(-2.2, 8.2)
axB.axvline(0, color="#7F7F7F", lw=0.6, ls="--", zorder=0)
for v in (-0.5, 0.5):
    axB.axvline(v, color="#BDBDBD", lw=0.6, ls=":", zorder=0)
srcB = []
for y, kind, lab in rows:
    if kind == "hdr":
        axB.text(-0.31, y, lab, fontsize=FS_SM, fontweight="bold", va="center", ha="left", color="#444444",
                 transform=axB.get_yaxis_transform(), clip_on=False, zorder=6,
                 bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none"))
        continue
    g = lab; c = GCLS[g]
    nsrc = int((S10b.gene == g).sum())
    glab = g + (f" ({nsrc})" if nsrc > 1 else "")
    axB.text(-0.03, y, glab, fontsize=FS_SM, style="italic", va="center", ha="right", transform=axB.get_yaxis_transform())
    if g not in mt.index:
        axB.text(0.3, y, "not measured", fontsize=FS_SM, va="center", ha="left", color="#555555")
        srcB.append(dict(gene=g, GC_class=c, meta_log2FC=np.nan, ci_lo=np.nan, ci_hi=np.nan, k=np.nan, n_sources=nsrc))
        continue
    r = mt.loc[g]
    col = CLS_COL[c] if c != "ambiguous" else AMBIGUOUS_EDGE
    axB.plot([r.ci_lo, r.ci_hi], [y, y], color=col, lw=0.9, zorder=3, solid_capstyle="butt")
    cls_marker(axB, r.meta_log2FC, y, c, g)
    s10 = S10b[S10b.gene == g].iloc[0]
    assert f"{r.meta_log2FC:.2f}" == f"{s10.meta_log2FC:.2f}"
    axB.text(1.02, y, f"{r.meta_log2FC:+.2f}".replace("-", "−"), fontsize=FS_SM, va="center", ha="left",
             transform=axB.get_yaxis_transform())
    srcB.append(dict(gene=g, GC_class=c, meta_log2FC=r.meta_log2FC, ci_lo=r.ci_lo, ci_hi=r.ci_hi, k=int(r.k), n_sources=nsrc))
axB.text(1.02, -1.25, "log2FC", fontsize=FS_SM, va="center", ha="left", transform=axB.get_yaxis_transform(),
         color="#555555")
axB.set_yticks([]); axB.spines["left"].set_visible(False)
axB.set_xticks([-2, 0, 2, 4, 6, 8])
axB.set_xlabel("GC vs vehicle log2FC (95% CI)", fontsize=FS_AX, labelpad=2)
axB.plot([-1.7], [-2.45], "o", ms=5.2, mfc=GC_CLASS["induced"], mec=TEXT, mew=1.4, clip_on=False)
axB.text(-1.1, -2.45, "positive control (by construction)", fontsize=FS_SM, va="center", ha="left")
SRC["6B"] = pd.DataFrame(srcB)
for g in ["ZBTB16", "PDK4"]:
    reg("B", f"{g} {mt.loc[g].meta_log2FC:+.2f}", round(mt.loc[g].meta_log2FC, 4), "S10b / meta_all_genes", "meta_log2FC", "yes", "T")
for g in ["FKBP5", "DDIT4", "KLF9", "TSC22D3", "TXNIP", "FOXO3", "FGD4", "SGK1", "MYOF", "ZMYND8", "PTPN13"]:
    reg("B", f"{g} {mt.loc[g].meta_log2FC:+.2f}", round(mt.loc[g].meta_log2FC, 4), "S10b / meta_all_genes", "meta_log2FC",
        "yes" if g in ("FKBP5", "SGK1", "MYOF") and False else "new", "S")
reg("B", "group sizes (9), (5), source counts (4), (2)", "9;5;4;2", "S10b", "counts", "new", "N")

# =====================================================================
# C  per-dataset heatmap aligned to B
# =====================================================================
xC, yC, wC, hC = 76.0, 62.0, 46.0, 100.0
header("C", xC + 1.0, yC + 1.0, "Per dataset")
fig.text((xC + 5.6) / W, 1 - (yC + 6.0) / H, "rows as in B", fontsize=FS_SM, va="top", ha="left", **POSTHOC_TAG)
D4 = [("GSE85333", "PRED Mφ"), ("GSE112101", "MP Mo"), ("E-MEXP-2598", "DEX Mo"), ("E-MTAB-7572", "DEX Mφ")]
axC = axmm(xC + 1.5, TOPY, 4 * 10.3, NROW * ROW_H)
cmapC = diverging_cmap(); normC = mpl.colors.TwoSlopeNorm(vmin=-4, vcenter=0, vmax=4)
pv = pdl[pdl.dataset.isin([d for d, _ in D4])].pivot_table(index="gene", columns="dataset", values="logFC")
srcC = []
for y, kind, lab in rows:
    if kind == "hdr":
        continue
    for j, (d, _) in enumerate(D4):
        v = pv.loc[lab, d] if (lab in pv.index and d in pv.columns) else np.nan
        if np.isnan(v):
            axC.add_patch(Rectangle((j, y - 0.5), 1, 1, fc=GC_CLASS["not_measured"], ec="#BDBDBD", lw=0, hatch="//////"))
        else:
            axC.add_patch(Rectangle((j, y - 0.5), 1, 1, fc=cmapC(normC(np.clip(v, -4, 4))), ec="white", lw=0.6))
        srcC.append(dict(gene=lab, dataset=d, logFC=v))
axC.set_xlim(0, 4); axC.set_ylim(NROW - 0.5, -0.5)
axC.set_yticks([]); axC.set_xticks([])
for j, lab in enumerate(["PRED\nMac", "MP\nMono", "DEX\nMono", "DEX\nMac"]):
    axC.text(j + 0.5, 0.42, lab, fontsize=FS_SM, ha="center", va="bottom", linespacing=0.95)
for s in axC.spines.values():
    s.set_visible(False)
caxC = axmm(xC + 1.5, TOPY + NROW * ROW_H + 1.3, 28.0, 2.0)
cb = mpl.colorbar.ColorbarBase(caxC, cmap=cmapC, norm=normC, extend="both", ticks=[-4, 0, 4], orientation="horizontal")
cb.outline.set_linewidth(0.4); caxC.tick_params(labelsize=FS_SM, length=2, pad=1)
cb.set_label("log2FC", fontsize=FS_SM, labelpad=1)
axCk = axmm(xC + 32.0, TOPY + NROW * ROW_H + 0.1, 13.5, 4.5); axCk.axis("off")
axCk.legend(handles=[Rectangle((0, 0), 1, 1, fc=GC_CLASS["not_measured"], ec="#BDBDBD", lw=0.4, hatch="//////", label="n.d.")],
            loc="upper left", fontsize=FS_SM, handlelength=1.0, handletextpad=0.3, borderaxespad=0)
SRC["6C"] = pd.DataFrame(srcC)

# =====================================================================
# D  genome-wide rank (logit percentile)
# =====================================================================
xD, yD, wD, hD = 122.0, 62.0, 61.0, 100.0
header("D", xD + 1.0, yD + 1.0, "Position in the\ngenome-wide GC response")
m = meta.sort_values("meta_log2FC", ascending=False).reset_index(drop=True)
N = len(m); m["rank"] = np.arange(1, N + 1); m["pct"] = (m["rank"] - 0.5) / N
axD = axmm(xD + 12.0, yD + 13.0, wD - 13.5, 50.0)
axD.set_xscale("logit")
axD.plot(m.pct, m.meta_log2FC, color="#9A9A9A", lw=1.0, zorder=2)
for c, col in [("GC_repressed", GC_CLASS["repressed"]), ("GC_induced", GC_CLASS["induced"])]:
    d = m[m["class"] == c]
    axD.scatter(d.pct, d.meta_log2FC, s=2.5, c=col, lw=0, zorder=3, rasterized=True)
axD.axhline(0, color="#7F7F7F", lw=0.6, ls="--", zorder=0)
lit = [g for g in GENES if g in set(m.gene)]
mm = m.set_index("gene")
# label placement: high genes labelled to the right, low genes to the left
upper = ["ZBTB16", "FKBP5", "TSC22D3", "FGD4"]   # extremes + positive controls; other GC-induced genes as in B
LAB_Y = dict(zip(upper, [6.6, 5.0, 3.4, 1.6]))
for g in upper:
    r = mm.loc[g]
    axD.annotate(g, xy=(r.pct, r.meta_log2FC), xytext=(0.30, LAB_Y[g]), fontsize=FS_SM, style="italic", ha="left",
                 va="center", arrowprops=dict(arrowstyle="-", lw=0.4, color="#888888", shrinkA=0.5, shrinkB=2.0), zorder=6)
lower = {"ZMYND8": (0.04, -0.85), "PTPN13": (0.30, -2.0), "SGK1": (0.70, -3.15), "MYOF": (0.95, -4.3)}
for g, (xt, yt) in lower.items():
    r = mm.loc[g]
    axD.annotate(g, xy=(r.pct, r.meta_log2FC), xytext=(xt, yt), fontsize=FS_SM, style="italic", ha="center", va="center",
                 arrowprops=dict(arrowstyle="-", lw=0.4, color="#888888", shrinkA=0.5, shrinkB=2.0), zorder=6)
for g in lit:
    r = mm.loc[g]
    cls_marker(axD, r.pct, r.meta_log2FC, r["class"], g, ms=4.2, z=5)
axD.set_xlim(0.5 / N * 0.6, 1 - 0.5 / N * 0.6); axD.set_ylim(-5.0, 7.4)
axD.set_xlabel("Rank percentile (%), 13,368\ngenes (logit scale)", fontsize=FS_AX, labelpad=34, linespacing=1.0)
axD.set_ylabel("Meta log2FC", fontsize=FS_AX, labelpad=1)
# reference rugs: resident efferocytic core and MHC-II core
rug = axmm(xD + 12.0, yD + 68.0, wD - 13.5, 8.6)
rug.set_xscale("logit"); rug.set_xlim(axD.get_xlim()); rug.set_ylim(0, 2); rug.axis("off")
for i, (gs, lab, xt, ha) in enumerate([("RES_EFF_core", "RES", 0.3, "left"), ("MHCII_core", "MHC-II", 0.75, "right")]):
    gl = [g for g in S3[S3.gene_set == gs].gene if g in mm.index]
    rug.vlines(mm.loc[gl, "pct"], 1.12 - i, 1.88 - i, color="#444444", lw=0.7)
    rug.text(xt, 1.5 - i, lab, fontsize=FS_SM, ha=ha, va="center", color="#444444")
axD.set_xticks([0.001, 0.1, 0.9, 0.999])
axD.set_xticklabels(["0.1", "10", "90", "99.9"])
axD.xaxis.set_minor_locator(mpl.ticker.NullLocator())
axD.set_yticks([-2, 0, 2, 4, 6])
SRC["6D"] = m[["gene", "meta_log2FC", "class", "rank", "pct"]]
reg("D", "13,368 genes", 13368, "S5_GC_classification / meta_all_genes", "row count", "yes", "T")

# =====================================================================
# E  per-donor strip grid
# =====================================================================
xE, yE, wE, hE = 0, 162.0, 122.0, 78.0
header("E", xE, yE + 1.0, "PDK4, TXNIP and FOXO3 are lower in most HFpEF cases")
genesE = ["PDK4", "TXNIP", "FOXO3", "SGK1"]   # LVAD-signature genes not shown in Fig 3F
colsE = [("SCP3342_HFpEF", "Fibroblast", "SCP3342\nFib", scp), ("SCP3342_HFpEF", "Cardiomyocyte", "SCP3342\nCM", scp),
         ("Hahn2021_RVS_HFpEF", "bulk", "Hahn\nbiopsy", hahn), ("MAGNet_exclOverlap_DCM", "bulk", "MAGNet\nbulk", mag),
         ("GSE290367_DCM", "Endothelial", "GSE290367\nEC", gse)]
cellw, cellh, gapx = 20.0, 9.2, 1.6
x0E, y0E = xE + 15.5, yE + 25.5
srcE = []
rng = np.random.default_rng(7)
for j, (ds, ct, lab, sc) in enumerate(colsE):
    meta_ = sc[sc.cell_type == ct][["unit_id", "case"] + (["study"] if "study" in sc else [])].copy()
    meta_["unit_id"] = meta_.unit_id.astype(str)
    sub = ge[(ge.dataset == ds) & (ge.cell_type == ct) & ge.gene.isin(genesE)].copy()
    sub["unit_id"] = sub.unit_id.astype(str)
    sub = sub.merge(meta_, on="unit_id", how="left")
    assert sub.case.notna().all()
    vals = sub.z_ref.dropna()
    lo, hi = np.floor(vals.min()) - 0.5, np.ceil(vals.max()) + 0.5
    coh = "SCP3342" if ds.startswith("SCP") else "Hahn" if ds.startswith("Hahn") else "MAGNet" if ds.startswith("MAG") else "GSE290367"
    fig.text((x0E + j * (cellw + gapx) + cellw / 2) / W, 1 - (y0E - 1.0) / H, lab, fontsize=FS_SM, ha="center", va="bottom",
             color=TEXT, linespacing=1.0)
    fig.patches.append(Rectangle(((x0E + j * (cellw + gapx)) / W, 1 - (y0E - 0.6) / H), cellw / W, 0.6 / H,
                                 transform=fig.transFigure, fc=COHORT[coh], ec="none"))
    for i, g in enumerate(genesE):
        ax = axmm(x0E + j * (cellw + gapx), y0E + i * (cellh + 0.6), cellw, cellh)
        d = sub[sub.gene == g]
        ax.set_xlim(lo, hi); ax.set_ylim(-0.5, 1.5)
        ax.set_yticks([]); ax.spines["left"].set_visible(False)
        if d.z_ref.notna().sum() == 0:
            ax.add_patch(Rectangle((lo, -0.5), hi - lo, 2, fc=GC_CLASS["not_measured"], ec="#BDBDBD", lw=0, hatch="//////"))
            ax.text((lo + hi) / 2, 0.5, "not measured", fontsize=FS_SM, ha="center", va="center", color="#555555",
                    bbox=dict(fc=GC_CLASS["not_measured"], ec="none", pad=0.6))
        else:
            ax.axvline(0, color="#BBBBBB", lw=0.5, zorder=0)
            for case, yv in [(0, 1.0), (1, 0.0)]:
                dd = d[(d.case == case) & d.z_ref.notna()]
                if case == 1 and coh == "GSE290367":
                    cols = dd.study.map(lambda s: GSE_STUDY[s.split()[0]]).values
                else:
                    cols = CONTROL if case == 0 else COHORT[coh]
                ax.scatter(dd.z_ref, yv + rng.uniform(-0.28, 0.28, len(dd)), s=3.2, c=cols, lw=0.25,
                           edgecolors="#6F6F6F" if case == 0 else ("#555555" if coh == "Hahn" else "none"), zorder=3,
                           rasterized=True)
                med = dd.z_ref.median()
                ax.plot([med, med], [yv - 0.42, yv + 0.42], color=TEXT, lw=1.0, zorder=4)
        stp = 5 if hi - lo > 7 else 2
        tks = [t for t in np.arange(-20, 21, stp) if lo + 0.3 <= t <= hi - 0.3]
        ax.set_xticks(tks)
        if i < len(genesE) - 1:
            ax.set_xticklabels([])
            ax.tick_params(axis="x", length=1.5)
        if j == 0:
            ax.text(-0.04, 0.5, g, fontsize=FS_SM, style="italic", ha="right", va="center", transform=ax.transAxes)
            ax.plot(-0.74, 0.5, "o", ms=3.6, transform=ax.transAxes, clip_on=False,
                    mfc=CLS_COL[GCLS[g]], mec="none")
        d2 = d.copy(); d2.insert(0, "column", lab); srcE.append(d2)
fig.text((x0E + 2.5 * (cellw + gapx) - gapx / 2) / W, 1 - (y0E + len(genesE) * (cellh + 0.6) + 4.7) / H,
         "Expression, z vs control group (per donor)", fontsize=FS_AX, ha="center", va="top")
axEk = axmm(xE + 5.6, yE + 6.6, 116.0, 4.0); axEk.axis("off")
axEk.legend(handles=[Line2D([], [], marker="o", ls="none", mfc=CONTROL, mec="#6F6F6F", mew=0.3, ms=3.5, label="Ctrl (upper row)"),
                     Line2D([], [], marker="o", ls="none", mfc=COHORT["SCP3342"], mec="none", ms=3.5, label="Case (lower row)"),
                     Line2D([], [], color=TEXT, lw=1.0, label="median")],
            loc="center left", ncol=3, fontsize=FS_SM, handletextpad=0.3, columnspacing=1.2, borderaxespad=0, handlelength=0.9)
axEk2 = axmm(xE + 5.6, yE + 11.0, 116.0, 4.0); axEk2.axis("off")
axEk2.legend(handles=[Line2D([], [], ls="none", label="GSE290367 cases:"),
                     Line2D([], [], marker="o", ls="none", mfc=GSE_STUDY["Chaffin"], mec="none", ms=3.5, label="Chaffin"),
                     Line2D([], [], marker="o", ls="none", mfc=GSE_STUDY["Koenig"], mec="none", ms=3.5, label="Koenig"),
                     Line2D([], [], marker="o", ls="none", mfc=GSE_STUDY["Reichart"], mec="none", ms=3.5, label="Reichart")],
            loc="center left", ncol=4, fontsize=FS_SM, handletextpad=0.3, columnspacing=1.2, borderaxespad=0, handlelength=0.9)
SRC["6E"] = pd.concat(srcE)[["column", "dataset", "cell_type", "unit_id", "gene", "z_ref", "log_expr", "case"] +
                           (["study"] if "study" in pd.concat(srcE) else [])]

# =====================================================================
# F  schematic: two readings
# =====================================================================
xF, yF, wF, hF = 122.0, 162.0, 61.0, 78.0
header("F", xF + 1.0, yF + 1.0, "Two readings of an LVAD\n'GC receptor signature'")
axF = axmm(xF + 1.0, yF + 11.0, wF - 1.5, hF - 11.0); axF.set_xlim(0, wF - 1.5); axF.set_ylim(hF - 11.0, 0); axF.axis("off")
bw = (wF - 1.5 - 1.6) / 2
BOXH = 42.5
for i, (t1, t2) in enumerate([("LVAD unloading", "(as published)\nmetabolic and\nanti-\ninflammatory\nadaptation"),
                              ("GC exposure", "(alternative;\nnot recorded)\nperioperative\ncortico-\nsteroids at\ntransplant, or\nGC-like imprint\nin donor\ncontrols")]):
    bx = i * (bw + 1.6)
    axF.add_patch(FancyBboxPatch((bx, 0.2), bw, BOXH, boxstyle="round,pad=0,rounding_size=1.2", fc="#F4F4F4", ec="#7F7F7F", lw=0.7))
    axF.text(bx + bw / 2, 1.0, t1, fontsize=FS_SM, fontweight="bold", ha="center", va="top")
    axF.text(bx + bw / 2, 5.3, t2, fontsize=FS_SM, ha="center", va="top", linespacing=1.0)
    axF.add_patch(FancyArrowPatch((bx + bw / 2, BOXH + 0.5), ((wF - 1.5) / 2 + (i - 0.5) * 8, BOXH + 3.6), arrowstyle="-|>",
                                  mutation_scale=7, lw=0.8, color=TEXT))
axF.add_patch(FancyBboxPatch((1.0, BOXH + 3.8), wF - 3.5, 9.4, boxstyle="round,pad=0,rounding_size=1.2", fc="white",
                             ec=GC_CLASS["induced"], lw=1.0))
axF.text((wF - 1.5) / 2, BOXH + 8.5, "Signature up at transplant\n(8 of 9 GC-induced in vitro)", fontsize=FS_SM,
         ha="center", va="center", linespacing=1.0)
axF.text((wF - 1.5) / 2, BOXH + 14.2, "Our data cannot decide\nbetween these.", fontsize=FS_SM, linespacing=1.0, style="italic", va="top",
         ha="center", color="#444444")
SRC["6F"] = pd.DataFrame({"note": ["Schematic; no data. 'Our data cannot decide between these explanations' (Discussion)."]})
reg("F", "8 of 9 GC-induced in vitro", "8/9", "S10b", "count", "yes", "T")


for p, t in [("B", "x -2..8"), ("C", "colour bar -4, 0, 4"), ("D", "x 0.1, 10, 90, 99.9 (%); y -2..6"),
             ("E", "x -10, -5, 0 per column")]:
    reg(p, "axis tick labels: " + t, "", "axis", "axis scale only", "new", "N")
for t, v, src_, st, lk in [("4 recent studies", 4, "S10b sources", "yes", "T"), ("4 human monocyte/macrophage datasets", 4, "S2", "yes", "T"),
                          ("13,368 classified genes", 13368, "S5", "yes", "T"), ("SCP3342 24 controls; 18 Fib, 19 CM", "24;18;19", "perdonor_scores_SCP3342", "yes", "T"),
                          ("Hahn 24/41", "24;41", "perdonor_scores_Hahn", "yes", "T"), ("MAGNet 102/128", "102;128", "perdonor_scores_MAGNet", "yes", "T"),
                          ("GSE290367 EC 45/67", "45;67", "perdonor_scores_GSE290367 / S9b", "new", "S")]:
    reg("legend", t, v, src_, "legend text", st, lk)
# ---------------- geometric QA ----------------
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
    small = sorted({round(t.get_fontsize(), 1) for t in T if t.get_fontsize() < 8})
    return out, ov, small


_o, _v, _s = qa(fig)
print("outside:", _o); print("overlaps:", _v)
_mf = min_font_audit(fig, 10); print("min_font_audit:", _mf)
assert _mf == [], _mf
save_figure(fig, "Figure_6", tight=True)
from PIL import Image  # JACC: RGB TIFF, 600 dpi, LZW
Image.open("Figure_6.tif").convert("RGB").save("Figure_6.tif", compression="tiff_lzw", dpi=(600, 600))
pd.DataFrame(REG).to_csv("Figure_6_number_registry.csv", index=False)
with pd.ExcelWriter("Figure_6_source_data.xlsx") as xw:
    for k, v in SRC.items():
        v.to_excel(xw, sheet_name=k[:31], index=False)
print("saved Figure_6", FONT)
