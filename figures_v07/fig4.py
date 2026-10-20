"""Figure 4 (v0.7, reader-facing labels only) (R3) - Circadian and cAMP transcriptional scores do not account for the imprint.
Standalone: reproduces Figure_4.{pdf,tif,png}, Figure_4_source_data.xlsx and Figure_4_number_registry.csv
from the input artifacts (paths below). Display-level derivations only (ROC curves, ratios, counts); no new tests.
"""
import os, tarfile, numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.colors import to_rgb

import os  # public package: inputs resolved relative to this script (override with JBTS_INPUTS)
HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.environ.get("JBTS_INPUTS", os.path.join(HERE, "..", "figures_v05", "inputs"))
STYLE = os.path.join(HERE, "jbts_style.py")
TAR_FI = os.path.join(INP, "figure_inputs_v05.tar.gz")
TAR_PD = os.path.join(INP, "figure_inputs_v05_perdonor.tar.gz")
XLS = os.path.join(INP, "Supplemental_Tables_S1-S14.xlsx")

exec(open(STYLE).read())
FAM = apply_jbts_style()
mpl.rcParams["font.family"] = [FAM, "DejaVu Sans"]   # glyph fallback for superscript minus/zero
mpl.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": FAM, "mathtext.it": FAM + ":italic",
                     "mathtext.bf": FAM + ":bold", "mathtext.default": "regular", "hatch.linewidth": 0.4})
FS_T, FS_S, FS_G = SIZES["axis_label"], SIZES["tick"], SIZES["min_any"]   # v0.6: 10.5 / 10 / 10 pt; titles 11 pt
os.makedirs("_inputs", exist_ok=True)
for t in (TAR_FI, TAR_PD):
    root = os.path.basename(t).split("_", 1)[1].replace(".tar.gz", "")
    if not os.path.isdir(f"_inputs/{root}"):
        with tarfile.open(t) as tf:
            tf.extractall("_inputs")
PD = "_inputs/figure_inputs_v05_perdonor/v13_rerun/"
FI = "_inputs/figure_inputs_v05/"

# ------------------------------------------------------------------ helpers
FIG_W, FIG_H = 183.0, 240.0
REG, SRC = [], {}
CT = ["Cardiomyocyte", "Fibroblast", "Endothelial", "Pericyte", "Myeloid"]
CTS = {"Cardiomyocyte": "CM", "Fibroblast": "Fib", "Endothelial": "EC", "Pericyte": "Peri", "Myeloid": "Mye"}
STUDIES = ["Chaffin 2022", "Koenig 2022", "Reichart 2022"]
SNAME = {s: s.split()[0] for s in STUDIES}


def m(v, nd=2, sign=True):
    s = f"{v:+.{nd}f}" if sign else f"{v:.{nd}f}"
    return s.replace("-", "\u2212")


SUP = str.maketrans("0123456789-", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b")


def fmt_p(p):
    """P in the text's style; exponents as inline Unicode superscripts at the base size (no smaller text)."""
    if p < 1e-3:
        e = int(np.floor(np.log10(p))); mant = p / 10 ** e
        if round(mant, 1) >= 10: mant, e = mant / 10, e + 1
        return f"{mant:.1f}\u00d710" + str(e).translate(SUP)
    if p >= 0.1: return f"{p:.2f}"
    d = int(-np.floor(np.log10(p))) + 1
    return f"{p:.{d}f}"


def _fmt_p_old(p):
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
    return ax.text(x, y, fmt_p(p), fontsize=fs, ha=ha, va=va, transform=transform or ax.transData, **kw)


def save_exact(fig, stem):
    """Exact canvas (183 mm wide), no tight cropping: PDF, PNG preview (300 dpi), TIFF 600 dpi LZW."""
    fig.savefig(f"{stem}.pdf")
    fig.savefig(f"{stem}.png", dpi=300)
    fig.savefig(f"{stem}.tif", dpi=600, pil_kwargs={"compression": "tiff_lzw"})


def reg(panel, text, value, source, computation, status, lock):
    REG.append(dict(figure="Figure 4", panel=panel, printed_text=text.replace("$", ""), value=value,
                    source_file=source, computation=computation, lock_status=status, matches_locked_text=lock))


def ax_mm(x, y, w, h, **kw):
    """axes from mm coordinates measured from the top-left corner of the canvas"""
    return fig.add_axes([x / FIG_W, 1 - (y + h) / FIG_H, w / FIG_W, h / FIG_H], **kw)


def panel_frame(letter, x, y, w, h, title, tag=None):
    fr = ax_mm(x, y + 6, w, h - 6); fr.axis("off")
    panel_label(fr, letter, x=0, y=1.0)
    fr.text(6.2 / w, 1.0 + 0.4 / (h - 6), title, transform=fr.transAxes, fontsize=SIZES["panel_title"], fontweight="bold",
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
s3 = pd.read_excel(XLS, "S3_gene_sets")
s7b = pd.read_excel(XLS, "S7b_G2_G5")
scp = pd.read_csv(PD + "perdonor_scores_SCP3342_HFpEF.csv")
mag = pd.read_csv(PD + "perdonor_scores_MAGNet_exclOverlap_DCM.csv")
gse = pd.read_csv(PD + "perdonor_scores_GSE290367_DCM.csv")
ge = pd.read_csv(PD + "perdonor_gene_expression_GC8_clock10.csv.gz")
import json
dec = json.load(open(FI + "04_results/validation_v1.3/v13_decisions.json"))

fig = plt.figure(figsize=(FIG_W * MM, FIG_H * MM))
CW = (FIG_W - 3 - 11 * 3) / 12.0
def cx(c): return 1.5 + c * (CW + 3)
def cw(n): return n * CW + (n - 1) * 3
ROW_Y = [0, 74, 152]; ROW_H = [74, 78, 88]
def ftext(xmm, ymm, t, **kw): fig.text(xmm / FIG_W, 1 - ymm / FIG_H, t, **kw)
CLS_COL = {"GC_induced": GC_CLASS["induced"], "GC_repressed": GC_CLASS["repressed"], "GC_insensitive": GC_CLASS["insensitive"],
           "ambiguous": GC_CLASS["ambiguous"], np.nan: GC_CLASS["not_measured"]}
ADJ_MK = {"base": "o", "+CLOCK": "^", "+CAMP": "s", "+CLOCK+CAMP": "D"}
ADJ_LAB = {"base": "Base", "+CLOCK": "+ circadian", "+CAMP": "+ cAMP", "+CLOCK+CAMP": "+ both"}
DS = {"SCP3342": "SCP3342_HFpEF", "MAGNet": "MAGNet_exclOverlap_DCM", "GSE290367": "GSE290367_DCM"}

def adj_mark(ax, x, y, coh, adj, ms=4.0, zorder=4):
    """Deck grammar: shape = cohort; fill = adjustment (filled base, open +circadian, half +cAMP, dotted outline both)."""
    mk = COHORT_MARKER[coh]; colr = COHORT[coh]
    if mk == "^": ms = ms * 1.1
    if adj == "base":
        ax.plot(x, y, marker=mk, ms=ms, mfc=colr, mec=colr, mew=0.8, ls="none", zorder=zorder)
    elif adj == "+CLOCK":
        ax.plot(x, y, marker=mk, ms=ms, mfc="white", mec=colr, mew=0.9, ls="none", zorder=zorder)
    elif adj == "+CAMP":
        ax.plot(x, y, marker=mk, ms=ms, fillstyle="left", mfc=colr, markerfacecoloralt="white", mec=colr, mew=0.9, ls="none", zorder=zorder)
    else:
        ax.scatter([x], [y], s=(ms * 1.05) ** 2, marker=mk, facecolors="white", edgecolors=colr, linewidths=0.9,
                   linestyles=(0, (1.0, 1.0)), zorder=zorder)

def s8row(ds, test, ct, mod, adj):
    q = s8[(s8.dataset == ds) & (s8.test == test) & (s8.cell_type == ct) & (s8.module == mod) & (s8.adjustment == adj)]
    assert len(q) == 1, (ds, test, ct, mod, adj, len(q)); return q.iloc[0]

def swatch(ax, x, y, cls, s=2.0, tr=None):
    fc = CLS_COL.get(cls, GC_CLASS["not_measured"]) if isinstance(cls, str) else GC_CLASS["not_measured"]
    ec = AMBIGUOUS_EDGE if cls == "ambiguous" else ("#9A9A9A" if not isinstance(cls, str) else fc)
    p = Rectangle((x, y - s / 2), s, s, fc=fc, ec=ec, lw=0.5, hatch="////" if not isinstance(cls, str) else None,
                  transform=tr or ax.transData, clip_on=False)
    ax.add_patch(p)

# ================================================================== A  schematic
x0, y0, w, h = 1.5, ROW_Y[0], 90.0, ROW_H[0]
panel_frame("A", x0, y0, w, h, "Three competing programs")
sa = ax_mm(x0, y0 + 6, w, h - 6); sa.axis("off"); sa.set_xlim(0, w); sa.set_ylim(h - 6, 0)   # axes in mm, y down
progs = [("GC", "GC", "proxy for\nGC exposure", [("", "GC_specific3")]),
         ("CLK", "Circadian", "proxy for\nsampling time", [("morning", "CLOCK_morning"), ("evening", "CLOCK_evening")]),
         ("CAMP", "cAMP", "proxy for\ncatecholamines", [("", "CAMP_IEG")])]
colw = [26.5, 36.0, 26.0]; G_ = 0.75
colx = [0.0, 26.5 + G_, 26.5 + G_ + 36.0 + G_]
LN = 3.75; BOXH = 46.5
srcA = []
for (key, name, proxy, sets), xc, wc in zip(progs, colx, colw):
    sa.add_patch(FancyBboxPatch((xc, 0.4), wc, BOXH, boxstyle="round,pad=0,rounding_size=1.2", fc="#F4F4F4", ec="#BDBDBD", lw=0.6))
    ix, iy = xc + 3.0, 3.6
    if key == "GC":      # steroid nucleus: three hexagons + pentagon
        for (cx_, cy_) in [(ix - 1.3, iy + 0.5), (ix + 0.3, iy - 0.4), (ix + 1.9, iy + 0.5)]:
            sa.add_patch(mpl.patches.RegularPolygon((cx_, cy_), 6, radius=0.95, orientation=0, fc="none", ec=GC_CLASS["induced"], lw=0.8))
        sa.add_patch(mpl.patches.RegularPolygon((ix + 3.4, iy - 0.4), 5, radius=0.8, orientation=np.pi, fc="none", ec=GC_CLASS["induced"], lw=0.8))
        tx_ = xc + 8.8
    elif key == "CLK":   # clock face
        sa.add_patch(mpl.patches.Circle((ix, iy), 1.9, fc="white", ec=TEXT, lw=0.8))
        sa.plot([ix, ix], [iy, iy - 1.4], color=TEXT, lw=0.8); sa.plot([ix, ix + 1.0], [iy, iy + 0.45], color=TEXT, lw=0.8)
        tx_ = xc + 6.0
    else:                # infusion bag with line (vasopressor / inotrope)
        sa.add_patch(FancyBboxPatch((ix - 1.3, iy - 2.0), 2.6, 3.0, boxstyle="round,pad=0,rounding_size=0.6", fc="white", ec=TEXT, lw=0.8))
        sa.plot([ix, ix, ix + 1.8], [iy + 1.0, iy + 2.0, iy + 2.0], color=TEXT, lw=0.8)
        tx_ = xc + 6.0
    sa.text(tx_, 3.4, name, fontsize=FS_S, fontweight="bold", va="center", ha="left")
    sa.text(xc + 1.0, 7.2, proxy, fontsize=FS_G, va="top", ha="left", style="italic", color="#555555", linespacing=1.0)
    for k, (lab, gs) in enumerate(sets):
        gx = xc + 1.0 + k * 19.2
        gy = 17.4
        if lab:
            sa.text(gx, gy, lab, fontsize=FS_G, va="center", ha="left", color="#555555"); gy += LN
        for gene, cls in s3[s3.gene_set == gs][["gene", "GC_class"]].values:
            swatch(sa, gx, gy, cls, s=2.3)
            sa.text(gx + 3.0, gy, gene, fontsize=FS_G, style="italic", va="center", ha="left")
            srcA.append(dict(program=name, gene_set=gs, gene=gene, GC_class=cls if isinstance(cls, str) else "not measured"))
            gy += LN
# in vitro class key inside the GC box (below its 3 genes)
sa.text(colx[0] + 1.0, 30.0, "in vitro class", fontsize=FS_G, va="center", ha="left", color="#555555")
for k, (cls, lab) in enumerate([("GC_induced", "induced"), ("GC_insensitive", "insensitive"), ("ambiguous", "ambiguous"), (np.nan, "not measured")]):
    yk = 30.0 + LN * (k + 1)
    swatch(sa, colx[0] + 1.0, yk, cls, s=2.3)
    sa.text(colx[0] + 4.0, yk, lab, fontsize=FS_G, va="center", ha="left")
sa.add_patch(mpl.patches.FancyArrowPatch((2, BOXH + 3.4), (w - 2, BOXH + 3.4), arrowstyle="-|>", mutation_scale=8, lw=0.9, color=TEXT))
sa.text(w / 2, BOXH + 5.2, "add circadian or cAMP score as a covariate:\nGC difference retained (E); GC-induced\nefferocytic genes stay lower (G)", fontsize=FS_G,
        ha="center", va="top", linespacing=1.05)
SRC["4A"] = pd.DataFrame(srcA)

# ================================================================== B  clock genes heatmap
x0, y0, w, h = 94.5, ROW_Y[0], 87.0, ROW_H[0]
panel_frame("B", x0, y0, w, h, "Clock genes, per donor")
CLK = ["NR1D2", "DBP", "TEF", "PER2", "BHLHE41", "CIART", "ARNTL", "NPAS2", "CLOCK", "CRY1"]
blocksB = [("SCP3342", "SCP3342_HFpEF", "Fibroblast", scp, "SCP3342\nFib", "HFpEF"),
           ("MAGNet", "MAGNet_exclOverlap_DCM", "bulk", mag, "MAGNet\nbulk", "DCM")]
gl, gr = 24.0, 12.0
bw_tot = w - gl - gr; gap = 2.5
fr_ = [0.45, 0.55]          # block widths (display)
ytr = y0 + 6 + 8.2; trh = 6.0
hm_top = ytr + 3.3 + 0.6 + trh + 0.8; hm_h = h - (hm_top - y0) - 5.0
srcB = []
bx = x0 + gl
for k, (coh, ds, ctp, sdf, lab, cl) in enumerate(blocksB):
    bw = (bw_tot - gap) * fr_[k]
    sd = sdf[sdf.cell_type == ctp][["unit_id", "case", "score_CLOCK"]].copy()
    sd = sd.sort_values(["case", "score_CLOCK"], ascending=[True, True]); sd["unit_id"] = sd.unit_id.astype(str)
    sub = ge[(ge.dataset == ds) & (ge.cell_type == ctp)].copy(); sub["unit_id"] = sub.unit_id.astype(str)
    sub = sub.pivot(index="gene", columns="unit_id", values="z_ref")
    Z = sub.loc[CLK, sd.unit_id].values
    nc = int((sd.case == 0).sum())
    ah = ax_mm(bx, hm_top, bw, hm_h)
    im = ah.imshow(np.clip(Z, -CAP, CAP), cmap=GCDIV, vmin=-CAP, vmax=CAP, aspect="auto", interpolation="nearest")
    for i in range(len(CLK)):
        js = np.where(np.isnan(Z[i]))[0]
        if len(js):
            st = pv = js[0]
            for jj in list(js[1:]) + [None]:
                if jj is None or jj != pv + 1:
                    ah.add_patch(Rectangle((st - .5, i - .5), pv - st + 1, 1, fc=GC_CLASS["not_measured"], ec="#9A9A9A", hatch="//////", lw=0))
                    if jj is not None: st = jj
                if jj is not None: pv = jj
    ah.axhline(5.5, color="white", lw=1.6); ah.axvline(nc - 0.5, color="white", lw=1.2)
    ah.set_xticks([]); [s.set_visible(False) for s in ah.spines.values()]
    if k == 0:
        ah.set_yticks(range(len(CLK))); ah.set_yticklabels(CLK, fontsize=FS_G, fontstyle="italic"); ah.tick_params(axis="y", length=0, pad=1.5)
    else:
        ah.set_yticks([])
    ag = ax_mm(bx, ytr, bw, 3.3)
    ag.imshow([[to_rgb(CONTROL) if c == 0 else to_rgb(COHORT[coh]) for c in sd.case]], aspect="auto", interpolation="nearest")
    ag.axvline(nc - 0.5, color="white", lw=1.2)
    ab = ax_mm(bx, ytr + 3.9, bw, trh)
    ab.bar(np.arange(len(sd)), sd.score_CLOCK.values, width=0.9, color=[CONTROL if c == 0 else COHORT[coh] for c in sd.case], lw=0)
    ab.axhline(0, color="#7F7F7F", lw=0.5); ab.set_xlim(-0.5, len(sd) - 0.5); ab.set_ylim(-4.5, 4.5); ab.axvline(nc - 0.5, color="white", lw=1.2)
    for a in (ag, ab):
        a.set_xticks([]); [s.set_visible(False) for s in a.spines.values()]
        a.set_yticks([])
    if k == 0:
        ag.set_yticks([0]); ag.set_yticklabels(["Group"], fontsize=FS_G); ag.tick_params(axis="y", length=0, pad=1.5)
        ab.set_yticks([0]); ab.set_yticklabels(["Circadian"], fontsize=FS_G); ab.tick_params(axis="y", length=0, pad=1.5)
    fig.text((bx + bw / 2) / FIG_W, 1 - (ytr - 0.7) / FIG_H, lab, fontsize=FS_S, fontweight="bold", ha="center", va="bottom",
             color=COHORT[coh], linespacing=1.0)
    fig.text((bx + bw / 2) / FIG_W, 1 - (y0 + h - 0.3) / FIG_H, f"{nc}/{len(sd) - nc}", fontsize=FS_G, ha="center", va="bottom")
    reg("B", f"{nc}/{len(sd) - nc}", f"{nc}/{len(sd) - nc}", ds, "n per group in block", "T" if coh == "MAGNet" else "S",
        "yes" if coh == "MAGNet" else "new")
    long = sub.loc[CLK, sd.unit_id].reset_index().melt(id_vars="gene", var_name="unit_id", value_name="z_ref").merge(sd, on="unit_id")
    long.insert(0, "block", f"{coh}_{ctp}"); srcB.append(long)
    bx += bw + gap
ftext(x0 + gl - 1.5, y0 + h - 0.3, "n", fontsize=FS_G, ha="right", va="bottom", style="italic")
# arm brackets
aa = ax_mm(x0 + 0.5, hm_top, 4.0, hm_h); aa.axis("off"); aa.set_ylim(9.5, -0.5); aa.set_xlim(0, 1)
for (a0, a1, lab) in ((0, 5, "morning"), (6, 9, "evening")):
    aa.plot([0.95, 0.75, 0.75, 0.95], [a0 - .35, a0 - .35, a1 + .35, a1 + .35], color=TEXT, lw=0.6, clip_on=False)
    aa.text(0.35, (a0 + a1) / 2, lab, rotation=90, fontsize=FS_G, ha="center", va="center")
cax = ax_mm(x0 + w - gr + 1.5, hm_top, 2.2, hm_h)
cb = fig.colorbar(im, cax=cax, ticks=[-3, 0, 3]); cb.outline.set_linewidth(0.4)
cb.ax.set_yticklabels(["\u2264\u22123", "0", "\u2265+3"]); cb.ax.tick_params(labelsize=FS_G, length=2)
cb.ax.set_title("z, SD", fontsize=FS_G, pad=6, loc="left")
SRC["4B"] = pd.concat(srcB)

# ================================================================== C  circadian and cAMP by group
x0, y0, w, h = cx(0), ROW_Y[1], cw(8), ROW_H[1]
panel_frame("C", x0, y0, w, h, "Circadian and cAMP scores by group")
facets = [("SCP3342", scp, ct, CTS[ct]) for ct in CT] + [("MAGNet", mag, "bulk", "MAGNet"), ("GSE290367", gse, "Endothelial", "GSE EC")]
lw_ = 18.5; gapC = 1.4; extra = 2.2
fwC = (w - lw_ - 6 * gapC - 2 * extra) / 7
rowsC = [("score_CLOCK", "CLOCK", "Circadian"), ("score_CAMP", "CAMP", "cAMP")]
topC = 12.5; rh = (h - 6 - topC - 5.0 - 5.5) / 2
srcC = []
for ri, (col_, mod, rname) in enumerate(rowsC):
    for fi, (coh, d0, ct, lab) in enumerate(facets):
        fx = x0 + lw_ + fi * (fwC + gapC) + (extra if fi >= 5 else 0) + (extra if fi >= 6 else 0)
        fy = y0 + 6 + topC + ri * (rh + 5.5)
        a = ax_mm(fx, fy, fwC, rh)
        d = d0[d0.cell_type == ct]
        for g, colr in ((0, CONTROL), (1, COHORT[coh])):
            v = d.loc[d.case == g, col_].dropna().values
            a.boxplot(v, positions=[g], widths=0.6, showfliers=False, medianprops=dict(color=TEXT, lw=1.0),
                      boxprops=dict(color=TEXT, lw=0.5), whiskerprops=dict(color=TEXT, lw=0.5), capprops=dict(color=TEXT, lw=0.5), zorder=1)
            a.scatter(g + rng.uniform(-0.22, 0.22, len(v)), v, s=3.5 if coh == "MAGNet" else 5, marker=COHORT_MARKER[coh], fc=colr, ec="none",
                      alpha=0.8, zorder=2)
            srcC.append(pd.DataFrame(dict(score=rname, cohort=coh, cell_type=ct, unit_id=d.loc[d.case == g].dropna(subset=[col_]).unit_id.values,
                                          case=g, value=v)))
        r = s8row(DS[coh], "C1/descriptive", ct, mod, "base")
        a.text(0.5, 1.02, m(r.beta), transform=a.transAxes, fontsize=FS_G, ha="center", va="bottom")
        if fi == 0:
            a.text(-0.12, 1.02, "case\u2212ctrl", transform=a.transAxes, fontsize=FS_G, ha="right", va="bottom", color="#555555")
        lockT = (mod == "CLOCK" and coh == "MAGNet") or (mod == "CAMP" and coh == "GSE290367")
        reg("C", m(r.beta), round(r.beta, 4), f"Table S8a (descriptive circadian/cAMP score, {DS[coh]}, {ct}, {mod})", "case - control, SD",
            "T" if (mod == "CLOCK" and coh == "MAGNet") else "S", "yes" if (mod == "CLOCK" and coh == "MAGNet") else "new")
        zero_line(a, "h")
        a.set_xlim(-0.6, 1.6); a.set_ylim(-6.5, 6.5); a.set_yticks([-5, 0, 5])
        a.set_xticks([]); a.spines["bottom"].set_visible(True)
        if fi == 0:
            a.set_yticklabels(["\u22125", "0", "5"]); a.set_ylabel(rname + "\nscore (SD)", fontsize=FS_T, labelpad=1, linespacing=1.0)
        else:
            a.set_yticklabels([])
        if ri == 0:
            a.text(0.5, 1.0 + 4.0 / rh, lab, transform=a.transAxes, fontsize=FS_S, ha="center", va="bottom",
                   fontweight="bold" if coh != "SCP3342" else "normal", color=COHORT[coh] if coh != "SCP3342" else TEXT)
# cohort headers
def hdr_span(f0, f1, text, colr):
    xa = x0 + lw_ + f0 * (fwC + gapC) + (extra if f0 >= 5 else 0) + (extra if f0 >= 6 else 0)
    xb = x0 + lw_ + f1 * (fwC + gapC) + fwC + (extra if f1 >= 5 else 0) + (extra if f1 >= 6 else 0)
    yy = y0 + 6 + 1.6
    fig.text(((xa + xb) / 2) / FIG_W, 1 - yy / FIG_H, text, fontsize=FS_S, fontweight="bold", color=colr, ha="center", va="center")
    fig.add_artist(mpl.lines.Line2D([xa / FIG_W, xb / FIG_W], [1 - (yy + 2.2) / FIG_H] * 2, color=colr, lw=0.8))
hdr_span(0, 4, "SCP3342 (HFpEF vs ctrl)", COHORT["SCP3342"])
fig.text((x0 + w) / FIG_W, 1 - (y0 + h - 0.2) / FIG_H, "Left box, control (gray); right box, case (cohort color)",
         fontsize=FS_G, ha="right", va="bottom", color="#555555")
SRC["4C"] = pd.concat(srcC)

# ================================================================== D  ROC + AUC strip
x0, y0, w, h = cx(8), ROW_Y[1], cw(4), ROW_H[1]
panel_frame("D", x0, y0, w, h, "Separation by each score")
SC = [("score_GC3", "GC3", "GC", "-", "o"), ("score_CLOCK", "CLOCK", "Circadian", (0, (3, 1.5)), "^"), ("score_CAMP", "CAMP", "cAMP", (0, (1, 1.2)), "s")]
ar = ax_mm(x0 + 11.0, y0 + 6 + 2.5, 20.0, 24.0)
d = scp[scp.cell_type == "Fibroblast"]
srcD = []
for col_, mod, nm, ls, mk in SC:
    y_true = (d.case == 0).astype(int).values; s_ = d[col_].values
    thr = np.r_[np.inf, np.sort(np.unique(s_))[::-1]]
    tpr = [((s_ >= t) & (y_true == 1)).sum() / y_true.sum() for t in thr]
    fpr = [((s_ >= t) & (y_true == 0)).sum() / (1 - y_true).sum() for t in thr]
    ar.plot(fpr, tpr, color=COHORT["SCP3342"], ls=ls, lw=1.1, drawstyle="steps-post")
    srcD.append(pd.DataFrame(dict(curve=f"ROC_Fibroblast_{nm}", fpr=fpr, tpr=tpr)))
ar.plot([0, 1], [0, 1], color="#BDBDBD", lw=0.6)
ar.set_xlim(-0.02, 1.02); ar.set_ylim(-0.02, 1.02); ar.set_xticks([0, 0.5, 1]); ar.set_yticks([0, 0.5, 1])
ar.set_xticklabels(["0", "0.5", "1"]); ar.set_yticklabels(["0", "0.5", "1"])
ar.set_xlabel("1 \u2212 specificity", fontsize=FS_T, labelpad=0.5); ar.set_ylabel("Sensitivity", fontsize=FS_T, labelpad=0.5)
# line key to the right of the ROC (fibroblasts)
tx = ax_mm(x0 + 33.0, y0 + 6 + 2.5, w - 33.0, 24.0); tx.axis("off"); tx.set_xlim(0, 1); tx.set_ylim(1, 0)
tx.text(0.0, 0.0, "SCP3342\nfibroblasts", fontsize=FS_G, va="top", color=COHORT["SCP3342"], fontweight="bold", linespacing=1.05)
for k, (col_, mod, nm, ls, mk) in enumerate(SC):
    yk = 0.45 + k * 0.17
    tx.plot([0.0, 0.24], [yk, yk], color=COHORT["SCP3342"], ls=ls, lw=1.1)
    tx.text(0.29, yk, nm, fontsize=FS_G, va="center")
# AUC strip: rows = scores, points = 5 cell types
aS = ax_mm(x0 + 17.0, y0 + 6 + 44.0, 22.0, h - 6 - 44.0 - 9.0)
aS.axvspan(0.9, 1.02, color="#EDEDED", lw=0); aS.axvspan(-0.02, 0.1, color="#EDEDED", lw=0)
aS.axvline(0.5, color="#7F7F7F", lw=0.6, ls=(0, (3, 2)), zorder=0)
rowsD = []
SCS = SC[1:]                      # GC AUCs are shown in Fig 5B; strip keeps the two alternative scores
for j, (col_, mod, nm, ls, mk) in enumerate(SCS):
    for i, ct in enumerate(CT):
        r = s8row("SCP3342_HFpEF", "K2_AUC", ct, mod, "none")
        rowsD.append(dict(cell_type=ct, score=nm, AUC=r.beta))
    v = [x["AUC"] for x in rowsD if x["score"] == nm]
    aS.plot(v, [j] * len(v), marker="o", ms=3.6, mfc=COHORT["SCP3342"], mec="white", mew=0.4, ls="none", alpha=0.9, zorder=3)
for i, ct in enumerate(CT):
    r = s8row("SCP3342_HFpEF", "K2_AUC", ct, "GC3", "none"); rowsD.append(dict(cell_type=ct, score="GC (not in strip; Fig 5B)", AUC=r.beta))
DA = pd.DataFrame(rowsD)
aS.set_ylim(1.6, -0.6); aS.set_yticks(range(2)); aS.set_yticklabels([x[2] for x in SCS], fontsize=FS_G); aS.tick_params(axis="y", length=0, pad=1)
aS.set_xlim(-0.02, 1.02); aS.set_xticks([0, 0.5, 1]); aS.set_xticklabels(["0", "0.5", "1"])
aS.set_xlabel("AUC, P(control > HFpEF)", fontsize=FS_T, labelpad=0.5)
aS.set_title("gray: AUC \u22650.9 or \u22640.1", fontsize=FS_G, pad=1.5, color="#555555", loc="center")
rngs = {nm: (DA[DA.score == nm].AUC.min(), DA[DA.score == nm].AUC.max()) for nm in ["Circadian", "cAMP"]}
for j, (col_, mod, nm, ls, mk) in enumerate(SCS):
    t = "{:.2f}".format(rngs[nm][0]) + "\u2013" + "{:.2f}".format(rngs[nm][1])
    aS.text(1.04, j, t, transform=aS.get_yaxis_transform(), fontsize=FS_G, va="center", ha="left")
    st = "S" if nm == "cAMP" else "T"
    reg("D", t, f"{rngs[nm][0]:.4f}..{rngs[nm][1]:.4f}", "XL:S8a K2_AUC SCP3342 (reproduced from per-donor scores)",
        f"range of AUC over 5 cell types ({nm})", st, "yes" if st == "T" else "new")
reg("D", "AUC >=0.9 or <=0.1", "0.9/0.1", "Manuscript Methods (identifiability thresholds)", "illustrative threshold shading", "T", "yes")
SRC["4D"] = pd.concat([DA.assign(curve="AUC")] + srcD)

# ================================================================== E  dumbbell retention
x0, y0, w, h = 1.5, ROW_Y[2], 70.0, ROW_H[2]
panel_frame("E", x0, y0, w, h, "GC difference is mostly retained")
rowsE = [("SCP3342", ct) for ct in CT] + [("MAGNet", "bulk")] + [("GSE290367", ct) for ct in CT]
yv, yy, prev, hdrs = [], 0.0, None, {}
for coh, ct in rowsE:
    if coh != prev:
        if prev is not None: yy += 0.35
        if coh != "MAGNet": hdrs[coh] = yy; yy += 1.0
    yv.append(yy); yy += 1; prev = coh
TOPE = 9.0
ae = ax_mm(x0 + 21.0, y0 + 6 + TOPE, 26.0, h - 6 - TOPE - 9.5)
al = ax_mm(x0, y0 + 6 + TOPE, 20.5, h - 6 - TOPE - 9.5); al.axis("off")
at = ax_mm(x0 + 48.0, y0 + 6 + TOPE, w - 48.0, h - 6 - TOPE - 9.5); at.axis("off")
for a in (ae, al, at): a.set_ylim(max(yv) + 0.6, -0.6)
al.set_xlim(0, 1); at.set_xlim(0, 1)
zero_line(ae)
DEC_CT = ["Fibroblast", "Pericyte", "Myeloid"]
srcE = []
DOD = {"base": -0.27, "+CLOCK": 0.0, "+CAMP": 0.27}
for (coh, ct), y in zip(rowsE, yv):
    test = "C/K_bulk" if coh == "MAGNet" else "C/K"
    vals = {adj: -s8row(DS[coh], test, ct, "GC3", adj).beta for adj in ["base", "+CLOCK", "+CAMP"]}
    colr = COHORT[coh]
    for adj in ["+CLOCK", "+CAMP"]:
        ae.plot([vals["base"], vals[adj]], [y + DOD["base"], y + DOD[adj]], color=colr, lw=0.6, zorder=1)
    for adj in ["base", "+CLOCK", "+CAMP"]:
        adj_mark(ae, vals[adj], y + DOD[adj], coh, adj, ms=3.8, zorder=3 if adj == "base" else 4)
    lab = CTS[ct] if ct != "bulk" else "MAGNet"
    al.text(0.97, y, lab, fontsize=FS_G, ha="right", va="center",
            fontweight="bold" if (coh == "GSE290367" and ct in DEC_CT) or coh == "MAGNet" else "normal",
            color=COHORT["MAGNet"] if coh == "MAGNet" else TEXT)
    srcE.append(dict(cohort=coh, cell_type=ct, base=vals["base"], plus_circadian=vals["+CLOCK"], plus_cAMP=vals["+CAMP"],
                     ret_circadian=vals["+CLOCK"] / vals["base"], ret_cAMP=vals["+CAMP"] / vals["base"],
                     in_decision_rule=(coh == "SCP3342") or (coh == "GSE290367" and ct in DEC_CT)))
EE = pd.DataFrame(srcE)
for coh, yh_ in hdrs.items():
    al.text(0.97, yh_ + 0.1, coh, fontsize=FS_G, fontweight="bold", color=COHORT[coh], ha="right", va="center")
for coh in ["SCP3342", "MAGNet", "GSE290367"]:
    ys = [y for (c, _), y in zip(rowsE, yv) if c == coh]
    sub = EE[EE.cohort == coh]
    if coh == "GSE290367":
        sub = sub[sub.cell_type.isin(DEC_CT)]
        mc, mk_ = dec["C2"]["retention"]["GSE290367_DCM"]["median"], dec["K1"]["retention"]["GSE290367_DCM"]["median"]
        assert abs(np.median(sub.ret_circadian) - mc) < 1e-9 and abs(np.median(sub.ret_cAMP) - mk_) < 1e-9
    elif coh == "SCP3342":
        mc, mk_ = dec["C2"]["retention"]["SCP3342_HFpEF"]["median"], dec["K1"]["retention"]["SCP3342_HFpEF"]["median"]
        assert abs(np.median(sub.ret_circadian) - mc) < 1e-9 and abs(np.median(sub.ret_cAMP) - mk_) < 1e-9
    else:
        mc, mk_ = sub.ret_circadian.iloc[0], sub.ret_cAMP.iloc[0]
    yc = np.mean(ys)
    at.text(0.22, yc, f"{mc * 100:.0f}%", fontsize=FS_S, ha="center", va="center")
    at.text(0.79, yc, f"{mk_ * 100:.0f}%", fontsize=FS_S, ha="center", va="center")
    src_ = ("04_results/v13_decisions.json (circadian and cAMP tests: median retention)" if coh != "MAGNet"
            else "Table S8a (MAGNet bulk, 3-gene GC score: +circadian/base, +cAMP/base)")
    reg("E", f"{mc * 100:.0f}%", round(mc, 4), src_, f"median retention, circadian, {coh}", "T", "yes")
    reg("E", f"{mk_ * 100:.0f}%", round(mk_, 4), src_, f"median retention, cAMP, {coh}", "T", "yes")
# headers (figure coords, above the axes)
ftext(x0 + 48.0 + (w - 48.0) / 2, y0 + 6 + 0.3, "median", fontsize=FS_G, ha="center", va="top", color="#555555")
ftext(x0 + 48.0 + (w - 48.0) * 0.22, y0 + 6 + 4.2, "+circ.", fontsize=FS_G, ha="center", va="top", fontweight="bold")
ftext(x0 + 48.0 + (w - 48.0) * 0.79, y0 + 6 + 4.2, "+cAMP", fontsize=FS_G, ha="center", va="top", fontweight="bold")
for k, (fs_, lab, fc, xk) in enumerate([("full", "base", "#555555", 1.2), ("none", "+ circadian", "white", 12.8), ("left", "+ cAMP", "#555555", 34.0)]):
    fig.add_artist(mpl.lines.Line2D([(x0 + xk) / FIG_W], [1 - (y0 + 6 + 2.2) / FIG_H], marker="o", ms=3.8, mfc=fc, mec="#555555", mew=0.8,
                                    ls="none", fillstyle=fs_ if fs_ != "none" else "full", markerfacecoloralt="white"))
    ftext(x0 + xk + 1.6, y0 + 6 + 2.2, lab, fontsize=FS_G, va="center", ha="left")
yG = [y for (c, _), y in zip(rowsE, yv) if c == "GSE290367"]
ae.set_xlim(-1.6, 7.4); ae.set_xticks([0, 3, 6]); ae.set_yticks([]); ae.spines["left"].set_visible(False)
ae.set_xlabel("Control \u2212 case GC score (SD)", fontsize=FS_T, labelpad=0.5)
SRC["4E"] = EE

# ================================================================== F  GC vs circadian scatter
x0, y0, w, h = 74.5, ROW_Y[2], 42.0, ROW_H[2]
panel_frame("F", x0, y0, w, h, "GC vs circadian")
srcF = []
fh = (h - 6 - 4.0 - 5.0 - 10.5) / 2
for k, (coh, d0, ct, lab) in enumerate([("SCP3342", scp, "Fibroblast", "SCP3342 Fib"), ("MAGNet", mag, "bulk", "MAGNet")]):
    a = ax_mm(x0 + 11.0, y0 + 6 + 4.0 + k * (fh + 5.0), w - 12.0, fh)
    d = d0[d0.cell_type == ct]
    for g, colr in ((0, CONTROL), (1, COHORT[coh])):
        dd = d[d.case == g]
        a.scatter(dd.score_CLOCK, dd.score_GC3, s=5 if coh == "SCP3342" else 3.5, marker=COHORT_MARKER[coh], fc=colr, ec="none", alpha=0.85)
        srcF.append(dd[["unit_id", "case", "score_CLOCK", "score_GC3"]].assign(cohort=coh, cell_type=ct))
    a.set_title(lab, fontsize=FS_S, pad=1.0, color=COHORT[coh], fontweight="bold")
    zero_line(a, "h"); zero_line(a, "v")
    a.set_ylabel("GC score", fontsize=FS_T, labelpad=0.5)
    a.set_xlim(-5, 5); a.set_xticks([-4, 0, 4])
    a.set_xticklabels(["\u22124", "0", "4"] if k == 1 else [])
    yl = (-10.5, 3.5) if coh == "SCP3342" else (-6.5, 4.5)
    a.set_ylim(*yl); a.set_yticks([-8, -4, 0] if coh == "SCP3342" else [-4, 0, 4])
    a.set_yticklabels([m(v, 0, sign=False) for v in a.get_yticks()])
    if k == 1: a.set_xlabel("Circadian score", fontsize=FS_T, labelpad=0.5)
c3 = s8[(s8.test == "C3_rho_in_ref") & (s8.module == "GC3_vs_CLOCK")]
rS = c3[c3.dataset == "SCP3342_HFpEF"].beta; rG = c3[c3.dataset == "GSE290367_DCM"].beta
gmax = rG.max(); gmax_t = "0.00" if abs(gmax) < 0.005 else m(gmax, sign=False)
tF = (f"Within-control \u03c1 (5 cell types):\nSCP3342 {m(rS.min(), sign=False)} to {m(rS.max(), sign=False)}\n"
      f"GSE290367 {m(rG.min(), sign=False)} to {gmax_t}")
reg("F (legend)", f"SCP3342 {m(rS.min(), sign=False)} to {m(rS.max(), sign=False)}", f"{rS.min():.4f}..{rS.max():.4f}", "Table S8a (within-control Spearman rho, 3-gene GC vs circadian score)", "range Spearman rho, controls", "S", "new")
reg("F (legend)", f"GSE290367 {m(rG.min(), sign=False)} to {gmax_t}", f"{rG.min():.4f}..{rG.max():.4f}", "Table S8a (within-control Spearman rho, 3-gene GC vs circadian score)", "range Spearman rho, controls", "S", "new")
SRC["4F"] = pd.concat(srcF)

# ================================================================== G  G4 forest
x0, y0, w, h = 119.5, ROW_Y[2], 62.0, ROW_H[2]
frG = panel_frame("G", x0, y0, w, h, "GC-induced efferocytic genes")
ftext(x0 + 6.2, y0 + 6 + 0.1, "stay lower in DCM", fontsize=SIZES["panel_title"], fontweight="bold", ha="left", va="top")
G4 = [("GSE290367", "GSE290367_Resident_DCM", "Myeloid:Resident", "GSE290367\nresident mac."), ("MAGNet", "MAGNet_exclOverlap_DCM", "bulk", "MAGNet\nbulk")]
yv, rowsG, yy = [], [], 0.0
for coh, ds, ct, lab in G4:
    for adj in ["base", "+CLOCK", "+CAMP", "+CLOCK+CAMP"]:
        r = s8row(ds, "C/K_G4", ct, "G4", adj)
        rowsG.append(dict(cohort=coh, dataset=ds, adjustment=adj, beta=r.beta, ci_lo=r.ci_lo, ci_hi=r.ci_hi, p=r.p, n=int(r.n), n_case=int(r.n_case), y=yy))
        yy += 1
    yy += 0.8
GG = pd.DataFrame(rowsG)
topG = 22.4
GB = 9.5
ag4 = ax_mm(x0 + 20.5, y0 + 6 + topG, 13.5, h - 6 - topG - GB)
agl = ax_mm(x0, y0 + 6 + topG, 20.0, h - 6 - topG - GB); agl.axis("off")
agt = ax_mm(x0 + 35.0, y0 + 6 + topG, w - 35.0, h - 6 - topG - GB); agt.axis("off")
for a in (ag4, agl, agt): a.set_ylim(GG.y.max() + 0.6, -1.5)
agl.set_xlim(0, 1); agt.set_xlim(0, 1)
zero_line(ag4)
for _, r in GG.iterrows():
    colr = COHORT[r.cohort]
    ag4.plot([r.ci_lo, r.ci_hi], [r.y, r.y], color=colr, lw=0.8)
    adj_mark(ag4, r.beta, r.y, r.cohort, r.adjustment, ms=4.0)
    agl.text(0.98, r.y, ADJ_LAB[r.adjustment], fontsize=FS_G, ha="right", va="center")
    agt.text(0.0, r.y, m(r.beta), fontsize=FS_G, ha="left", va="center")
    ptext(agt, 1.0, r.y, r.p, FS_G, ha="right")
    isT = r.adjustment in ("+CLOCK", "+CAMP")
    isTp = isT or (r.adjustment == "base" and r.cohort == "GSE290367")
    reg("G", m(r.beta), round(r.beta, 4), f"Table S8a (GC-induced efferocytic genes, {r.dataset}, {r.adjustment})", "DCM - control, GC-induced efferocytic genes (SD)", "T" if isT else "S",
        "yes" if isT else "new")
    reg("G", fmt_p(r.p), r.p, f"Table S8a (GC-induced efferocytic genes, {r.dataset}, {r.adjustment})", "nominal P (HC3)", "T" if isT else "S", "yes" if isT else "new")
for coh, ds, ct, lab in G4:
    ys = GG[GG.cohort == coh]
    nn = ys.iloc[0]
    agl.text(0.0, ys.y.min() - 0.85, lab.replace("\n", " ") + f" ({nn.n - nn.n_case}/{nn.n_case})", fontsize=FS_G, fontweight="bold",
             color=COHORT[coh], ha="left", va="center", bbox=dict(fc="white", ec="none", pad=0.4), zorder=5)
    reg("G", f"({nn.n - nn.n_case}/{nn.n_case})", f"{nn.n - nn.n_case}/{nn.n_case}", ds, "n control / DCM", "T", "yes")
ftext(x0 + 35.0, y0 + 26.0, "SD", fontsize=FS_G, fontweight="bold", ha="left", va="center")
ftext(x0 + w, y0 + 26.0, "P", fontsize=FS_G, fontweight="bold", ha="right", va="center")
ag4.set_xlim(-2.2, 0.5); ag4.set_xticks([-2, 0]); ag4.set_xticklabels(["\u22122", "0"])
ag4.set_yticks([]); ag4.spines["left"].set_visible(False)
ag4.set_xlabel("DCM \u2212 control score (SD)", fontsize=FS_T, labelpad=0.5)
fr = fig.axes[0]
# outcome badges for the two planned explanation tests (both met), stacked under the 2-line title
_hG = h - 6 - topG - GB
badge(agl, 0.03, 1.0 + (topG - 7.8) / _hG, "Circadian test: met", True, fs=FS_G)
badge(agl, 0.03, 1.0 + (topG - 14.4) / _hG, "cAMP test: met", True, fs=FS_G)
# legend-only cross-reference: Figure 7C MAGNet base estimate (8-gene-adjusted models, score standardized to all 166 NF controls)
v12 = s7b[(s7b.dataset == "MAGNet_excl_overlap") & (s7b.module == "G4") & (s7b.contrast == "DCM_vs_NF") & (s7b.adjustment == "age+sex")].iloc[0]
reg("G (legend)", m(v12.beta), round(v12.beta, 4), "Table S7b (MAGNet_excl_overlap, DCM_vs_NF, GC-induced efferocytic genes, age+sex)",
    "MAGNet base estimate shown in Figure 7C (cross-reference in legend)", "S", "new")
# legend: reason for the difference (reference controls used to standardize the score)
_mt = pd.read_csv(FI + "04_results/validation_v1.1/bulk/magnet_sample_meta.csv")
_n_nf_all = int((_mt.etiology == "NF").sum())
_n_nf_here = int((mag.case == 0).sum())
reg("G (legend)", str(_n_nf_here), _n_nf_here, "perdonor_scores_MAGNet_exclOverlap_DCM.csv", "non-failing controls used to standardize the score here", "N", "new")
reg("G (legend)", str(_n_nf_all), _n_nf_all, "04_results/validation_v1.1/bulk/magnet_sample_meta.csv", "all MAGNet non-failing controls (standardization reference in Figure 7C)", "N", "new")
SRC["4G"] = GG.drop(columns="y")

# ------------------------------------------------------------------ save
assert min_font_audit(fig, 10) == [], min_font_audit(fig, 10)
save_figure(fig, "Figure_4", tight=False)

pd.DataFrame(REG).to_csv("Figure_4_number_registry.csv", index=False)
with pd.ExcelWriter("Figure_4_source_data.xlsx") as xw:
    for k, d in SRC.items():
        d.to_excel(xw, sheet_name=f"Fig{k}", index=False)
