# Figures for manuscript v0.4 (JACC: BTS), labels updated to v0.4 terminology. All values read from saved result tables; no numbers typed by hand.
import json, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
apply_figure_style(sizes=(8, 7, 6))
W = 7.0  # in, JACC double column
CLS = {"GC_induced": "#D55E00", "GC_repressed": "#0072B2", "GC_insensitive": "#999999", "ambiguous": "#DDDDDD", "not measured": "white"}
CLS_LAB = {"GC_induced": "GC-induced", "GC_repressed": "GC-repressed", "GC_insensitive": "GC-insensitive", "ambiguous": "Ambiguous", "not measured": "Not measured"}
COH = {"SCP3342": "#7B3294", "Hahn": "#C2A5CF", "GSE290367": "#1B9E77", "MAGNet": "#A6761D", "GTEx": "#525252"}
CT = ["Cardiomyocyte", "Fibroblast", "Endothelial", "Pericyte", "Myeloid"]
PL = lambda ax, l, **k: panel_letter(ax, l, case="upper", **k)
TB = lambda ax: mpl.transforms.blended_transform_factory(ax.transAxes, ax.transData)
TD = lambda ax: mpl.transforms.blended_transform_factory(ax.transData, ax.transAxes)

gcc = pd.read_csv("v12/steroid/gc_gene_classification.csv"); meta = pd.read_csv("v12/steroid/meta_all_genes.csv.gz").set_index("gene")
S13 = pd.read_csv("v13/res/v13_score_tests.csv"); COLD = pd.read_csv("v13/res/v13_collinearity_posthoc.csv"); AU = pd.read_csv("v13/res/v13_impact_audit.csv")
D13 = json.load(open("v13/res/v13_decisions.json")); A3d = pd.read_csv("v13/res/v13_A3_positive_controls_descriptive.csv"); A4 = pd.read_csv("v13/A4_literature_gene_GC_classification.csv")
g25 = pd.read_csv("v12/G2_G5_results.csv"); sub = pd.read_csv("v13/dedup/res/substates_GSE290367_effects.csv")
SRC = {}

def forest(ax, y, b, lo, hi, color, marker="o", fill=True, size=4.5, lw=1.0):
    ax.hlines(y, lo, hi, color=color, lw=lw, zorder=2)
    ax.plot(b, y, marker, ms=size, mfc=color if fill else "white", mec=color, mew=1.0, ls="none", zorder=3)

def clsmark(k):
    amb = k in ("ambiguous", "not measured"); return dict(mfc="white" if amb else CLS[k], mec="#8C8C8C" if amb else CLS[k]), ("#8C8C8C" if amb else CLS[k])

def fig1():
    fig = plt.figure(figsize=(W, 4.3)); gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 0.9], wspace=1.15)
    ax = fig.add_subplot(gs[0])
    order = [("GC_specific3", ["FKBP5", "TSC22D3", "ZBTB16"]),
             ("RES_EFF_core", ["CD163", "LYVE1", "C1QB", "C1QA", "C1QC", "MERTK", "MARCO", "F13A1", "LRP1"]),
             ("MHCII_core", ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DMA", "HLA-DMB", "CD74", "CIITA"])]
    rows = []; y = 0; yt = []; ylab = []; grp = []
    for g, genes in order:
        ys = []
        for gene in genes:
            if gene not in meta.index: continue
            r = meta.loc[gene]; rows.append(dict(group=g, gene=gene, y=y, **r[["k", "meta_log2FC", "ci_lo", "ci_hi", "P", "class"]].to_dict()))
            yt.append(y); ylab.append(gene); ys.append(y); y -= 1
        grp.append((g, np.mean(ys))); y -= 0.8
    d = pd.DataFrame(rows); SRC["Fig1A"] = d
    for _, r in d.iterrows():
        mk, lc = clsmark(r["class"]); ax.hlines(r.y, r.ci_lo, r.ci_hi, color=lc, lw=1.0); ax.plot(r.meta_log2FC, r.y, "o", ms=4.5, mew=1.0, **mk)
    ax.axvline(0, color="k", lw=0.6); [ax.axvline(v, color="#BBBBBB", lw=0.5, ls=":") for v in (-0.5, 0.5)]
    ax.set_yticks(yt); ax.set_yticklabels(ylab, style="italic"); ax.set_xlabel("Glucocorticoid effect, meta log$_2$FC (95% CI)")
    ax.set_ylim(min(yt) - 0.8, max(yt) + 0.8)
    for g, yc in grp: ax.text(1.03, yc, {"GC_specific3": "GC score\ngenes", "RES_EFF_core": "Resident\nefferocytic", "MHCII_core": "MHC\nclass II"}[g], transform=TB(ax), va="center", fontsize=7)
    ax.set_title("Resident efferocytic genes are GC-induced;\nmost MHC-II genes are GC-repressed", loc="left", fontsize=8)
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", label=CLS_LAB[k], **clsmark(k)[0]) for k in ["GC_induced", "GC_repressed", "GC_insensitive", "ambiguous"]], loc="lower right", fontsize=6, handletextpad=0.2)
    PL(ax, "a", dx=-0.25)
    ax2 = fig.add_subplot(gs[1]); mods = ["RES_EFF_core", "RES_EFF_lit", "MHCII_core", "ISG_IFN", "MAC_PAN_ID"]
    names = {"RES_EFF_core": "Resident efferocytic\n(core)", "RES_EFF_lit": "Efferocytosis\n(literature)", "MHCII_core": "MHC class II", "ISG_IFN": "Interferon", "MAC_PAN_ID": "Pan-macrophage\nidentity"}
    ct = gcc[gcc.module.isin(mods)].assign(cls=lambda x: x["class"].fillna("not measured")).groupby(["module", "cls"]).size().unstack(fill_value=0).reindex(mods).fillna(0)
    SRC["Fig1B"] = ct.reset_index(); left = np.zeros(len(mods))
    for k in ["GC_induced", "GC_repressed", "GC_insensitive", "ambiguous", "not measured"]:
        if k not in ct: continue
        v = ct[k].values; ax2.barh(range(len(mods)), v, left=left, color=CLS[k], edgecolor="#8C8C8C" if k in ("ambiguous", "not measured") else CLS[k], lw=0.5, label=CLS_LAB[k]); left += v
    ax2.set_yticks(range(len(mods))); ax2.set_yticklabels(names[m] for m in mods); ax2.invert_yaxis(); ax2.set_xlabel("Genes (n)")
    ax2.set_title("Classification by module", loc="left", fontsize=8); ax2.legend(loc="lower right", fontsize=6)
    PL(ax2, "b", dx=-0.75)
    return fig

def fig2():
    fig = plt.figure(figsize=(W, 5.9)); gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.75, wspace=0.5)
    ax = fig.add_subplot(gs[0, :]); base = S13[S13.test.isin(["C/K", "C/K_bulk"]) & (S13.adjustment == "base") & (S13.module == "GC3")].copy()
    base["eff"] = -base.beta; base["lo"] = -base.ci_hi; base["hi"] = -base.ci_lo; SRC["Fig2A"] = base
    x = 0; xt = []; xl = []
    for ds, coh, head in [("SCP3342_HFpEF", "SCP3342", "SCP3342 HFpEF (snRNA)"), ("GSE290367_DCM", "GSE290367", "GSE290367 DCM (snRNA)"), ("MAGNet_exclOverlap_DCM", "MAGNet", "MAGNet (bulk)")]:
        dd = base[base.dataset == ds]; cts = CT if coh != "MAGNet" else ["bulk"]
        for c_ in cts:
            r = dd[dd.cell_type == c_].iloc[0]; sig = r.p < 0.05
            ax.vlines(x, r.lo, r.hi, color=COH[coh], lw=1.1); ax.plot(x, r.eff, "o", ms=5, mfc=COH[coh] if sig else "white", mec=COH[coh])
            xt.append(x); xl.append(c_ if c_ != "bulk" else "LV bulk"); x += 1
        ax.text(x - (len(cts) + 1) / 2, 1.02, head, transform=TD(ax), ha="center", va="bottom", fontsize=7, color=COH[coh]); x += 0.8
    ax.axhline(0, color="k", lw=0.6); ax.set_xticks(xt); ax.set_xticklabels(xl, rotation=30, ha="right"); ax.set_xlim(-0.6, x - 0.2)
    ax.set_ylabel("GC score, control − case\n(SD, 95% CI)"); ax.set_title("Control hearts carry a higher GC score in two of three cohorts (open circle: P ≥ 0.05)", loc="left", fontsize=8, pad=22)
    PL(ax, "a", dx=-0.08, dy=1.12)
    ax2 = fig.add_subplot(gs[1, 0]); T = S13[S13.test.isin(["C/K", "C/K_bulk"]) & (S13.module == "GC3")].copy(); T["eff"] = -T.beta
    pv = T.pivot_table(index=["dataset", "cell_type"], columns="adjustment", values="eff").reset_index(); SRC["Fig2B"] = pv
    for adj, mk in [("+CLOCK", "o"), ("+CAMP", "s")]:
        for ds, coh in [("SCP3342_HFpEF", "SCP3342"), ("GSE290367_DCM", "GSE290367"), ("MAGNet_exclOverlap_DCM", "MAGNet")]:
            dd = pv[pv.dataset == ds]; ax2.plot(dd["base"], dd[adj], mk, ms=4.5, mfc=COH[coh] if adj == "+CLOCK" else "white", mec=COH[coh], ls="none")
    lim = [-1.5, 7]; ax2.plot(lim, lim, color="#999999", lw=0.6, ls="--"); ax2.plot([0, 7], [0, 3.5], color="#BBBBBB", lw=0.6, ls=":")
    ax2.text(6.9, 2.9, "50% retained", fontsize=6, color="#8C8C8C", ha="right"); ax2.text(3.2, 4.6, "fully retained", fontsize=6, color="#8C8C8C", ha="right")
    ax2.set_xlim(lim); ax2.set_ylim(lim); ax2.set_xlabel("GC score difference, base model (SD)"); ax2.set_ylabel("After adding covariate (SD)")
    ax2.legend(handles=[Line2D([], [], marker="o", ls="none", mfc="k", mec="k", label="+ circadian score"), Line2D([], [], marker="s", ls="none", mfc="white", mec="k", label="+ cAMP score")], loc="upper left", fontsize=6)
    ax2.set_title("Circadian and cAMP scores leave\nthe GC difference largely intact", loc="left", fontsize=8); PL(ax2, "b", dx=-0.25)
    ax3 = fig.add_subplot(gs[1, 1]); SRC["Fig2C"] = COLD
    for _, r in COLD.iterrows():
        coh = "SCP3342" if r.dataset.startswith("SCP") else "Hahn" if r.dataset.startswith("Hahn") else "GSE290367" if r.dataset.startswith("GSE") else "MAGNet"
        ax3.plot(r.AUC_GC3_ctrl_gt_case, r.VIF_case_coef, "o", ms=5, mfc=COH[coh], mec="k" if coh == "Hahn" else COH[coh], mew=0.5)
    ax3.axvline(0.9, color="#8C8C8C", lw=0.6, ls="--"); ax3.axhline(5, color="#8C8C8C", lw=0.6, ls="--"); ax3.axvspan(0.9, 1.04, color="#F2E6F7", zorder=0)
    ax3.set_xlim(0.1, 1.04); ax3.set_yscale("log"); ax3.set_yticks([1, 2, 5, 10]); ax3.set_yticklabels(["1", "2", "5", "10"]); ax3.minorticks_off()
    ax3.set_xlabel("AUC of GC score (control > case)"); ax3.set_ylabel("Variance inflation of\ndisease effect")
    ax3.text(0.965, 2.2, "not identifiable", fontsize=6, ha="center", va="center", rotation=90, color="#7B3294")
    ax3.legend(handles=[Line2D([], [], marker="o", ls="none", mfc=COH[k], mec="k" if k == "Hahn" else COH[k], mew=0.5, label=k) for k in ["SCP3342", "Hahn", "MAGNet", "GSE290367"]], loc="upper left", fontsize=6)
    ax3.set_title("GC score and diagnosis are aliased\nin the HFpEF program (post hoc)", loc="left", fontsize=8); PL(ax3, "c", dx=-0.3)
    return fig

def fig3():
    fig = plt.figure(figsize=(W, 4.0)); gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1], wspace=0.75)
    ax = fig.add_subplot(gs[0]); mods = [("G3", "GC-insensitive\nreceptors"), ("G4", "GC-induced\nsubset"), ("RES_EFF_core", "Resident\nefferocytic"), ("MHCII_core", "MHC class II")]
    spec = [("GSE290367_Resident", "DCM_vs_ND", "age+sex+study", "GSE290367 resident", COH["GSE290367"], True, "o"),
            ("GSE290367_Resident", "DCM_vs_ND", "age+sex+study+GC", "GSE290367 resident + 8-gene GC score", COH["GSE290367"], False, "o"),
            ("MAGNet_excl_overlap", "DCM_vs_NF", "age+sex", "MAGNet bulk", COH["MAGNet"], True, "s"),
            ("MAGNet_excl_overlap", "DCM_vs_NF", "age+sex+GC", "MAGNet bulk + 8-gene GC score", COH["MAGNet"], False, "s")]
    rows = []; yt = []
    for i, (m, lab) in enumerate(mods):
        for j, (ds, con, adj, sl, col, fill, mk) in enumerate(spec):
            r = g25[(g25.dataset == ds) & (g25.contrast == con) & (g25.module == m) & (g25.adjustment == adj)].iloc[0]
            y = -i * 5.2 - j; forest(ax, y, r.beta, r.ci_lo, r.ci_hi, col, marker=mk, fill=fill); rows.append(dict(module=m, series=sl, y=y, beta=r.beta, ci_lo=r.ci_lo, ci_hi=r.ci_hi, p=r.p, n=r.n))
        yt.append(-i * 5.2 - 1.5)
    SRC["Fig3A"] = pd.DataFrame(rows); ax.axvline(0, color="k", lw=0.6); ax.set_yticks(yt); ax.set_yticklabels([l for _, l in mods])
    ax.set_xlabel("DCM − control (SD, 95% CI)"); ax.set_title("Macrophage-state shifts shrink only partly\nafter GC adjustment", loc="left", fontsize=8)
    ax.legend(handles=[Line2D([], [], marker=s[6], ls="none", mfc=s[4] if s[5] else "white", mec=s[4], label=s[3]) for s in spec], loc="upper left", fontsize=6, bbox_to_anchor=(0.0, -0.17), ncol=2)
    PL(ax, "a", dx=-0.3)
    ax2 = fig.add_subplot(gs[1]); rows = []; y = 0; yt = []; yl = []
    for m, lab in [("GC_response", "8-gene GC\nresponse score"), ("RES_EFF_core", "Resident\nefferocytic"), ("MHCII_core", "MHC class II")]:
        for tis, mk in [("GTEx_LV", "o"), ("GTEx_AA", "D")]:
            r = g25[(g25.dataset == tis) & (g25.test == "G2") & (g25.contrast == "ventilator_vs_other") & (g25.module == m)].iloc[0]
            forest(ax2, y, r.beta, r.ci_lo, r.ci_hi, COH["GTEx"], marker=mk, fill=r.p < 0.05); rows.append(dict(module=m, tissue=tis, y=y, beta=r.beta, ci_lo=r.ci_lo, ci_hi=r.ci_hi, p=r.p, n=r.n)); y -= 1
        yt.append(y + 1.5); yl.append(lab); y -= 0.8
    SRC["Fig3B"] = pd.DataFrame(rows); ax2.axvline(0, color="k", lw=0.6); ax2.set_yticks(yt); ax2.set_yticklabels(yl)
    ax2.set_xlabel("Ventilator-managed − other deaths\n(SD, 95% CI)"); ax2.set_title("GTEx ventilator-managed deaths: macrophage\nscores differ, GC score does not", loc="left", fontsize=8)
    ax2.legend(handles=[Line2D([], [], marker="o", ls="none", color=COH["GTEx"], label="Left ventricle"), Line2D([], [], marker="D", ls="none", color=COH["GTEx"], label="Atrial appendage")], loc="upper left", fontsize=6, bbox_to_anchor=(0.0, -0.24), ncol=2)
    PL(ax2, "b", dx=-0.42)
    return fig

def fig4():
    fig = plt.figure(figsize=(W, 6.6)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.95], width_ratios=[1, 1.2], hspace=0.95, wspace=0.75)
    ax = fig.add_subplot(gs[0, 0])
    lab = {("SCP3342_HFpEF", "Myeloid"): "SCP3342 myeloid", ("SCP3342_HFpEF", "VascStromal"): "SCP3342 vasc.-stromal", ("SCP3342_HFpEF", "Fibroblast"): "SCP3342 fibroblast", ("SCP3342_HFpEF", "Endothelial"): "SCP3342 endothelial",
           ("Hahn2021_RVS_HFpEF", "bulk"): "Hahn HFpEF bulk", ("MAGNet_exclOverlap_DCM", "bulk"): "MAGNet DCM bulk", ("GSE290367_DCM", "Myeloid"): "GSE290367 myeloid", ("GSE290367_DCM", "Fibroblast"): "GSE290367 fibroblast"}
    cohc = lambda ds: COH["SCP3342"] if ds.startswith("SCP") else COH["Hahn"] if ds.startswith("Hahn") else COH["MAGNet"] if ds.startswith("MAG") else COH["GSE290367"]
    rows = []
    for i, ((ds, ct), l) in enumerate(lab.items()):
        b = AU[(AU.dataset == ds) & (AU.cell_type == ct) & (AU.model == "base")].iloc[0]; a = AU[(AU.dataset == ds) & (AU.cell_type == ct) & (AU.model == "GC3_adjusted")].iloc[0]
        ax.plot([b.top100down_frac_GCinduced * 100, a.top100down_frac_GCinduced * 100], [-i, -i], color="#BBBBBB", lw=0.8, zorder=1)
        ax.plot(b.top100down_frac_GCinduced * 100, -i, "o", ms=5, color=cohc(ds), zorder=3); ax.plot(a.top100down_frac_GCinduced * 100, -i, "o", ms=5, mfc="white", mec=cohc(ds), zorder=3)
        rows.append(dict(contrast=l, base=b.top100down_frac_GCinduced, adjusted=a.top100down_frac_GCinduced, background=b.bg_frac_GCinduced, n_fdr05_base=b.n_fdr05, n_fdr05_adj=a.n_fdr05, OR_down_base=b.OR_down_GCinduced))
    SRC["Fig4A"] = pd.DataFrame(rows); bg = SRC["Fig4A"].background.mean() * 100
    ax.axvline(bg, color="#8C8C8C", lw=0.6, ls="--"); ax.text(bg + 1.5, 0.3, "background", fontsize=6, color="#8C8C8C", va="center")
    ax.set_yticks(-np.arange(len(lab))); ax.set_yticklabels(list(lab.values())); ax.set_xlabel("GC-induced genes among top 100\ndisease-down-regulated genes (%)")
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", color="k", label="standard model"), Line2D([], [], marker="o", ls="none", mfc="white", mec="k", label="+ GC score")], loc="lower right", fontsize=6)
    ax.set_title("GC targets are over-represented\namong top disease genes", loc="left", fontsize=8); ax.set_ylim(-len(lab) + 0.4, 0.6); PL(ax, "a", dx=-0.75)
    ax2 = fig.add_subplot(gs[0, 1]); genes = ["NPPA", "NPPB", "MYH6", "POSTN", "ASPN", "SFRP4", "FRZB"]
    mag = pd.DataFrame(D13["A3"]["table"]).assign(dataset="MAGNet_exclOverlap_DCM", cell_type="bulk")
    allA = pd.concat([mag[["gene", "model", "beta", "fdr", "dataset", "cell_type"]], A3d[["gene", "model", "beta", "fdr", "dataset", "cell_type"]]])
    cols = [("MAGNet_exclOverlap_DCM", "bulk", "MAGNet\nbulk"), ("GSE290367_DCM", "Fibroblast", "GSE\nfibro."), ("GSE290367_DCM", "Cardiomyocyte", "GSE\ncardio."), ("SCP3342_HFpEF", "Fibroblast", "SCP\nfibro."), ("Hahn2021_RVS_HFpEF", "bulk", "Hahn\nbulk")]
    M = np.full((len(genes), 2 * len(cols)), np.nan); S = np.zeros_like(M, bool)
    for j, (ds, ct, _) in enumerate(cols):
        for k, md in enumerate(["base", "GC3_adjusted"]):
            for i, g in enumerate(genes):
                x = allA[(allA.dataset == ds) & (allA.cell_type == ct) & (allA.model == md) & (allA.gene == g)]
                if len(x): M[i, 2 * j + k] = x.beta.iloc[0]; S[i, 2 * j + k] = x.fdr.iloc[0] < 0.05
    SRC["Fig4B"] = allA; im = ax2.imshow(M, cmap="RdBu_r", vmin=-3.5, vmax=3.5, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isnan(M[i, j]): ax2.text(j, i, "n.d.", ha="center", va="center", fontsize=6, color="#666666")
            elif S[i, j]: ax2.text(j, i + 0.15, "*", ha="center", va="center", fontsize=8, color="white" if abs(M[i, j]) > 2.3 else "k")
    ax2.set_yticks(range(len(genes))); ax2.set_yticklabels(genes, style="italic"); ax2.set_xticks(range(M.shape[1])); ax2.set_xticklabels(["base", "+GC"] * len(cols), rotation=90, fontsize=6)
    for j, (_, _, l) in enumerate(cols): ax2.text(2 * j + 0.5, -0.8, l, ha="center", va="bottom", fontsize=6)
    [ax2.axvline(2 * j - 0.5, color="white", lw=2) for j in range(1, len(cols))]
    cb = fig.colorbar(im, ax=ax2, fraction=0.04, pad=0.03); cb.set_label("Disease effect (SD)\n* FDR < 0.05", fontsize=6); cb.ax.tick_params(labelsize=6)
    ax2.set_title("Known HF genes survive adjustment only\nwhere the design is identifiable", loc="left", fontsize=8, pad=30); PL(ax2, "b", dx=-0.25, dy=1.3)
    ax3 = fig.add_subplot(gs[1, :]); a4 = A4.copy()
    a4["src"] = a4.source.map(lambda s: "JACC BTS 2026 (LVAD)" if s.startswith("JACC") else "Amrute 2022 (recovery)" if s.startswith("Amrute") else "Nat Commun\n2025" if s.startswith("Nat Commun") else "Circ Res\n2026 (HFpEF)")
    a4 = a4.dropna(subset=["meta_log2FC"]); SRC["Fig4C"] = a4; x = 0; xt = []; xl = []
    for s_ in ["Circ Res\n2026 (HFpEF)", "Nat Commun\n2025", "JACC BTS 2026 (LVAD)", "Amrute 2022 (recovery)"]:
        dd = a4[a4.src == s_].sort_values("meta_log2FC", ascending=False); x0 = x
        for _, r in dd.iterrows():
            lo, hi = [float(v) for v in r.ci.strip("[]").split(",")]; mk, lc = clsmark(r.GC_class)
            ax3.vlines(x, lo, hi, color=lc, lw=1.0); ax3.plot(x, r.meta_log2FC, "o", ms=5, mew=1.0, **mk); xt.append(x); xl.append(r.gene); x += 1
        ax3.text((x0 + x - 1) / 2, 1.03, s_, transform=TD(ax3), ha="center", va="bottom", fontsize=6); x += 1.6
    ax3.axhline(0, color="k", lw=0.6); ax3.set_xticks(xt); ax3.set_xticklabels(xl, rotation=45, ha="right", style="italic"); ax3.set_xlim(-0.7, x - 0.2)
    ax3.set_ylabel("GC effect in vitro\n(meta log$_2$FC, 95% CI)")
    ax3.legend(handles=[Line2D([], [], marker="o", ls="none", label=CLS_LAB[k], **clsmark(k)[0]) for k in ["GC_induced", "GC_repressed", "ambiguous"]], loc="upper right", fontsize=6)
    ax3.set_title("Most genes read as heart-failure or recovery biology in recent studies are GC-induced in vitro", loc="left", fontsize=8, pad=26); PL(ax3, "c", dx=-0.1, dy=1.22)
    return fig

def fig5():
    fig, ax = plt.subplots(figsize=(W * 0.66, 3.9)); rows = []; y = 0; yt = []; yl = []; centers = []
    spec = [("GTEx_LV", "age_per_decade_non_ventilator_DTHHRDY1-4", "GTEx LV, not ventilator-managed"), ("GTEx_AA", "age_per_decade_non_ventilator_DTHHRDY1-4", "GTEx AA, not ventilator-managed"),
            ("GTEx_LV", "age_per_decade_ventilator_only", "GTEx LV, ventilator-managed"), ("GTEx_AA", "age_per_decade_ventilator_only", "GTEx AA, ventilator-managed")]
    for m, mk in [("RES_EFF_core", "o"), ("MHCII_core", "s")]:
        y0 = y
        for ds, con, l in spec:
            r = g25[(g25.dataset == ds) & (g25.test == "G5") & (g25.contrast == con) & (g25.module == m)].iloc[0]
            forest(ax, y, r.beta, r.ci_lo, r.ci_hi, COH["GTEx"], marker=mk, fill=r.p_tost < 0.05); rows.append(dict(module=m, series=l, y=y, beta=r.beta, ci_lo=r.ci_lo, ci_hi=r.ci_hi, p=r.p, p_tost=r.p_tost, n=r.n)); yt.append(y); yl.append(l); y -= 1
        for subset, l in [("Resident", "GSE290367 resident"), ("Monocyte_derived", "GSE290367 monocyte-derived")]:
            r = sub[(sub.analysis == "ND_ageing_state") & (sub.subset == subset) & (sub.module == m) & (sub.contrast == "per_decade_ND") & (sub.model == "pooled_study_FE")].iloc[0]
            forest(ax, y, r.beta, r.ci_lo, r.ci_hi, COH["GSE290367"], marker=mk, fill=bool(r.p_tost < 0.05)); rows.append(dict(module=m, series=l, y=y, beta=r.beta, ci_lo=r.ci_lo, ci_hi=r.ci_hi, p=r.p, p_tost=r.p_tost, n=r.n)); yt.append(y); yl.append(l); y -= 1
        centers.append((y0 + y + 1) / 2); y -= 0.8
    SRC["Fig5"] = pd.DataFrame(rows); ax.axvspan(-0.1, 0.1, color="#EEEEEE", zorder=0); ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(yt); ax.set_yticklabels(yl); ax.set_xlabel("Change per decade of age (SD, 95% CI)")
    for c_, l in zip(centers, ["Resident\nefferocytic", "MHC class II"]): ax.text(1.03, c_, l, transform=TB(ax), va="center", fontsize=7)
    ax.text(0, 0.75, "±0.1 SD equivalence", ha="center", va="bottom", fontsize=6, color="#8C8C8C", bbox=dict(fc="#EEEEEE", ec="none", pad=0.5))
    ax.set_title("Efferocytic genes: equivalent in 3 of 4 GTEx strata, small\ndecline in GSE290367 donors (filled: equivalence P < 0.05)", loc="left", fontsize=8)
    ax.set_ylim(y + 0.5, 1.8)
    return fig

FIGS = {"Fig1": fig1, "Fig2": fig2, "Fig3": fig3, "Fig4": fig4, "Fig5": fig5}
