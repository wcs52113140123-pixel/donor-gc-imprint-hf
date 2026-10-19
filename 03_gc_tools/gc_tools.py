"""gc_tools — flag and adjust a glucocorticoid (GC)-like transcriptional imprint in human cardiac (and other tissue) omics.

Reference: in-vitro GC response of human primary monocytes / monocyte-derived macrophages (random-effects meta of
GSE85333, GSE112101, E-MEXP-2598, E-MTAB-7572; QC: FKBP5 and TSC22D3 log2FC>1, FDR<0.05), see gc_gene_classification_all.csv.
Classes: GC_induced (meta log2FC>=0.5, P<0.01, >=2/3 datasets concordant), GC_repressed (symmetric), GC_insensitive
(95% CI within (-0.5, 0.5)), ambiguous (otherwise). Pre-registered in ANALYSIS_PLAN_v1.2/v1.3 (sha256 26db5075.../9e7ab16c...).

Functions
---------
gc_score(Z, genes, ref_mask, which="GC_specific3") -> per-sample score (reference-standardised mean z)
gc_flag(df, gene_col="gene", table=None) -> df with columns GC_class, GC_meta_log2FC
gc_adjusted_de(Z, genes, design, score) -> gene-wise OLS with the GC score added as covariate (BH-FDR)
Use: (1) always run gc_flag on published/obtained DE lists; (2) report DE with and without gc_adjusted_de;
(3) run gc_collinearity first — adjusted DE is only interpretable when the GC score is not aliased with diagnosis;
(4) treat genes flagged GC_induced/GC_repressed whose effect disappears after adjustment as not interpretable as disease biology.
"""
import numpy as np, pandas as pd
from scipy import stats
GC_SPECIFIC3 = ["FKBP5", "TSC22D3", "ZBTB16"]
GC_RESPONSE8 = ["FKBP5", "TSC22D3", "ZBTB16", "PER1", "KLF9", "DDIT4", "KLF15", "MT2A"]

def load_table(path=None):
    import os
    return pd.read_csv(path or os.path.join(os.path.dirname(os.path.abspath(__file__)), "gc_gene_classification_all.csv"))

def gc_score(Z, genes, ref_mask, which="GC_specific3"):
    """Z: samples x genes matrix of reference-standardised expression (e.g. (logCPM - mean_ref)/sd_ref). Returns z-scored mean."""
    gl = GC_SPECIFIC3 if which == "GC_specific3" else GC_RESPONSE8
    genes = np.asarray(genes); idx = [i for i, g in enumerate(genes) if g in set(gl)]
    if len(idx) < 2: raise ValueError(f"only {len(idx)} of {gl} measured")
    s = np.nanmean(np.asarray(Z, float)[:, idx], 1); ref = np.asarray(ref_mask, bool)
    return (s - s[ref].mean()) / s[ref].std(ddof=1)

def gc_flag(df, gene_col="gene", table=None):
    t = load_table() if table is None else table
    t = t.set_index("gene")[["class", "meta_log2FC"]].rename(columns={"class": "GC_class", "meta_log2FC": "GC_meta_log2FC"})
    out = df.merge(t, left_on=gene_col, right_index=True, how="left"); out["GC_class"] = out.GC_class.fillna("not_in_reference"); return out

def gc_collinearity(design, score, case_col=1):
    """Identifiability check before adjustment. Returns AUC (P(score_ctrl > score_case)), R2 of case on other
    covariates + score, and the variance inflation of the case coefficient. Post-hoc rule (v1.3 report, not pre-registered):
    if AUC >= 0.9 or <= 0.1, or VIF >= 5, do NOT interpret gc_adjusted_de; report the design as non-identifiable."""
    D0 = np.asarray(design, float); D1 = np.column_stack([D0, np.asarray(score, float)]); cs = D0[:, case_col].astype(bool)
    vif = np.linalg.pinv(D1.T @ D1)[case_col, case_col] / np.linalg.pinv(D0.T @ D0)[case_col, case_col]
    u = stats.mannwhitneyu(np.asarray(score)[~cs], np.asarray(score)[cs]); auc = u.statistic / ((~cs).sum() * cs.sum())
    oth = np.delete(D1, case_col, 1); f = oth @ np.linalg.lstsq(oth, D0[:, case_col], rcond=None)[0]; y = D0[:, case_col]
    r2 = 1 - ((y - f) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return dict(AUC=auc, R2_case=r2, VIF=vif, identifiable=bool((0.1 < auc < 0.9) and vif < 5))

def gc_adjusted_de(Z, genes, design, score, coef=1, check=True):
    """design: samples x p matrix with intercept in column 0 and the effect of interest in column `coef`.
    Returns DataFrame gene, beta, t, p, fdr for the model design + GC score."""
    from statsmodels.stats.multitest import multipletests
    if check:
        cc = gc_collinearity(design, score, coef)
        if not cc["identifiable"]: import warnings; warnings.warn(f"GC score nearly aliased with the contrast: {cc}; adjusted DE not interpretable")
    D = np.column_stack([np.asarray(design, float), np.asarray(score, float)]); Z = np.asarray(Z, float)
    B = np.linalg.lstsq(D, Z, rcond=None)[0]; R = Z - D @ B; df = len(Z) - D.shape[1]
    se = np.sqrt((R ** 2).sum(0) / df * np.linalg.pinv(D.T @ D)[coef, coef]); t = B[coef] / se; p = 2 * stats.t.sf(np.abs(t), df)
    out = pd.DataFrame({"gene": genes, "beta": B[coef], "t": t, "p": p}); out["fdr"] = multipletests(out.p, method="fdr_bh")[1]; return out
