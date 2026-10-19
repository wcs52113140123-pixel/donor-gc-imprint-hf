"""pbtools v1.1 (2026-10-03). Imports v1.0 unchanged and overrides two functions:
B1 ols_hc3: when any leverage h>=0.999 (singleton blocking level), fall back to HC2 and flag (v1.0 returned NaN).
B6 competitive: random genes drawn WITHOUT replacement within each expression bin and excluding module genes.
Policy B2 (thresholds) is enforced by callers: per donor x cell type; stromal/endothelial >=50 nuclei & >=10k UMI; myeloid >=20 & >=5k (both datasets)."""
import numpy as np, pandas as pd
from scipy import stats
from pbtools import tmm_factors, lognorm, design, ref_z, module_score   # noqa: F401 (re-export, unchanged)
def ols_hc3(D, y, j=1, return_kind=False):
    X = np.asarray(D, float); y = np.asarray(y, float); XtXi = np.linalg.pinv(X.T @ X); b = XtXi @ X.T @ y; e = y - X @ b
    h = np.einsum("ij,jk,ik->i", X, XtXi, X); kind = "HC3"
    if (h >= 0.999).any(): w = e ** 2 / np.clip(1 - h, 1e-12, None); w[h >= 0.999] = 0.0; kind = "HC2_singleton"
    else: w = (e / (1 - h)) ** 2
    V = XtXi @ (X.T * w) @ X @ XtXi; se = np.sqrt(V[j, j]); p = 2 * stats.t.sf(abs(b[j] / se), len(y) - X.shape[1])
    return (b[j], se, p, kind) if return_kind else (b[j], se, p)
def competitive(Z, gnames, gl, expr, D, j=1, n=1000, rng=None, bins=20, is_ref=None):
    rng = rng or np.random.default_rng(1); gnames = np.asarray(gnames); ok = ~np.isnan(Z).any(0); inmod = np.isin(gnames, list(gl))
    b = pd.qcut(pd.Series(expr).rank(method="first"), bins, labels=False).values
    idx = np.flatnonzero(inmod & ok); pool = {k: np.flatnonzero((b == k) & ok & ~inmod) for k in range(bins)}
    need = pd.Series(b[idx]).value_counts().to_dict()
    def draw():
        out = []
        for k, m in need.items(): out.extend(rng.choice(pool[k], m, replace=False))
        return np.array(out)
    obs = Z[:, idx].mean(1); nul = np.stack([Z[:, draw()].mean(1) for _ in range(n)], 1)
    Y = np.column_stack([obs, nul]); ref = np.ones(len(Y), bool) if is_ref is None else np.asarray(is_ref, bool)
    Y = (Y - Y[ref].mean(0)) / Y[ref].std(0); B = np.linalg.lstsq(np.asarray(D, float), Y, rcond=None)[0][j]; b0, bn = B[0], B[1:]
    return dict(beta_obs=b0, null_mean=bn.mean(), null_sd=bn.std(), z_comp=(b0 - bn.mean()) / bn.std(),
                p_comp=(np.sum(np.abs(bn - bn.mean()) >= abs(b0 - bn.mean())) + 1) / (n + 1))
