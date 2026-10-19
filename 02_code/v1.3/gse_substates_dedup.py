#!/usr/bin/env python3
"""Track Substates (v1.1 validation) - GSE290367 adult LV myeloid L2 subsets.
H3: composition (subset fractions) vs state (within-subset module scores) for DCM/HCM/ICM vs ND; H2-within-subset ND ageing (TOST +/-0.1 SD/decade, HKSJ).
Kitagawa-type decomposition of pooled-myeloid effect into composition-only and state-only counterfactual pseudobulks.
Frozen modules: 01_config/modules/modules_v1.1.tsv (sha256 checked). Unit = donor. DCM/HCM pooled-myeloid effects = discovery (circular)."""
import os, sys, hashlib, json, numpy as np, pandas as pd
from scipy import stats
ROOT = os.environ["PROJ_ROOT"]  # [v12-wrap] was hard-coded "/home/weichangsheng/projects/cardiac_ageing_target_clearance"; os.chdir(ROOT)
sys.path.insert(0, f"{ROOT}/03_code/04_models"); import pbtools_v11 as pb
OUT = "04_results/validation_v1.3/substates_dedup"  ; os.makedirs(OUT, exist_ok=True)  # [v1.3 dedup] new output dir; v1.1 outputs untouched
DROP = set(pd.read_csv(os.environ["DUP_FILE"]).record_drop)  # [v1.3 dedup] 10 duplicated individuals, keep record with more nuclei
MODF = "01_config/modules/modules_v1.1.tsv"
assert hashlib.sha256(open(MODF, "rb").read()).hexdigest() == "3a2abdad2b609002752b4d4833d8219114fd305b2339c7ddf406da86448591be"
mods = pd.read_csv(MODF, sep="\t"); MODS = {k: sorted(set(v)) for k, v in mods.groupby("module")["gene"]}
MLIST = ["RES_EFF_core", "MHCII_core", "RES_EFF_lit", "ISG_IFN", "MAC_PAN_ID"]
PB = "02_data/processed/GSE290367_pseudobulk"
meta = pd.read_csv(f"{PB}/pseudobulk_L2_meta.csv"); genes = pd.read_csv(f"{PB}/genes.csv").gene.astype(str).values
CNT = np.load(f"{PB}/pseudobulk_L2_counts.npy", mmap_mode="r"); assert CNT.shape == (len(meta), len(genes)), CNT.shape
meta["row"] = np.arange(len(meta))
log = {"dedup_dropped_records": sorted(DROP), "L2_shape": list(CNT.shape), "sex_levels": meta.sex.value_counts(dropna=False).to_dict()}
DIS = ["DCM", "HCM", "ICM"]; MINC, MINU = 20, 5000; MIN_REF_STUDY = 3
A = meta[(meta.region == "LV") & (meta.age_status == "postnatal") & (meta.age >= 20) & meta.disease.isin(["ND"] + DIS) & ~meta.donor_id.isin(DROP)].copy()  # [v1.3 dedup]
A["age10"] = A.age / 10; A["female"] = (A.sex == "female").astype(float); A["isND"] = A.disease == "ND"
log["adult_LV_donors_by_disease"] = A.drop_duplicates("donor_id").disease.value_counts().to_dict()
MY = sorted(c for c in A.cell_type.unique() if c.startswith("Myeloid:")); log["myeloid_L2_labels"] = MY
log["myeloid_rows_by_label_all"] = A[A.cell_type.isin(MY)].cell_type.value_counts().to_dict()
DON = A.drop_duplicates("donor_id").set_index("donor_id")[["age", "age10", "female", "sex", "study", "disease", "isND"]]
tq = lambda df_, a: stats.t.ppf(1 - a / 2, df_)

def fit(df, y, term, covs, cats=("study",)):
    d = df.assign(_y=y).dropna(subset=["_y", term] + list(covs)); D = pb.design(d, [term] + list(covs), cats)
    b, se, p, kind = pb.ols_hc3(D.values, d._y.values, j=1, return_kind=True); dfree = len(d) - D.shape[1]
    return dict(beta=b, se=se, ci_lo=b - tq(dfree, .05) * se, ci_hi=b + tq(dfree, .05) * se, p=p, df_resid=dfree, se_kind=kind, n=len(d)), d

def tost(b, se, dfree, bound=0.1):
    pl = stats.t.sf((b + bound) / se, dfree); pu = stats.t.cdf((b - bound) / se, dfree)
    return dict(ci90_lo=b - tq(dfree, .10) * se, ci90_hi=b + tq(dfree, .10) * se, p_tost=max(pl, pu), equivalent=bool(max(pl, pu) < .05))

def hksj(df, y, min_n=8):
    """Within-study slope (y ~ age10 + female, HC3), DL random effects, HKSJ variance (untruncated), t_{k-1}."""
    rows = []
    for st, d in df.assign(_y=y).dropna(subset=["_y"]).groupby("study"):
        if len(d) < min_n or d.age.std() < 5: continue
        cov = ["female"] if d.female.std() > 0 else []
        D = pb.design(d, ["age10"] + cov); b, se, p = pb.ols_hc3(D.values, d._y.values, 1); rows.append((st, len(d), b, se))
    if len(rows) < 2: return dict(k=len(rows)), rows
    yb = np.array([r[2] for r in rows]); v = np.array([r[3] for r in rows]) ** 2; w = 1 / v; k = len(yb)
    mu_f = (w * yb).sum() / w.sum(); Q = (w * (yb - mu_f) ** 2).sum(); C = w.sum() - (w ** 2).sum() / w.sum(); tau2 = max(0, (Q - (k - 1)) / C)
    ws = 1 / (v + tau2); mu = (ws * yb).sum() / ws.sum(); q = (ws * (yb - mu) ** 2).sum() / (k - 1); se = np.sqrt(q / ws.sum())
    out = dict(k=k, n=int(sum(r[1] for r in rows)), beta=mu, se=se, ci_lo=mu - tq(k - 1, .05) * se, ci_hi=mu + tq(k - 1, .05) * se,
               p=2 * stats.t.sf(abs(mu / se), k - 1), tau2=tau2, I2=max(0, (Q - (k - 1)) / Q) if Q > 0 else 0.0)
    out.update({f"{kk}": vv for kk, vv in tost(mu, se, k - 1).items()}); return out, rows

def subset_rows(labels, thr=True):
    s = A[A.cell_type.isin(labels)]
    g = s.groupby("donor_id").agg(n_cells=("n_cells", "sum"), umi=("total_umi", "sum"), rows=("row", list))
    if thr: g = g[(g.n_cells >= MINC) & (g.umi >= MINU)]
    return g

def scores_for(g):
    M = np.stack([np.asarray(CNT[r]).sum(0) for r in g.rows]).astype(float)
    df = DON.loc[g.index].copy(); df["n_cells"] = g.n_cells.values
    nref = df[df.isND].groupby("study").size(); okst = nref[nref >= MIN_REF_STUDY].index; keep = df.study.isin(okst).values
    df, M = df[keep], M[keep]
    lc, gk = pb.lognorm(M, genes, "CPM"); Z = pb.ref_z(lc, df.study.values, df.isND.values)
    S, ng = {}, {}
    for m in MLIST: S[m], ng[m] = pb.module_score(Z, gk, MODS[m], df.isND.values)
    S["EFF_MHCII_index"] = S["RES_EFF_core"] - S["MHCII_core"]; ng["EFF_MHCII_index"] = ng["RES_EFF_core"] + ng["MHCII_core"]
    for m in S: df[m] = S[m]
    return df, ng

SUBSETS = {lab.split(":")[1]: [lab] for lab in MY}; SUBSETS["Myeloid_pooled"] = MY
EFF = []; SCORES = []
for sname, labs in SUBSETS.items():
    g = subset_rows(labs); df, ng = scores_for(g); df["subset"] = sname; SCORES.append(df.reset_index())
    for m in ["RES_EFF_core", "MHCII_core", "EFF_MHCII_index", "RES_EFF_lit", "ISG_IFN", "MAC_PAN_ID"]:
        for dz in DIS:
            cnt = df.groupby("study").apply(lambda d: pd.Series({"nd": int(d.isND.sum()), "dz": int((d.disease == dz).sum())}))
            sts = cnt[(cnt.nd >= 3) & (cnt.dz >= 3)].index
            d = df[df.study.isin(sts) & df.disease.isin(["ND", dz])].assign(isD=lambda x: (x.disease == dz).astype(float))
            if d.isD.sum() < 3: continue
            for adj, covs in [("primary", ["age10", "female"]), ("adj_MAC_PAN_ID", ["age10", "female", "MAC_PAN_ID"])]:
                if adj != "primary" and m == "MAC_PAN_ID": continue
                r, dd = fit(d, d[m].values, "isD", covs)
                EFF.append(dict(dataset="GSE290367", analysis="disease_state", subset=sname, module=m, contrast=f"{dz}_vs_ND", model=adj,
                                n_disease=int(dd.isD.sum()), n_ref=int((dd.isD == 0).sum()), n_studies=dd.study.nunique(), studies=";".join(sorted(dd.study.unique())),
                                n_genes=ng[m], **r))
        nd = df[df.isND]
        r, dd = fit(nd, nd[m].values, "age10", ["female"]); r.update(tost(r["beta"], r["se"], r["df_resid"]))
        EFF.append(dict(dataset="GSE290367", analysis="ND_ageing_state", subset=sname, module=m, contrast="per_decade_ND", model="pooled_study_FE",
                        n_disease=0, n_ref=r["n"], n_studies=dd.study.nunique(), studies=";".join(sorted(dd.study.unique())), n_genes=ng[m], **r))
        h, rows = hksj(nd, nd[m].values)
        if "beta" in h:
            EFF.append(dict(dataset="GSE290367", analysis="ND_ageing_state", subset=sname, module=m, contrast="per_decade_ND", model="HKSJ_within_study_meta",
                            n_disease=0, n_ref=h.pop("n"), n_studies=h["k"], studies=";".join(r_[0] for r_ in rows), n_genes=ng[m], **h))
SC = pd.concat(SCORES); SC.to_csv(f"{OUT}/GSE290367_subset_donor_scores.csv", index=False)

# ---------- composition: subset fractions within myeloid; myeloid share of all nuclei ----------
L1 = pd.read_csv(f"{PB}/pseudobulk_L1_meta.csv"); L1 = L1[(L1.region == "LV") & L1.donor_id.isin(DON.index)]
tot_all = L1.groupby("donor_id").n_cells.sum()
myc = A[A.cell_type.isin(MY)].pivot_table(index="donor_id", columns="cell_type", values="n_cells", aggfunc="sum", fill_value=0)
myc.columns = [c.split(":")[1] for c in myc.columns]; myc["TOTAL_MYELOID"] = myc.sum(1)
FR = DON.loc[myc.index].join(myc).join(tot_all.rename("TOTAL_NUCLEI")); FR = FR[FR.TOTAL_MYELOID >= 50]
nref = FR[FR.isND].groupby("study").size(); FR = FR[FR.study.isin(nref[nref >= 3].index)]
for k in [c.split(":")[1] for c in MY] + ["MYELOID_SHARE"]:
    num, den = (FR.TOTAL_MYELOID, FR.TOTAL_NUCLEI) if k == "MYELOID_SHARE" else (FR[k], FR.TOTAL_MYELOID)
    FR[f"logit_{k}"] = np.log((num + .5) / (den - num + .5)); FR[f"frac_{k}"] = num / den
FR.reset_index().to_csv(f"{OUT}/GSE290367_myeloid_fractions_by_donor.csv", index=False)
for k in [c.split(":")[1] for c in MY] + ["MYELOID_SHARE"]:
    y = FR[f"logit_{k}"]; sd = y[FR.isND].std()
    for dz in DIS:
        cnt = FR.groupby("study").apply(lambda d: pd.Series({"nd": int(d.isND.sum()), "dz": int((d.disease == dz).sum())}))
        sts = cnt[(cnt.nd >= 3) & (cnt.dz >= 3)].index
        d = FR[FR.study.isin(sts) & FR.disease.isin(["ND", dz])].assign(isD=lambda x: (x.disease == dz).astype(float))
        r, dd = fit(d, d[f"logit_{k}"].values, "isD", ["age10", "female"])
        EFF.append(dict(dataset="GSE290367", analysis="disease_composition", subset=k, module="logit_fraction", contrast=f"{dz}_vs_ND", model="primary",
                        n_disease=int(dd.isD.sum()), n_ref=int((dd.isD == 0).sum()), n_studies=dd.study.nunique(), studies=";".join(sorted(dd.study.unique())),
                        mean_frac_disease=dd.loc[dd.isD == 1, f"frac_{k}"].mean(), mean_frac_ref=dd.loc[dd.isD == 0, f"frac_{k}"].mean(), **r))
    nd = FR[FR.isND]; r, dd = fit(nd, (nd[f"logit_{k}"] / sd).values, "age10", ["female"]); r.update(tost(r["beta"], r["se"], r["df_resid"]))
    EFF.append(dict(dataset="GSE290367", analysis="ND_ageing_composition", subset=k, module="logit_fraction_in_ND_SD", contrast="per_decade_ND", model="pooled_study_FE",
                    n_disease=0, n_ref=r["n"], n_studies=dd.study.nunique(), studies=";".join(sorted(dd.study.unique())), **r))

# ---------- Kitagawa decomposition of pooled-myeloid scores ----------
labs = [c.split(":")[1] for c in MY]
allg = subset_rows(MY); donors = allg.index[DON.loc[allg.index].study.isin(nref[nref >= 3].index)]
X = {}; NC = {}; PASS = {}
for lab in MY:
    s = A[A.cell_type == lab].set_index("donor_id")
    for dn in donors:
        if dn in s.index:
            rr = s.loc[[dn]]; cts = np.asarray(CNT[rr.row.values]).sum(0).astype(float); n = rr.n_cells.sum(); u = rr.total_umi.sum()
            X[(lab, dn)] = cts / n; NC[(lab, dn)] = n; PASS[(lab, dn)] = (n >= MINC) & (u >= MINU)
        else: NC[(lab, dn)] = 0; PASS[(lab, dn)] = False
dd0 = DON.loc[donors]; F = pd.DataFrame({lab: [NC[(lab, d)] for d in donors] for lab in MY}, index=donors); F = F.div(F.sum(1), axis=0)
xbar = {}; fbar = {}
for st, idx in dd0.groupby("study").groups.items():
    ndd = [d for d in idx if DON.loc[d, "isND"]]
    fbar[st] = F.loc[ndd].mean()
    for lab in MY:
        pas = [X[(lab, d)] for d in ndd if PASS[(lab, d)]]
        if len(pas) == 0: pas = [X[(lab, d)] for d in DON.index[DON.isND] if (lab, d) in PASS and PASS[(lab, d)]]
        xbar[(lab, st)] = np.mean(pas, 0)
P_obs, P_comp, P_state = [], [], []; imput = {lab: 0 for lab in MY}
for d in donors:
    st = DON.loc[d, "study"]
    P_obs.append(sum(F.loc[d, lab] * X[(lab, d)] for lab in MY if NC[(lab, d)] > 0))
    P_comp.append(sum(F.loc[d, lab] * xbar[(lab, st)] for lab in MY))
    xs = []
    for lab in MY:
        if PASS[(lab, d)]: xs.append(fbar[st][lab] * X[(lab, d)])
        else: xs.append(fbar[st][lab] * xbar[(lab, st)]); imput[lab] += 1
    P_state.append(sum(xs))
log["decomp_donors"] = len(donors); log["decomp_imputed_subset_profiles"] = imput
Mobs = np.stack(P_obs); cpm = Mobs / Mobs.sum(1, keepdims=True) * 1e6; keepg = (cpm > 1).mean(0) >= .2
DEC = dd0.copy()
# COMMON SCALING: all counterfactuals use the OBSERVED pooled profile's per-study ND gene means, pooled ND gene SD and ND score mean/SD,
# so composition-only and state-only effects are in the same observed-SD units and approximately add to the observed effect.
isr = dd0.isND.values; stv = dd0.study.values
lc_obs = np.log2(np.stack(P_obs)[:, keepg] / np.stack(P_obs).sum(1, keepdims=True) * 1e6 + 1)
Lo = pd.DataFrame(lc_obs); MU = Lo[isr].groupby(stv[isr]).mean(); SDg = (Lo - MU.reindex(stv).values)[isr].std().replace(0, np.nan).values
SPAR = {}
for tag, Mx in [("observed", np.stack(P_obs)), ("composition_only", np.stack(P_comp)), ("state_only", np.stack(P_state))]:
    lc = np.log2(Mx[:, keepg] / Mx.sum(1, keepdims=True) * 1e6 + 1); Z = (lc - MU.reindex(stv).values) / SDg
    for m in ["RES_EFF_core", "MHCII_core", "MAC_PAN_ID"]:
        idx = [i for i, x in enumerate(genes[keepg]) if x in set(MODS[m])]; raw = np.nanmean(Z[:, idx], 1)
        if tag == "observed": SPAR[m] = (raw[isr].mean(), raw[isr].std())
        DEC[f"{tag}:{m}"] = (raw - SPAR[m][0]) / SPAR[m][1]
    DEC[f"{tag}:EFF_MHCII_index"] = DEC[f"{tag}:RES_EFF_core"] - DEC[f"{tag}:MHCII_core"]
DEC.reset_index().to_csv(f"{OUT}/GSE290367_decomposition_donor_scores.csv", index=False)
for tag in ["observed", "composition_only", "state_only"]:
    for m in ["RES_EFF_core", "MHCII_core", "EFF_MHCII_index"]:
        for dz in DIS:
            cnt = DEC.groupby("study").apply(lambda d: pd.Series({"nd": int(d.isND.sum()), "dz": int((d.disease == dz).sum())}))
            sts = cnt[(cnt.nd >= 3) & (cnt.dz >= 3)].index
            d = DEC[DEC.study.isin(sts) & DEC.disease.isin(["ND", dz])].assign(isD=lambda x: (x.disease == dz).astype(float))
            r, dd = fit(d, d[f"{tag}:{m}"].values, "isD", ["age10", "female"])
            EFF.append(dict(dataset="GSE290367", analysis=f"decomposition_{tag}", subset="Myeloid_pooled_reconstructed", module=m, contrast=f"{dz}_vs_ND", model="primary",
                            n_disease=int(dd.isD.sum()), n_ref=int((dd.isD == 0).sum()), n_studies=dd.study.nunique(), studies=";".join(sorted(dd.study.unique())), **r))
        nd = DEC[DEC.isND]; r, dd = fit(nd, nd[f"{tag}:{m}"].values, "age10", ["female"]); r.update(tost(r["beta"], r["se"], r["df_resid"]))
        EFF.append(dict(dataset="GSE290367", analysis=f"decomposition_{tag}", subset="Myeloid_pooled_reconstructed", module=m, contrast="per_decade_ND", model="pooled_study_FE",
                        n_disease=0, n_ref=r["n"], n_studies=dd.study.nunique(), studies=";".join(sorted(dd.study.unique())), **r))

E = pd.DataFrame(EFF)
def fam(r):
    core = r.module in ("RES_EFF_core", "MHCII_core") and r.subset in ("Resident", "Monocyte_derived") and r.model in ("primary", "pooled_study_FE")
    if r.analysis == "disease_state": return "H3_state_primary" if core else "H3_state_secondary_exploratory"
    if r.analysis == "ND_ageing_state": return "H2_subset_ageing_primary" if core else "H2_subset_ageing_secondary"
    if r.analysis == "disease_composition": return "H3_composition"
    return r.analysis
E["family"] = E.apply(fam, axis=1)
def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); q = p[o] * len(p) / np.arange(1, len(p) + 1); q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty_like(p); out[o] = np.minimum(q, 1); return out
E["q_BH_family"] = np.nan; E["q_TOST_BH_family"] = np.nan
for (f, mdl), idx in E.groupby(["family", "model"]).groups.items():
    E.loc[idx, "q_BH_family"] = bh(E.loc[idx, "p"])
    if "p_tost" in E and E.loc[idx, "p_tost"].notna().all(): E.loc[idx, "q_TOST_BH_family"] = bh(E.loc[idx, "p_tost"])
E["circularity_note"] = np.where((E.subset.isin(["Myeloid_pooled", "Myeloid_pooled_reconstructed"])) & E.contrast.isin(["DCM_vs_ND", "HCM_vs_ND"]) & E.module.isin(["RES_EFF_core", "MHCII_core", "EFF_MHCII_index"]),
                                 "DISCOVERY (circular: module selected on this contrast)",
                                 np.where(E.contrast.isin(["DCM_vs_ND", "HCM_vs_ND"]) & E.module.isin(["RES_EFF_core", "MHCII_core", "EFF_MHCII_index", "RES_EFF_lit"]), "within-subset decomposition of discovery contrast (new, not independent validation)", ""))
E.to_csv(f"{OUT}/substates_GSE290367_effects.csv", index=False)
json.dump(log, open(f"{OUT}/GSE290367_run_log.json", "w"), indent=1, default=str)
print(json.dumps(log, default=str)[:3000])
pd.set_option("display.width", 250)
print(E[E.module.isin(["RES_EFF_core", "MHCII_core", "EFF_MHCII_index", "logit_fraction"]) & (E.model.isin(["primary", "pooled_study_FE", "HKSJ_within_study_meta"]))][["analysis", "subset", "module", "contrast", "model", "n_disease", "n_ref", "beta", "ci_lo", "ci_hi", "p"]].round(3).to_string())
