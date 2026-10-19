#!/usr/bin/env python3
"""v1.2 pre-registered test G1 (ANALYSIS_PLAN_v1.2.md, sha256 26db5075...02b2): is the control-vs-case GC_response difference
systemic (all major cell types) or myeloid-restricted?
Implementation notes (fixed before running):
- GC_response = FKBP5 TSC22D3 ZBTB16 PER1 KLF9 DDIT4 KLF15 MT2A; score = mean reference-standardised z (pbtools.ref_z/module_score).
- Thresholds per donor x cell type: Myeloid >=20 nuclei & >=5k UMI; others >=50 & >=10k.
- GSE290367: adult LV postnatal, de-duplicated (drop record_drop), DCM vs ND, studies with >=1 DCM and >=3 ND donors for that cell type;
  model score ~ case + age/10 + female + study; reference = ND, centred within study. Statistic = control - case = -beta_case.
- SCP3342: HFpEF vs control; model score ~ case + female + pool; reference = controls (overall centring). Statistic = -beta_case.
- HC3 via pbtools_v11.ols_hc3 (HC2 on singleton blocks). Holm across the 5 cell types within each dataset.
- Decision (plan): PASS if in BOTH SCP3342 and GSE290367-DCM: >=4/5 cell types with control>case AND >=3/5 Holm P<0.05
  (Holm-significant cell types counted only when direction is control>case — clarification written before running).
- Pre-specified sub-analysis: Chaffin 2022 only (DCM). Descriptive: HCM vs ND. Auxiliary: median pairwise Spearman rho of
  GC_response across cell types among reference donors (donors present in both cell types)."""
import sys, os, json, itertools, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, "audit/03_code/04_models"); import pbtools_v11 as pb
GC = ["FKBP5", "TSC22D3", "ZBTB16", "PER1", "KLF9", "DDIT4", "KLF15", "MT2A"]
THR = lambda ct: (20, 5000) if ct == "Myeloid" else (50, 10000)
OUT = "v12/out_G1"; os.makedirs(OUT, exist_ok=True)
rows, refscores = [], {}

def test(d, M, genes, case_col, groups, ref, terms, cats):
    lc, gk = pb.lognorm(M, genes, "CPM"); Z = pb.ref_z(lc, groups, ref); ok = ~np.isnan(Z).any(0)
    y, ng = pb.module_score(Z[:, ok], gk[ok], GC, ref)
    D = pb.design(d, [case_col] + terms, cats); b, se, p, kind = pb.ols_hc3(D.values, y, 1, return_kind=True)
    t = stats.t.ppf(.975, len(y) - D.shape[1])
    return dict(n_case=int(d[case_col].sum()), n_ref=int(ref.sum()), n_GC_genes=ng, ctrl_minus_case=-b, ci_lo=-b - t * se, ci_hi=-b + t * se, p=p, se_kind=kind), pd.Series(y, index=d.index)

# ---------------- GSE290367 ----------------
genes = pd.read_csv("pb/G/genes.csv").gene.astype(str).values
L1 = pd.read_csv("pb/G/pseudobulk_L1_meta.csv"); L1["age"] = pd.to_numeric(L1.age, errors="coerce"); L1["row"] = np.arange(len(L1))
X1 = np.load("pb/G/pseudobulk_L1_counts.npy", mmap_mode="r")
drop = set(pd.read_csv("gse290367_duplicate_individuals.csv").record_drop)
base = L1[(L1.region == "LV") & (L1.age_status == "postnatal") & (L1.age >= 20) & ~L1.donor_id.isin(drop)]
CTS = ["Cardiomyocyte", "Fibroblast", "Endothelial", "Pericyte", "Myeloid"]
def gse(dis, ct, studies=None):
    s = base[(base.cell_type == ct) & base.disease.isin(["ND", dis])]; thr = THR(ct)
    s = s[(s.n_cells >= thr[0]) & (s.total_umi >= thr[1])].groupby("donor_id").first()
    if studies is not None: s = s[s.study.isin(studies)]
    cnt = s.groupby("study").disease.agg(lambda x: ((x == dis).sum(), (x == "ND").sum()))
    keep = [st for st, (a, b) in cnt.items() if a >= 1 and b >= 3]; s = s[s.study.isin(keep)]
    d = s[["age", "sex", "study", "disease"]].copy(); d["case"] = (d.disease == dis).astype(float); d["female"] = (d.sex == "female").astype(float); d["age10"] = d.age / 10
    M = np.stack([np.asarray(X1[r]).astype(float) for r in s.row.values])
    r, y = test(d, M, genes, "case", d.study.values, (d.disease == "ND").values, ["age10", "female"], ["study"] if d.study.nunique() > 1 else [])
    r["studies"] = ";".join(sorted(keep)); return r, y, d
for label, dis, studies in [("GSE290367_DCM", "DCM", None), ("GSE290367_DCM_Chaffin", "DCM", ["Chaffin 2022"]), ("GSE290367_HCM_descriptive", "HCM", None)]:
    for ct in CTS:
        r, y, d = gse(dis, ct, studies); rows.append(dict(dataset=label, cell_type=ct, **r))
        if label == "GSE290367_DCM": refscores[("GSE290367", ct)] = y[(d.disease == "ND").values]
# ---------------- SCP3342 ----------------
sg = pd.read_csv("pb/S/genes.csv").gene.astype(str).values; sm = pd.read_csv("pb/S/pseudobulk_meta.csv", dtype={"donor_id": str})
C = np.load("pb/S/pseudobulk_raw_count_cellbender.npy", mmap_mode="r")
SMAP = {"Cardiomyocyte": ["Cardiomyocyte"], "Fibroblast": ["Fibroblast"], "Endothelial": ["Endothelial1", "Endothelial2"], "Pericyte": ["Pericyte"],
        "Myeloid": ["Macrophage", "Proliferating_macrophage"]}
for ct, labs in SMAP.items():
    ss = sm[sm.cell_type.isin(labs)].copy(); ss["row"] = ss.index
    g = ss.groupby("donor_id").agg(n=("n_cells", "sum"), umi=("total_umi_raw_count_cellbender", "sum"), disease=("disease", "first"), sex=("sex", "first"), pool=("pool", "first"))
    thr = THR(ct); g = g[(g.n >= thr[0]) & (g.umi >= thr[1])]
    M = np.stack([np.asarray(C[ss[ss.donor_id == dn].row.values]).astype(float).sum(0) for dn in g.index])
    g["case"] = (g.disease == "hfpef").astype(float); g["female"] = (g.sex == "female").astype(float)
    r, y = test(g, M, sg, "case", np.zeros(len(g)), (g.case == 0).values, ["female"], ["pool"])
    r["studies"] = "SCP3342"; rows.append(dict(dataset="SCP3342_HFpEF", cell_type=ct, **r)); refscores[("SCP3342", ct)] = y[(g.case == 0).values]
R = pd.DataFrame(rows)
R["p_holm"] = R.groupby("dataset").p.transform(lambda p: multipletests(p, method="holm")[1])
dec = {}
for ds in ["SCP3342_HFpEF", "GSE290367_DCM", "GSE290367_DCM_Chaffin"]:
    s = R[R.dataset == ds]; npos = int((s.ctrl_minus_case > 0).sum()); nsig = int(((s.p_holm < 0.05) & (s.ctrl_minus_case > 0)).sum())
    dec[ds] = dict(n_positive=npos, n_holm_sig_positive=nsig, n_holm_sig_any=int((s.p_holm < 0.05).sum()), criterion_met=bool(npos >= 4 and nsig >= 3))
dec["G1_PASS"] = bool(dec["SCP3342_HFpEF"]["criterion_met"] and dec["GSE290367_DCM"]["criterion_met"])
aux = []
for dsn in ["GSE290367", "SCP3342"]:
    for a, b in itertools.combinations(CTS, 2):
        x, y = refscores[(dsn, a)], refscores[(dsn, b)]; j = x.index.intersection(y.index)
        if len(j) >= 8: rho, p = stats.spearmanr(x[j], y[j]); aux.append(dict(dataset=dsn, ct_a=a, ct_b=b, n=len(j), rho=rho, p=p))
A = pd.DataFrame(aux); dec["aux_median_rho"] = A.groupby("dataset").rho.median().round(3).to_dict()
R.to_csv(f"{OUT}/G1_results.csv", index=False); A.to_csv(f"{OUT}/G1_aux_crosscelltype_rho.csv", index=False)
json.dump(dec, open(f"{OUT}/G1_decision.json", "w"), indent=1)
pd.set_option("display.width", 250); print(R.drop(columns=["studies"]).round(4).to_string()); print(json.dumps(dec, indent=1))
