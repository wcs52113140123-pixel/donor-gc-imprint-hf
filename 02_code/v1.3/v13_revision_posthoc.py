#!/usr/bin/env python3
"""v1.3 pre-registered tests (ANALYSIS_PLAN_v1.3.md, sha256 9e7ab16cbef7532bd5d3479a9fa961b0c4928c083e69a6f78010efb9e22929d6):
C (circadian) and K (cAMP/catecholamine) decision tests, A (impact audit, descriptive + A3 positive-control decision), inputs for T.
Gene sets are fixed in the plan (selected on in-vitro GC reference only). Gene-wise models: OLS on reference-standardised gene z
(classical SE, vectorised), BH-FDR. Score-level models: HC3 (pbtools_v11, HC2 on singleton blocks)."""
import os, sys, io, gzip, tarfile, json, hashlib, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
ROOT = "/home/weichangsheng/projects/cardiac_ageing_target_clearance"; sys.path.insert(0, f"{ROOT}/03_code/04_models"); import pbtools_v11 as pb
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True); INP = os.environ.get("INP", ".")
assert hashlib.sha256(open(f"{ROOT}/01_config/ANALYSIS_PLAN_v1.3.md", "rb").read()).hexdigest() == "9e7ab16cbef7532bd5d3479a9fa961b0c4928c083e69a6f78010efb9e22929d6"
GC3 = ["FKBP5", "TSC22D3", "ZBTB16"]; GCR = ["FKBP5", "TSC22D3", "ZBTB16", "PER1", "KLF9", "DDIT4", "KLF15", "MT2A"]
MORN = ["NR1D2", "DBP", "TEF", "PER2", "BHLHE41", "CIART"]; EVE = ["ARNTL", "NPAS2", "CLOCK", "CRY1"]
CAMP = ["NR4A1", "NR4A2", "CREM", "FOSB", "JUNB", "PDE4B", "PDE4D"]
POS = {"NPPA": 1, "NPPB": 1, "POSTN": 1, "ASPN": 1, "SFRP4": 1, "FRZB": 1, "MYH6": -1}
cls = pd.read_csv(f"{INP}/gc_gene_classification.csv"); lit = cls[cls.module == "RES_EFF_lit"]; G4 = sorted(lit.gene[lit["class"] == "GC_induced"])
mt = pd.read_csv(f"{INP}/meta_all_genes.csv.gz"); IND = set(mt.gene[mt["class"] == "GC_induced"]); REP = set(mt.gene[mt["class"] == "GC_repressed"]); MEAS = set(mt.gene)
ROWS, DEC, AUD, GENEW, COVG = [], {}, [], [], {}
RES_L = ['C1QA', 'C1QB', 'C1QC', 'CD163', 'F13A1', 'LRP1', 'LYVE1', 'MARCO', 'MERTK']; MHC_L = ['CD74', 'CIITA', 'HLA-DMA', 'HLA-DMB', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB1']; MAC_L = ['AIF1', 'CD14', 'CD68', 'CYBB', 'FCGR1A', 'MS4A6A', 'SPI1', 'TYROBP']; G3_L = ['ITGAV', 'ITGB5', 'LRP1', 'SCARB1']
REV = []; NULL = []; LAST = {}
def scores(lc, gk, groups, ref):
    LAST['mean'] = dict(zip(gk, np.asarray(lc).mean(0)))
    Z = pb.ref_z(lc, groups, ref); ok = ~np.isnan(Z).any(0); Z, gk = Z[:, ok], gk[ok]
    S = {}
    for k, v in {"GC3": GC3, "GCR": GCR, "MORN": MORN, "EVE": EVE, "CAMP": CAMP, "G4": G4, "RES": RES_L, "MHC": MHC_L, "MAC": MAC_L, "G3": G3_L}.items():
        S[k], COVG[k] = pb.module_score(Z, gk, v, ref)
    S = pd.DataFrame(S); S["CLOCK"] = S.MORN - S.EVE; S["CLOCK"] = (S.CLOCK - S.CLOCK[ref].mean()) / S.CLOCK[ref].std(); return S, Z, gk
def fit(ds, test, ct, module, adj, d, y, terms, cats, extra=None):
    D = pb.design(d, terms, cats); b, se, p, kind = pb.ols_hc3(D.values, y, 1, return_kind=True); t = stats.t.ppf(.975, len(y) - D.shape[1])
    r = dict(dataset=ds, test=test, cell_type=ct, module=module, adjustment=adj, n=len(y), n_case=int(d[terms[0]].sum()), beta=b, ci_lo=b - t * se, ci_hi=b + t * se, p=p, se_kind=kind)
    if extra: r.update(extra)
    ROWS.append(r); return r
def genewise(Z, gk, D):
    D = np.asarray(D, float); B = np.linalg.lstsq(D, Z, rcond=None)[0]; R = Z - D @ B; df = len(Z) - D.shape[1]
    se = np.sqrt((R ** 2).sum(0) / df * np.linalg.pinv(D.T @ D)[1, 1]); t = B[1] / se; p = 2 * stats.t.sf(np.abs(t), df)
    return pd.DataFrame({"gene": gk, "beta": B[1], "t": t, "p": p})
def audit(ds, ct, d, Z, gk, terms, cats, sc):
    out = {}
    for lab, tt in [("base", terms), ("GC3_adjusted", terms + ["GC3"])]:
        dd = d.assign(GC3=sc.GC3.values); g = genewise(Z, gk, pb.design(dd, tt, cats).values); g = g[~g.gene.isin(GC3)].copy()
        g["fdr"] = multipletests(g.p, method="fdr_bh")[1]; g["dataset"] = ds; g["cell_type"] = ct; g["model"] = lab; out[lab] = g
        u = g[g.gene.isin(MEAS)]; down = (u.fdr < 0.05) & (u.beta < 0); up = (u.fdr < 0.05) & (u.beta > 0)
        def orr(sig, S_):
            a = int((sig & u.gene.isin(S_)).sum()); b_ = int((sig & ~u.gene.isin(S_)).sum()); c_ = int((~sig & u.gene.isin(S_)).sum()); d_ = int((~sig & ~u.gene.isin(S_)).sum())
            o, p = stats.fisher_exact([[a, b_], [c_, d_]]) if (a + b_) > 0 else (np.nan, np.nan); return a, o, p
        a1, o1, p1 = orr(down, IND); a2, o2, p2 = orr(up, REP)
        t100d = u.nsmallest(100, "t"); t100u = u.nlargest(100, "t")
        AUD.append(dict(dataset=ds, cell_type=ct, model=lab, n_genes=len(g), n_genes_in_GCref=len(u), n_fdr05=int((g.fdr < 0.05).sum()), n_down=int(down.sum()), n_up=int(up.sum()),
                        down_GCinduced=a1, OR_down_GCinduced=o1, p_down_GCinduced=p1, up_GCrepressed=a2, OR_up_GCrepressed=o2, p_up_GCrepressed=p2,
                        top100down_frac_GCinduced=t100d.gene.isin(IND).mean(), top100up_frac_GCrepressed=t100u.gene.isin(REP).mean(),
                        bg_frac_GCinduced=u.gene.isin(IND).mean(), bg_frac_GCrepressed=u.gene.isin(REP).mean()))
    s0 = set(out["base"].gene[out["base"].fdr < 0.05]); s1 = set(out["GC3_adjusted"].gene[out["GC3_adjusted"].fdr < 0.05])
    AUD[-1]["jaccard_vs_base"] = len(s0 & s1) / max(1, len(s0 | s1)); AUD[-1]["n_lost"] = len(s0 - s1); AUD[-1]["n_gained"] = len(s1 - s0)
    GENEW.append(pd.concat(out.values()))
    return out
def block(ds, ct, d, sc, terms, cats, ref, do_auc=True):
    d = d.copy()
    for k in sc: d[k] = sc[k].values
    res = {}
    for adj, extra in [("base", []), ("+CLOCK", ["CLOCK"]), ("+CAMP", ["CAMP"]), ("+CLOCK+CAMP", ["CLOCK", "CAMP"])]:
        r = fit(ds, "C/K", ct, "GC3", adj, d, d.GC3.values, terms + extra, cats); res[adj] = -r["beta"]
    for m in ["CLOCK", "CAMP", "GCR"]: fit(ds, "C1/descriptive", ct, m, "base", d, d[m].values, terms, cats)
    if do_auc:
        for m in ["GC3", "CAMP", "CLOCK"]:
            xc, xk = d[m].values[ref], d[m].values[~ref]; u = stats.mannwhitneyu(xc, xk, alternative="two-sided")
            ROWS.append(dict(dataset=ds, test="K2_AUC", cell_type=ct, module=m, adjustment="none", n=len(d), n_case=int((~ref).sum()), beta=u.statistic / (len(xc) * len(xk)), ci_lo=np.nan, ci_hi=np.nan, p=u.pvalue, se_kind="AUC=P(ctrl>case)"))
        rr = {m: stats.spearmanr(d.GC3[ref], d[m][ref])[0] for m in ["CLOCK", "CAMP"]}
        ROWS.append(dict(dataset=ds, test="C3_rho_in_ref", cell_type=ct, module="GC3_vs_CLOCK", adjustment="none", n=int(ref.sum()), n_case=0, beta=rr["CLOCK"], ci_lo=np.nan, ci_hi=np.nan, p=np.nan, se_kind="spearman"))
        ROWS.append(dict(dataset=ds, test="C3_rho_in_ref", cell_type=ct, module="GC3_vs_CAMP", adjustment="none", n=int(ref.sum()), n_case=0, beta=rr["CAMP"], ci_lo=np.nan, ci_hi=np.nan, p=np.nan, se_kind="spearman"))
    return res
RET = {}

# ===================== [REVISION POST HOC] wrappers =====================
_audit0 = audit; _block0 = block
RNG = np.random.default_rng(20261003)
def matched_null(gkeep, top, nrep=2000):
    mean = pd.Series(LAST["mean"]).reindex(gkeep).dropna(); mean = mean[mean.index.isin(MEAS)]
    bins = pd.qcut(mean.rank(method="first"), 20, labels=False); tb = bins.reindex(top).dropna().astype(int).value_counts()
    pools = {b: mean.index[bins.values == b].values for b in tb.index}; obs = np.mean([g in IND for g in top])
    sims = np.array([np.mean([g in IND for b, k in tb.items() for g in RNG.choice(pools[b], k, replace=False)]) for _ in range(nrep)])
    return obs, sims.mean(), np.quantile(sims, 0.975), (1 + (sims >= obs).sum()) / (nrep + 1)
def topstats(ds, ct, lab, d, Z, gk, terms, cats):
    g = genewise(Z, gk, pb.design(d, terms, cats).values); g = g[~g.gene.isin(GC3)].copy(); g["fdr"] = multipletests(g.p, method="fdr_bh")[1]
    u = g[g.gene.isin(MEAS)]; top = list(u.nsmallest(100, "t").gene); down = (u.fdr < 0.05) & (u.beta < 0)
    a = int((down & u.gene.isin(IND)).sum()); b_ = int((down & ~u.gene.isin(IND)).sum()); c_ = int((~down & u.gene.isin(IND)).sum()); d_ = int((~down & ~u.gene.isin(IND)).sum())
    o = stats.fisher_exact([[a, b_], [c_, d_]])[0] if a + b_ else np.nan
    obs, nm, nhi, pe = matched_null(u.gene.values, top)
    NULL.append(dict(dataset=ds, cell_type=ct, model=lab, n_fdr05=int((g.fdr < 0.05).sum()), OR_down_GCinduced=o, top100_frac_GCinduced=obs, matched_null_mean=nm, matched_null_q975=nhi, p_empirical=pe))
def audit(ds, ct, d, Z, gk, terms, cats, sc):
    out = _audit0(ds, ct, d, Z, gk, terms, cats, sc)
    dd = d.assign(GC3=sc.GC3.values, MAC=sc.MAC.values)
    topstats(ds, ct, "base", dd, Z, gk, terms, cats)
    if ct == "bulk":
        topstats(ds, ct, "+MAC", dd, Z, gk, terms + ["MAC"], cats); topstats(ds, ct, "+MAC+GC3", dd, Z, gk, terms + ["MAC", "GC3"], cats)
    return out
def rfit(ds, test, ct, module, adj, d, y, terms, cats):
    D = pb.design(d, terms, cats); b, se, p, kind = pb.ols_hc3(D.values, y, 1, return_kind=True); t = stats.t.ppf(.975, len(y) - D.shape[1])
    REV.append(dict(dataset=ds, test=test, cell_type=ct, module=module, adjustment=adj, n=len(y), beta=b, ci_lo=b - t * se, ci_hi=b + t * se, p=p))
def block(ds, ct, d, sc, terms, cats, ref, do_auc=True):
    out = _block0(ds, ct, d, sc, terms, cats, ref, do_auc)
    dd = d.assign(**{k: sc[k].values for k in ["GC3", "GCR"]})
    if "age10" in dd.columns:
        r = dd[ref]; rfit(ds, "R3_GC3_age_in_reference", ct, "GC3", "sex+study", r, r.GC3.values, ["age10", "female"], ["study"] if r.study.nunique() > 1 else [])
        for stu, x in dd.groupby("study"):
            if x.case.sum() >= 3 and (x.case == 0).sum() >= 3:
                for m in ["GC3", "GCR"]: rfit(ds, "R5_per_study", ct + "|" + stu, m, "age+sex", x, -x[m].values, ["case", "age10", "female"], [])
    return out
# ===================== GSE290367 L1 (DCM vs ND, de-duplicated) =====================
PB = f"{ROOT}/02_data/processed/GSE290367_pseudobulk"; genes = pd.read_csv(f"{PB}/genes.csv").gene.astype(str).values
drop = set(pd.read_csv(f"{INP}/gse290367_duplicate_individuals.csv").record_drop)
def gse_load(level):
    m = pd.read_csv(f"{PB}/pseudobulk_{level}_meta.csv"); m["age"] = pd.to_numeric(m.age, errors="coerce"); m["row"] = np.arange(len(m))
    return m[(m.region == "LV") & (m.age_status == "postnatal") & (m.age >= 20) & ~m.donor_id.isin(drop)], np.load(f"{PB}/pseudobulk_{level}_counts.npy", mmap_mode="r")
def gse_sub(base, X, ct, dis, thr):
    s = base[(base.cell_type == ct) & base.disease.isin(["ND", dis]) & (base.n_cells >= thr[0]) & (base.total_umi >= thr[1])].groupby("donor_id").first()
    cnt = s.groupby("study").disease.agg(lambda x: ((x == dis).sum(), (x == "ND").sum())); s = s[s.study.isin([k for k, (a, b) in cnt.items() if a >= 1 and b >= 3])]
    d = s[["age", "sex", "study", "disease"]].copy(); d["case"] = (d.disease == dis).astype(float); d["female"] = (d.sex == "female").astype(float); d["age10"] = d.age / 10
    M = np.stack([np.asarray(X[r]).astype(float) for r in s.row.values]); lc, gk = pb.lognorm(M, genes, "CPM"); ref = (d.disease == "ND").values
    sc, Z, gk2 = scores(lc, gk, d.study.values, ref); return d, sc, Z, gk2, ref
b1, X1 = gse_load("L1")
for ct in ["Cardiomyocyte", "Fibroblast", "Endothelial", "Pericyte", "Myeloid"]:
    thr = (20, 5000) if ct == "Myeloid" else (50, 10000)
    d, sc, Z, gk, ref = gse_sub(b1, X1, ct, "DCM", thr); cats = ["study"] if d.study.nunique() > 1 else []
    RET[("GSE290367_DCM", ct)] = block("GSE290367_DCM", ct, d, sc, ["case", "age10", "female"], cats, ref)
    o = audit("GSE290367_DCM", ct, d, Z, gk, ["case", "age10", "female"], cats, sc)
    if ct in ["Cardiomyocyte", "Fibroblast"]:
        for g_, dirn in POS.items():
            for lab in ["base", "GC3_adjusted"]:
                x = o[lab][o[lab].gene == g_]
                if len(x): ROWS.append(dict(dataset="GSE290367_DCM", test="A3_descriptive", cell_type=ct, module=g_, adjustment=lab, n=len(d), n_case=int(d.case.sum()), beta=float(x.beta.iloc[0]), ci_lo=np.nan, ci_hi=np.nan, p=float(x.fdr.iloc[0]), se_kind=f"fdr;expected_sign={dirn}"))
b2, X2 = gse_load("L2")
d, sc, Z, gk, ref = gse_sub(b2, X2, "Myeloid:Resident", "DCM", (20, 5000)); d2 = d.copy()
for k in sc: d2[k] = sc[k].values
for adj, extra in [("base", []), ("+CLOCK", ["CLOCK"]), ("+CAMP", ["CAMP"]), ("+CLOCK+CAMP", ["CLOCK", "CAMP"])]:
    fit("GSE290367_Resident_DCM", "C/K_G4", "Myeloid:Resident", "G4", adj, d2, d2.G4.values, ["case", "age10", "female"] + extra, ["study"] if d2.study.nunique() > 1 else [])
# ===================== SCP3342 (HFpEF vs control) =====================
S_ = f"{ROOT}/02_data/processed/SCP3342_pseudobulk"; sg = pd.read_csv(f"{S_}/genes.csv").gene.astype(str).values; sm = pd.read_csv(f"{S_}/pseudobulk_meta.csv", dtype={"donor_id": str})
C_ = np.load(f"{S_}/pseudobulk_raw_count_cellbender.npy", mmap_mode="r")
SMAP = {"Cardiomyocyte": ["Cardiomyocyte"], "Fibroblast": ["Fibroblast"], "Endothelial": ["Endothelial1", "Endothelial2"], "Pericyte": ["Pericyte"],
        "Myeloid": ["Macrophage", "Proliferating_macrophage"], "VascStromal": ["Endothelial1", "Endothelial2", "Fibroblast", "Pericyte", "VSMC"]}
for ct, labs in SMAP.items():
    ss = sm[sm.cell_type.isin(labs)].copy(); ss["row"] = ss.index
    g = ss.groupby("donor_id").agg(n=("n_cells", "sum"), umi=("total_umi_raw_count_cellbender", "sum"), disease=("disease", "first"), sex=("sex", "first"), pool=("pool", "first"))
    thr = (20, 5000) if ct == "Myeloid" else (50, 10000); g = g[(g.n >= thr[0]) & (g.umi >= thr[1])]
    M = np.stack([np.asarray(C_[ss[ss.donor_id == dn].row.values]).astype(float).sum(0) for dn in g.index])
    g["case"] = (g.disease == "hfpef").astype(float); g["female"] = (g.sex == "female").astype(float); ref = (g.case == 0).values
    lc, gk = pb.lognorm(M, sg, "CPM"); sc, Z, gk2 = scores(lc, gk, np.zeros(len(g)), ref)
    if ct != "VascStromal": RET[("SCP3342_HFpEF", ct)] = block("SCP3342_HFpEF", ct, g, sc, ["case", "female"], ["pool"], ref)
    audit("SCP3342_HFpEF", ct, g, Z, gk2, ["case", "female"], ["pool"], sc)
# ===================== MAGNet (overlap-excluded DCM vs NF) =====================
gtx = pd.read_csv(f"{ROOT}/02_data/raw/GTEx_v11/rerelease_2026-05/gene_reads_adult_gtex_v11_heart_left_ventricle.gct.gz", sep="\t", skiprows=2, usecols=["Name", "Description"])
E2S = pd.Series(gtx.Description.values, index=gtx.Name.str.split(".").str[0].values); E2S = E2S[~E2S.index.duplicated()]
t = tarfile.open(f"{ROOT}/02_data/raw/GSE141910_MAGNet/GSE141910_RAW.tar"); cols = []
for mb in t.getmembers():
    s = pd.read_csv(io.BytesIO(gzip.decompress(t.extractfile(mb).read())), index_col=0); s.columns = [s.columns[0]]; cols.append(s)
X = pd.concat(cols, axis=1); ph = pd.read_csv(f"{INP}/MAGNet_phenoData.csv"); fl = pd.read_csv(f"{INP}/magnet_overlap_flags.csv")
ph = ph.merge(fl[["id", "flag_possible_overlap"]], left_on="sample_name", right_on="id", how="left").set_index("sample_name").loc[X.columns]; ph["male"] = (ph.gender == "Male").astype(float)
lin = 2 ** X; sym = lin.index.str.split(".").str[0].map(E2S); lin = lin[sym.notna()].groupby(sym[sym.notna()].values).sum()
keep = (~ph.flag_possible_overlap.fillna(False).astype(bool)) & ph.etiology.isin(["NF", "DCM"]); ph = ph[keep.values]; lin = lin.loc[:, keep.values]
Mm = lin.T.values; cpm = Mm / Mm.sum(1, keepdims=True) * 1e6; kg = (cpm > 1).mean(0) >= 0.2; lcM = np.log2(cpm[:, kg] + 1); gM = lin.index.values[kg]
ph["case"] = (ph.etiology == "DCM").astype(float); refM = (ph.etiology == "NF").values
scM, ZM, gM2 = scores(lcM, gM, np.zeros(len(ph)), refM); dM = ph.copy()
for k in scM: dM[k] = scM[k].values
for adj, extra in [("base", []), ("+CLOCK", ["CLOCK"]), ("+CAMP", ["CAMP"]), ("+CLOCK+CAMP", ["CLOCK", "CAMP"])]:
    fit("MAGNet_exclOverlap_DCM", "C/K_G4", "bulk", "G4", adj, dM, dM.G4.values, ["case", "age", "male"] + extra, [])
    fit("MAGNet_exclOverlap_DCM", "C/K_bulk", "bulk", "GC3", adj, dM, dM.GC3.values, ["case", "age", "male"] + extra, [])
for m in ["CLOCK", "CAMP", "GCR"]: fit("MAGNet_exclOverlap_DCM", "C1/descriptive", "bulk", m, "base", dM, dM[m].values, ["case", "age", "male"], [])
oM = audit("MAGNet_exclOverlap_DCM", "bulk", ph, ZM, gM2, ["case", "age", "male"], [], scM)
A3 = []
for g_, dirn in POS.items():
    for lab in ["base", "GC3_adjusted"]:
        x = oM[lab][oM[lab].gene == g_]
        A3.append(dict(gene=g_, expected=dirn, model=lab, beta=float(x.beta.iloc[0]) if len(x) else np.nan, fdr=float(x.fdr.iloc[0]) if len(x) else np.nan))
A3 = pd.DataFrame(A3); A3["ok"] = (np.sign(A3.beta) == A3.expected) & (A3.fdr < 0.05)
DEC["A3"] = dict(n_ok_adjusted=int(A3[A3.model == "GC3_adjusted"].ok.sum()), n_ok_base=int(A3[A3.model == "base"].ok.sum()), PASS=bool(A3[A3.model == "GC3_adjusted"].ok.sum() >= 6), table=A3.to_dict("records"))
# ===================== Hahn 2021 RVS (HFpEF vs control) =====================
HD = f"{ROOT}/02_data/raw/Hahn2021_HFpEF_EMB_bulk"
cl = pd.read_csv(f"{HD}/ClinicalCharacteristicsRNAseqSamples.csv", encoding="latin1"); cl.columns = [c.replace("\ufeff", "").replace("ï»¿", "").strip() for c in cl.columns]
st = pd.read_csv(f"{HD}/SamplesTable_ALLpatients_sex_race_age.csv"); st["sample"] = st.samplenames.astype(str)
mt_ = st[st.tissue == "RVS"].copy(); mt_["group"] = mt_.disease.map({"Normal": "Control", "HFpEF": "HFpEF", "HFrEF": "HFrEF"}); mt_ = mt_[mt_.group.isin(["Control", "HFpEF"])].set_index("sample")
raw = pd.read_excel(f"{HD}/RawReads_RVS_Ensembl_geneName.xlsx", index_col=0); raw.columns = raw.columns.astype(str); raw = raw[[c for c in mt_.index if c in raw.columns]]; mt_ = mt_.loc[raw.columns]
sym = raw.index.str.split(".").str[0].map(E2S); cnt = raw[sym.notna()].groupby(sym[sym.notna()].values).sum()
lc, gk = pb.lognorm(cnt.T.values.astype(float), cnt.index.values, "TMM"); mt_["case"] = (mt_.group == "HFpEF").astype(float); mt_["male"] = (mt_.sex == "Male").astype(float); refH = (mt_.group == "Control").values
scH, ZH, gH = scores(lc, gk, np.zeros(len(mt_)), refH); audit("Hahn2021_RVS_HFpEF", "bulk", mt_, ZH, gH, ["case", "age", "male"], [], scH)

# ===================== [REVISION POST HOC] macrophage modules with GC3 adjustment; age in reference; Hahn medications =====================
for ds, dd_, terms, cats in [("GSE290367_Resident_DCM", d2, ["case", "age10", "female"], ["study"] if d2.study.nunique() > 1 else []), ("MAGNet_exclOverlap_DCM", dM, ["case", "age", "male"], [])]:
    for m in ["G3", "G4", "RES", "MHC"]:
        for adj, extra in [("base", []), ("+GC3", ["GC3"]), ("+GCR", ["GCR"]), ("+GC3+MAC", ["GC3", "MAC"])]:
            rfit(ds, "R2_macrophage_GC3adj", "bulk" if ds.startswith("MAG") else "Myeloid:Resident", m, adj, dd_, dd_[m].values, terms + extra, cats)
rM = dM[refM]; rfit("MAGNet_exclOverlap_DCM", "R3_GC3_age_in_reference", "bulk", "GC3", "sex", rM, rM.GC3.values, ["age", "male"], [])
dH = mt_.copy()
for k in scH: dH[k] = scH[k].values
rH = dH[refH]; rfit("Hahn2021_RVS_HFpEF", "R3_GC3_age_in_reference", "bulk", "GC3", "sex", rH, rH.GC3.values, ["age", "male"], [])
rfit("Hahn2021_RVS_HFpEF", "R2_GC3_case_adj_age", "bulk", "GC3", "age+sex", dH, -dH.GC3.values, ["case", "age", "male"], [])
try:
    clr = cl[cl["Tissue Abbr"].astype(str).str.upper() == "RVS"].copy(); clr["sample"] = clr.SequenceRunId.astype(str)
    hx = dH[dH.case == 1].join(clr.set_index("sample")[["ACEi/ARB", "Beta Blocker", "Loop Diuretic"]], how="inner")
    REV.append(dict(dataset="Hahn2021_RVS_HFpEF", test="R4_med_match", cell_type="bulk", module="n_matched_cases", adjustment="", n=len(hx), beta=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan))
    for med in ["ACEi/ARB", "Beta Blocker", "Loop Diuretic"]:
        x = hx.dropna(subset=[med]).copy(); x["med"] = (x[med].astype(str) == "Yes").astype(float)
        if 3 <= x.med.sum() <= len(x) - 3: rfit("Hahn2021_RVS_HFpEF", "R4_med_within_HFpEF", "bulk", "GC3", med, x, x.GC3.values, ["med", "age", "male"], [])
except Exception as e:
    REV.append(dict(dataset="Hahn2021_RVS_HFpEF", test="R4_error", cell_type=str(e)[:80], module="", adjustment="", n=0, beta=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan))
pd.DataFrame(REV).to_csv(f"{OUT}/rev_posthoc_models.csv", index=False); pd.DataFrame(NULL).to_csv(f"{OUT}/rev_posthoc_audit_null.csv", index=False)
print(pd.DataFrame(REV).round(4).to_string()); print(pd.DataFrame(NULL).round(4).to_string()); raise SystemExit(0)
# ===================== decisions C2 / K1 =====================
R = pd.DataFrame(ROWS)
for nm, adj in [("C2", "+CLOCK"), ("K1", "+CAMP"), ("CK_sensitivity", "+CLOCK+CAMP")]:
    rat = {}
    for ds in ["SCP3342_HFpEF", "GSE290367_DCM"]:
        rr = [RET[(ds, ct)][adj] / RET[(ds, ct)]["base"] for ct in ["Cardiomyocyte", "Fibroblast", "Endothelial", "Pericyte", "Myeloid"] if RET[(ds, ct)]["base"] > 0]
        rat[ds] = dict(retention=[round(x, 3) for x in rr], median=float(np.median(rr)) if rr else np.nan, n_celltypes=len(rr))
    g4 = R[(R.test == "C/K_G4") & (R.adjustment == adj)][["dataset", "beta", "ci_lo", "ci_hi", "p"]].to_dict("records")
    ok = all(v["median"] >= 0.5 for v in rat.values()) and all((x["beta"] < 0) and (x["p"] < 0.05) for x in g4) and len(g4) == 2
    DEC[nm] = dict(retention=rat, G4=g4, PASS=bool(ok) if nm != "CK_sensitivity" else None, rule="median retention >=0.5 in both datasets AND G4 beta<0,P<0.05 in GSE Resident DCM and MAGNet exclOverlap DCM")
DEC["coverage_last"] = COVG
R.to_csv(f"{OUT}/v13_score_tests.csv", index=False); pd.DataFrame(AUD).to_csv(f"{OUT}/v13_impact_audit.csv", index=False)
pd.concat(GENEW).to_csv(f"{OUT}/v13_genewise_before_after.csv.gz", index=False)
json.dump(DEC, open(f"{OUT}/v13_decisions.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print(json.dumps({k: v for k, v in DEC.items() if k != "coverage_last"}, indent=1, default=str)[:5000])
print(pd.DataFrame(AUD).round(3).to_string())
