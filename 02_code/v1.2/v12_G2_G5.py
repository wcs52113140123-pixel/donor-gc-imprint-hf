#!/usr/bin/env python3
"""v1.2 pre-registered tests G3/G4 (decision), MHC-II secondary, G5 (decision), G2 (supportive) — ANALYSIS_PLAN_v1.2.md
(sha256 26db5075acb7d1b859e078690b4a8cfd536bc09f95c99b0c975e069a719602b2). Gene classes from the in-vitro GC reference
(gc_gene_classification.csv; 4 QC-passing datasets). G3-primary (RES_EFF_core ∩ insensitive) has 1 gene (LRP1) -> not testable;
per plan G3-secondary (RES_EFF_lit ∩ insensitive, >=3 genes) becomes the decision test. Thresholds are NOT changed after results.
Run on server: cd <root>; source env_analysis.sh; INP=<staged dir> python v12_G2_G5.py <outdir>"""
import os, sys, io, gzip, tarfile, json, hashlib, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
ROOT = "/home/weichangsheng/projects/cardiac_ageing_target_clearance"; sys.path.insert(0, f"{ROOT}/03_code/04_models"); import pbtools_v11 as pb
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True); INP = os.environ.get("INP", ".")
assert hashlib.sha256(open(f"{ROOT}/01_config/ANALYSIS_PLAN_v1.2.md", "rb").read()).hexdigest() == "26db5075acb7d1b859e078690b4a8cfd536bc09f95c99b0c975e069a719602b2"
MODF = f"{ROOT}/01_config/modules/modules_v1.1.tsv"; assert hashlib.sha256(open(MODF, "rb").read()).hexdigest().startswith("3a2abdad")
mods = pd.read_csv(MODF, sep="\t"); MOD = lambda m: sorted(set(mods.gene[mods.module == m]))
cls = pd.read_csv(f"{INP}/gc_gene_classification.csv")
GC = ["FKBP5", "TSC22D3", "ZBTB16", "PER1", "KLF9", "DDIT4", "KLF15", "MT2A"]
lit = cls[cls.module == "RES_EFF_lit"]; core = cls[cls.module == "RES_EFF_core"]
G3P = sorted(core.gene[core["class"] == "GC_insensitive"]); G3S = sorted(lit.gene[lit["class"] == "GC_insensitive"]); G4 = sorted(lit.gene[lit["class"] == "GC_induced"])
MH = cls[cls.module == "MHCII_core"]; MH_INS = sorted(MH.gene[MH["class"] == "GC_insensitive"])
SETS = {"G3": G3S if len(G3P) < 3 else G3P, "G4": G4, "GC_response": GC, "RES_EFF_core": MOD("RES_EFF_core"), "MHCII_core": MOD("MHCII_core"), "MAC_PAN_ID": MOD("MAC_PAN_ID")}
META = dict(G3_primary_genes=G3P, G3_primary_testable=len(G3P) >= 3, G3_used=("secondary" if len(G3P) < 3 else "primary"), G3_genes=SETS["G3"], G4_genes=G4,
            MHCII_insensitive=MH_INS, MHCII_criterion3_applicable=len(MH_INS) >= 3)
assert len(SETS["G3"]) >= 3 and len(G4) >= 3
ROWS = []
def fit(ds, test, contrast, module, adj, d, y, terms, cats, extra=None, tost=False):
    D = pb.design(d, terms, cats); b, se, p, kind = pb.ols_hc3(D.values, y, 1, return_kind=True); df = len(y) - D.shape[1]; t = stats.t.ppf(.975, df)
    r = dict(dataset=ds, test=test, contrast=contrast, module=module, adjustment=adj, n=len(y), n_case=int(d[terms[0]].sum()) if set(d[terms[0]].unique()) <= {0, 1} else np.nan,
             beta=b, ci_lo=b - t * se, ci_hi=b + t * se, se=se, p=p, se_kind=kind)
    if tost: r["p_tost"] = max(stats.t.sf((b + 0.1) / se, df), stats.t.cdf((b - 0.1) / se, df))
    if extra: r.update(extra)
    ROWS.append(r); return r
def scores(lc, gk, groups, ref):
    Z = pb.ref_z(lc, groups, ref); ok = ~np.isnan(Z).any(0); Z, gk = Z[:, ok], gk[ok]
    return pd.DataFrame({k: pb.module_score(Z, gk, v, ref)[0] for k, v in SETS.items()}), Z, gk, lc[:, ok].mean(0), {k: pb.module_score(Z, gk, v, ref)[1] for k, v in SETS.items()}
COV, COMP = {}, []
# ===================== GSE290367 myeloid subsets (L2), de-duplicated =====================
PB = f"{ROOT}/02_data/processed/GSE290367_pseudobulk"; genes = pd.read_csv(f"{PB}/genes.csv").gene.astype(str).values
L2 = pd.read_csv(f"{PB}/pseudobulk_L2_meta.csv"); L2["age"] = pd.to_numeric(L2.age, errors="coerce"); L2["row"] = np.arange(len(L2))
X2 = np.load(f"{PB}/pseudobulk_L2_counts.npy", mmap_mode="r"); drop = set(pd.read_csv(f"{INP}/gse290367_duplicate_individuals.csv").record_drop)
base = L2[(L2.region == "LV") & (L2.age_status == "postnatal") & (L2.age >= 20) & ~L2.donor_id.isin(drop)]
DEC = {}
for sub in ["Myeloid:Resident", "Myeloid:Monocyte_derived"]:
    for dis in ["DCM", "HCM"]:
        s = base[(base.cell_type == sub) & base.disease.isin(["ND", dis]) & (base.n_cells >= 20) & (base.total_umi >= 5000)].groupby("donor_id").first()
        cnt = s.groupby("study").disease.agg(lambda x: ((x == dis).sum(), (x == "ND").sum())); keep = [k for k, (a, b) in cnt.items() if a >= 1 and b >= 3]; s = s[s.study.isin(keep)]
        d = s[["age", "sex", "study", "disease"]].copy(); d["case"] = (d.disease == dis).astype(float); d["female"] = (d.sex == "female").astype(float); d["age10"] = d.age / 10
        M = np.stack([np.asarray(X2[r]).astype(float) for r in s.row.values]); lc, gk = pb.lognorm(M, genes, "CPM"); ref = (d.disease == "ND").values
        sc, Z, gk2, expr, cov = scores(lc, gk, d.study.values, ref); COV[f"GSE290367_{sub}_{dis}"] = cov
        for k in sc: d[k] = sc[k].values
        cats = ["study"] if d.study.nunique() > 1 else []; ds = f"GSE290367_{sub.split(':')[1]}"; contrast = f"{dis}_vs_ND"
        decision = (sub == "Myeloid:Resident" and dis == "DCM")
        for m in ["G3", "G4", "RES_EFF_core", "MHCII_core", "GC_response"]:
            fit(ds, "G3/G4/MHCII", contrast, m, "age+sex+study", d, d[m].values, ["case", "age10", "female"], cats)
            if m != "GC_response":
                r = fit(ds, "G3/G4/MHCII", contrast, m, "age+sex+study+GC", d, d[m].values, ["case", "age10", "female", "GC_response"], cats,
                        extra=dict(decision_test=decision and m == "G3"))
                lowref = ref & (d.GC_response.values <= np.median(d.GC_response.values[ref])); kk = lowref | (d.case.values == 1)
                fit(ds, "G3/G4/MHCII", contrast + "_GClow_ref", m, "age+sex+study", d[kk], d[m].values[kk], ["case", "age10", "female"], cats, extra=dict(n_ref_GClow=int(lowref.sum())))
        if decision:
            Dz = pb.design(d, ["case", "age10", "female", "GC_response"], cats).values
            c = pb.competitive(Z, gk2, SETS["G3"], expr, Dz, j=1, n=1000, rng=np.random.default_rng(20261003), is_ref=ref); COMP.append(dict(dataset=ds, contrast=contrast, module="G3", adjustment="age+sex+study+GC", **c))
            DEC["GSE_rho"] = {m: stats.spearmanr(d[m], d.GC_response)[:2] for m in ["G3", "G4"]}
            d.to_csv(f"{OUT}/GSE290367_Resident_DCM_scores.csv")
# ===================== MAGNet (overlap-excluded for decision) =====================
gtx = pd.read_csv(f"{ROOT}/02_data/raw/GTEx_v11/rerelease_2026-05/gene_reads_adult_gtex_v11_heart_left_ventricle.gct.gz", sep="\t", skiprows=2)
E2S = pd.Series(gtx.Description.values, index=gtx.Name.str.split(".").str[0].values); E2S = E2S[~E2S.index.duplicated()]
t = tarfile.open(f"{ROOT}/02_data/raw/GSE141910_MAGNet/GSE141910_RAW.tar"); cols = []
for mb in t.getmembers():
    s = pd.read_csv(io.BytesIO(gzip.decompress(t.extractfile(mb).read())), index_col=0); s.columns = [s.columns[0]]; cols.append(s)
X = pd.concat(cols, axis=1)
ph = pd.read_csv(f"{INP}/MAGNet_phenoData.csv"); fl = pd.read_csv(f"{INP}/magnet_overlap_flags.csv")
ph = ph.merge(fl[["id", "flag_possible_overlap"]], left_on="sample_name", right_on="id", how="left").set_index("sample_name").loc[X.columns]
ph["male"] = (ph.gender == "Male").astype(float)
lin = 2 ** X; sym = lin.index.str.split(".").str[0].map(E2S); lin = lin[sym.notna()].groupby(sym[sym.notna()].values).sum()
Mm = lin.T.values; cpm = Mm / Mm.sum(1, keepdims=True) * 1e6; keepg = (cpm > 1).mean(0) >= 0.2; lcM = np.log2(cpm[:, keepg] + 1); gM = lin.index.values[keepg]
isnf = (ph.etiology == "NF").values; scM, ZM, gM2, exprM, covM = scores(lcM, gM, np.zeros(len(ph)), isnf); COV["MAGNet"] = covM
for k in scM: ph[k] = scM[k].values
excl = ~ph.flag_possible_overlap.fillna(False).astype(bool).values
for aset, keepset in [("excl_overlap", excl), ("all", np.ones(len(ph), bool))]:
    for dis in ["DCM", "HCM"]:
        sub = keepset & ph.etiology.isin(["NF", dis]).values; d = ph[sub].copy(); d["case"] = (d.etiology == dis).astype(float)
        decision = (aset == "excl_overlap" and dis == "DCM"); ds = f"MAGNet_{aset}"; contrast = f"{dis}_vs_NF"
        for m in ["G3", "G4", "RES_EFF_core", "MHCII_core", "GC_response"]:
            fit(ds, "G3/G4/MHCII", contrast, m, "age+sex", d, d[m].values, ["case", "age", "male"], [])
            if m != "GC_response":
                fit(ds, "G3/G4/MHCII", contrast, m, "age+sex+GC+MAC", d, d[m].values, ["case", "age", "male", "GC_response", "MAC_PAN_ID"], [], extra=dict(decision_test=decision and m == "G3"))
                fit(ds, "G3/G4/MHCII", contrast, m, "age+sex+GC", d, d[m].values, ["case", "age", "male", "GC_response"], [])
        if decision:
            Dz = pb.design(d, ["case", "age", "male", "GC_response", "MAC_PAN_ID"], []).values
            c = pb.competitive(ZM[sub], gM2, SETS["G3"], exprM, Dz, j=1, n=1000, rng=np.random.default_rng(20261003), is_ref=isnf[sub]); COMP.append(dict(dataset=ds, contrast=contrast, module="G3", adjustment="age+sex+GC+MAC", **c))
            DEC["MAGNet_rho"] = {m: stats.spearmanr(d[m], d.GC_response)[:2] for m in ["G3", "G4"]}
            d.to_csv(f"{OUT}/MAGNet_exclOverlap_DCM_scores.csv")
# ===================== GTEx: G2 (supportive) and G5 (decision) =====================
subj = pd.read_csv(f"{ROOT}/02_data/raw/GTEx_v11/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt", sep="\t")
samp = pd.read_csv(f"{ROOT}/02_data/raw/GTEx_v11/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt", sep="\t", low_memory=False, usecols=["SAMPID", "SMCENTER", "SMRIN", "SMTSISCH"])
for tissue, f in [("LV", "left_ventricle"), ("AA", "atrial_appendage")]:
    g = gtx if tissue == "LV" else pd.read_csv(f"{ROOT}/02_data/raw/GTEx_v11/rerelease_2026-05/gene_reads_adult_gtex_v11_heart_{f}.gct.gz", sep="\t", skiprows=2)
    cnt = g.drop(columns=["Name"]).groupby("Description").sum().astype(float); lc, gk = pb.lognorm(cnt.T.values, cnt.index.values, "TMM")
    d = pd.DataFrame({"SAMPID": cnt.columns}); d["SUBJID"] = d.SAMPID.str.split("-").str[:2].str.join("-"); d = d.merge(samp, on="SAMPID", how="left").merge(subj, on="SUBJID", how="left")
    d["age_dec"] = (d.AGE.str[:2].astype(float) + 4.5) / 10; d["age"] = d.age_dec * 10; d["male"] = (d.SEX == 1).astype(float); d["vent"] = (d.DTHHRDY == 0).astype(float)
    d["hardy"] = d.DTHHRDY.fillna(-1).astype(int).astype(str)
    sc, Z, gk2, expr, cov = scores(lc, gk, np.zeros(len(d)), np.ones(len(d), bool)); COV[f"GTEx_{tissue}"] = cov
    for k in sc: d[k] = sc[k].values
    ok = d[["SMRIN", "SMTSISCH", "SMCENTER", "DTHHRDY"]].notna().all(1).values; dd = d[ok].copy()
    for m in ["GC_response", "RES_EFF_core", "MHCII_core"]:
        fit(f"GTEx_{tissue}", "G2", "ventilator_vs_other", m, "age+sex+RIN+isch+center+MAC", dd, dd[m].values, ["vent", "age", "male", "SMRIN", "SMTSISCH", "MAC_PAN_ID"], ["SMCENTER"])
    for m in ["RES_EFF_core", "MHCII_core"]:
        fit(f"GTEx_{tissue}", "G2", "per_SD_GC_response", m, "age+sex+RIN+isch+center+MAC", dd.assign(gcz=(dd.GC_response - dd.GC_response.mean()) / dd.GC_response.std()), dd[m].values,
            ["gcz", "age", "male", "SMRIN", "SMTSISCH", "MAC_PAN_ID"], ["SMCENTER"])
    for aset, kk in [("non_ventilator_DTHHRDY1-4", dd.DTHHRDY.isin([1, 2, 3, 4]).values), ("ventilator_only", (dd.DTHHRDY == 0).values)]:
        s2 = dd[kk]; cats = ["SMCENTER", "hardy"] if aset.startswith("non") else ["SMCENTER"]
        for m in ["RES_EFF_core", "MHCII_core"]:
            fit(f"GTEx_{tissue}", "G5", f"age_per_decade_{aset}", m, "sex+RIN+isch+center+hardy+MAC", s2, s2[m].values, ["age_dec", "male", "SMRIN", "SMTSISCH", "MAC_PAN_ID"], cats,
                tost=True, extra=dict(decision_test=(m == "RES_EFF_core" and aset.startswith("non"))))
R = pd.DataFrame(ROWS); C = pd.DataFrame(COMP)
# ---------------- decisions ----------------
dt = R[R.get("decision_test", False) == True]
g3 = dt[(dt.test == "G3/G4/MHCII") & (dt.module == "G3")].copy(); g3["p_holm"] = multipletests(g3.p, method="holm")[1]
G3_PASS = bool((g3.beta < 0).all() and (g3.p_holm < 0.05).all() and len(g3) == 2 and (C.p_comp < 0.05).any())
def att(ds, contrast, base_adj, gc_adj):
    a = R[(R.dataset == ds) & (R.contrast == contrast) & (R.module == "G4")].set_index("adjustment").beta; return float(1 - a[gc_adj] / a[base_adj])
G4 = dict(GSE_attenuation=att("GSE290367_Resident", "DCM_vs_ND", "age+sex+study", "age+sex+study+GC"), MAGNet_attenuation=att("MAGNet_excl_overlap", "DCM_vs_NF", "age+sex", "age+sex+GC"),
          rho=dict(GSE=DEC["GSE_rho"], MAGNet=DEC["MAGNet_rho"]))
G4["prediction_a_rho_G4_gt_G3_both"] = bool(abs(DEC["GSE_rho"]["G4"][0]) > abs(DEC["GSE_rho"]["G3"][0]) and abs(DEC["MAGNet_rho"]["G4"][0]) > abs(DEC["MAGNet_rho"]["G3"][0]))
G4["prediction_b_attenuation_ge50_both"] = bool(G4["GSE_attenuation"] >= 0.5 and G4["MAGNet_attenuation"] >= 0.5)
mh = R[(R.module == "MHCII_core") & (((R.dataset == "GSE290367_Resident") & (R.contrast == "DCM_vs_ND") & (R.adjustment == "age+sex+study+GC")) |
                                      ((R.dataset == "MAGNet_excl_overlap") & (R.contrast == "DCM_vs_NF") & (R.adjustment == "age+sex+GC")))]
MHC = dict(rows=mh[["dataset", "beta", "p"]].to_dict("records"), criterion3_applicable=META["MHCII_criterion3_applicable"], retained=bool((mh.p < 0.05).all() and (mh.beta > 0).all() and len(mh) == 2))
g5 = R[(R.test == "G5") & (R.get("decision_test", False) == True)].copy(); g5["p_tost_holm"] = multipletests(g5.p_tost, method="holm")[1]
G5_PASS = bool(g5.loc[g5.dataset == "GTEx_LV", "p_tost"].iloc[0] < 0.05)
DECISION = dict(meta=META, G3=dict(rows=g3[["dataset", "beta", "ci_lo", "ci_hi", "p", "p_holm", "n", "n_case"]].to_dict("records"), competitive=C.to_dict("records"), PASS=G3_PASS),
                G4=G4, MHCII=MHC, G5=dict(rows=g5[["dataset", "beta", "ci_lo", "ci_hi", "p_tost", "p_tost_holm", "n"]].to_dict("records"), PASS_LV_rawTOST=G5_PASS), coverage=COV)
R.to_csv(f"{OUT}/G2_G5_results.csv", index=False); C.to_csv(f"{OUT}/G3_competitive.csv", index=False)
json.dump(DECISION, open(f"{OUT}/G2_G5_decisions.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
pd.set_option("display.width", 250); print(R.drop(columns=[c for c in ["se"] if c in R]).round(4).to_string()); print(json.dumps(DECISION, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))[:6000])
