#!/usr/bin/env python3
# [v12-wrap] Post-processing step originally executed in the analysis sandbox (not on the server); code recovered from
# artifact lineage of GSE290367_decomposition_shares.csv (version e2c094dc-2676-4cb8-87a8-7f3b48f09dce). Only input/output paths were parameterised:
#   inputs  -> files under $IN_DIR (rerun copies of the upstream server outputs), outputs -> cwd ($OUT_DIR).
import os; IN = os.environ["IN_DIR"].rstrip("/") + "/"; os.chdir(os.environ["OUT_DIR"])
import numpy as np
import pandas as pd

def bh(p):
    p=np.asarray(p,float); ok=~np.isnan(p); out=np.full_like(p,np.nan); pp=p[ok]
    if len(pp)==0: return out
    o=np.argsort(pp); q=pp[o]*len(pp)/np.arange(1,len(pp)+1); q=np.minimum.accumulate(q[::-1])[::-1]; t=np.empty_like(pp); t[o]=np.minimum(q,1); out[ok]=t; return out

RULE='n_disease>=5 & n_ref>=8 & df_resid>=5 & finite CI (disease/composition); ageing: n>=15'

def estim(r):
    if str(r.get('contrast',''))=='per_decade_ND' or r.analysis.startswith('ND_ageing'):
        return bool(r.n_ref>=15 and np.isfinite(r.ci_lo))
    return bool(r.n_disease>=5 and r.n_ref>=8 and r.df_resid>=5 and np.isfinite(r.ci_lo) and np.isfinite(r.ci_hi) and r.se>0)

G=pd.read_csv('' + IN + 'substates_GSE290367_effects.csv')
G['estimable']=G.apply(estim,axis=1); G['estimability_rule']=RULE
G['q_BH_family']=np.nan; G['q_TOST_BH_family']=np.nan
for (f,m),idx in G[G.estimable].groupby(['family','model']).groups.items():
    G.loc[idx,'q_BH_family']=bh(G.loc[idx,'p'])
    if G.loc[idx,'p_tost'].notna().all(): G.loc[idx,'q_TOST_BH_family']=bh(G.loc[idx,'p_tost'])

d=G[G.analysis.str.startswith('decomposition')&G.contrast.ne('per_decade_ND')].pivot_table(index=['module','contrast'],columns='analysis',values='beta')
d.columns=[c.replace('decomposition_','') for c in d.columns]; d['interaction_residual']=d.observed-d.composition_only-d.state_only
d['composition_share']=d.composition_only/d.observed; d['state_share']=d.state_only/d.observed
d.round(4).to_csv('GSE290367_decomposition_shares.csv')