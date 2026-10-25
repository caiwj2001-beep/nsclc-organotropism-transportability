# -*- coding: utf-8 -*-
# Path configuration for the public repository.
#   CHORD_DIR : directory holding the cBioPortal msk_chord_2024 export (comma-delimited CSVs)
#   GENE_MAP  : entrezGeneId -> hugoGeneSymbol mapping (2 columns, exported from _gene_map.csv)
#   OUT_DIR   : working/output directory for intermediate and result tables
import os as _os
_CHORD = _os.environ.get("CHORD_DIR", "data/msk_chord_2024")
_GENEMAP = _os.environ.get("GENE_MAP", "data/gene_map.csv")
_OUT = _os.environ.get("OUT_DIR", "results")
_os.makedirs(_OUT, exist_ok=True)
"""STAGE 2C: burden-adjusted tropism — organ distribution among SINGLE-organ patients,
and number-of-organs distribution by histology."""
import os, numpy as np, pandas as pd
from scipy.stats import chi2_contingency
OUT = _os.path.abspath(_OUT)
pat=pd.read_csv(os.path.join(OUT,"stage2_patient_table.csv"),index_col=0,low_memory=False)
DISTANT=["BONE","LIVER","CNS_BRAIN","PLEURA","ADRENAL_GLANDS","INTRA_ABDOMINAL","REPRODUCTIVE_ORGANS","OTHER"]
print("n_org distribution by histology")
for g in ["LUAD","LUSC","NSCLC-NOS"]:
    d=pat[(pat.grp==g)&pat.any_distant]
    print(f"  {g}: n={len(d)} | n_org mean={d.n_distant.mean():.2f} | dist={d.n_distant.value_counts().sort_index().to_dict()}")
print("\nOrgan distribution among SINGLE-organ patients (the burden-free tropism comparison)")
rows={}
for g in ["LUAD","LUSC","NSCLC-NOS"]:
    d=pat[(pat.grp==g)&(pat.n_distant==1)]
    rows[g]={c:int(d[c+"_b"].sum()) for c in DISTANT}|{"n":len(d)}
t=pd.DataFrame(rows).T
pct=(t.drop(columns="n").div(t.n,axis=0)*100).round(1)
print(t.to_string()); print("\n%:"); print(pct.to_string())
# 7-organ chi-square (drop OTHER to avoid heterogeneity)
cols=["BONE","LIVER","CNS_BRAIN","PLEURA","ADRENAL_GLANDS","INTRA_ABDOMINAL","REPRODUCTIVE_ORGANS"]
ct=t.loc[["LUAD","LUSC"],cols].values
chi2,p,dof,_=chi2_contingency(ct)
print(f"\nLUAD vs LUSC single-organ 7-organ chi2={chi2:.1f} dof={dof} p={p:.3g}")
# pairwise per organ (BH)
from statsmodels.stats.multitest import multipletests
ps=[];names=[]
for c in cols:
    a=t.loc["LUAD",c]; n1=t.loc["LUAD","n"]; b=t.loc["LUSC",c]; n0=t.loc["LUSC","n"]
    from scipy.stats import fisher_exact
    orr,pv=fisher_exact([[a,n1-a],[b,n0-b]]); ps.append(pv); names.append(c)
q=multipletests(ps,method="fdr_bh")[1]
print("\nper-organ LUAD vs LUSC (single-organ patients):")
for c,pv,qv in zip(names,ps,q):
    a=t.loc["LUAD",c]; n1=t.loc["LUAD","n"]; b=t.loc["LUSC",c]; n0=t.loc["LUSC","n"]
    print(f"  {c:20s} LUAD {a:4d}/{n1} ({100*a/n1:4.1f}%) | LUSC {b:4d}/{n0} ({100*b/n0:4.1f}%) | p={pv:.3g} q={qv:.3g}")
t.to_csv(os.path.join(OUT,"stage2_single_organ_distribution.csv"))
pct.to_csv(os.path.join(OUT,"stage2_single_organ_pct.csv"))
