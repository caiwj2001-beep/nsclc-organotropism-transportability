# -*- coding: utf-8 -*-
"""V2 decisive analysis: leave-one-out burden covariate.
Burden N sums the eight organ indicators, one of which is the outcome Y_j.
Refit every principal association with log(1 + N - Y_j) and compare to log(N) and to no burden.
"""
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm

PROJ = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
D2 = os.path.join(PROJ, "02_data_probe")
BASE = r"F:\数据库汇总\msk_chord_2024"
OUT = os.path.join(PROJ, "04_multivariable")

pat = pd.read_csv(os.path.join(D2, "stage2_patient_table.csv"), index_col=0, low_memory=False)
pat = pat[pat.grp.isin(["LUAD", "LUSC"])]
core = pd.read_csv(os.path.join(D2, "panel_core_genes.csv"))
core_e = set(core.entrezGeneId.astype(int))
h2e = dict(zip(core.entrezGeneId.astype(int), core.hugoGeneSymbol))

mut = pd.read_csv(os.path.join(BASE, "data_mutations.csv"), sep=",",
                  usecols=["sampleId", "entrezGeneId"], low_memory=False)
s2p = dict(zip(pat.sampleId, pat.index))
mut = mut[mut.sampleId.isin(s2p) & mut.entrezGeneId.isin(core_e)].drop_duplicates()
mut["patientId"] = mut.sampleId.map(s2p)
M = mut.assign(v=1).pivot_table(index="patientId", columns="entrezGeneId", values="v", aggfunc="max")
M.columns = ["M_" + h2e.get(c, str(c)) for c in M.columns]

AMP = pd.read_csv(os.path.join(D2, "stage2_cna_amp.csv.gz"), index_col=0)
DEL = pd.read_csv(os.path.join(D2, "stage2_cna_deepdel.csv.gz"), index_col=0)
AMP.columns = ["A_" + c for c in AMP.columns]
DEL.columns = ["D_" + c for c in DEL.columns]

df = pat.join([M, AMP, DEL], how="inner").copy()
numcols = [c for c in df.columns if c.startswith(("M_", "A_", "D_"))]
df[numcols] = df[numcols].fillna(0).astype(int)
df = df[df.any_distant].copy()
df["n_distant"] = df.n_distant.astype(int)
df["panel505"] = (df.GENE_PANEL == "IMPACT505").astype(int)
df["stage4"] = df.stage4.astype(int)
smk = df.SMOKING_PREDICTIONS_3_CLASSES.astype(str)
df["smk_fc"] = np.where(smk.str.strip() == "Former/Current Smoker", 1,
                        np.where(smk.str.strip() == "Never", 0, np.nan))
df["logtmb"] = np.log1p(pd.to_numeric(df.tmb, errors="coerce"))
df["hist"] = (df.grp == "LUSC").astype(int)
df["lgorg"] = np.log(df.n_distant)
print("frame distant+:", len(df), "| LUAD", int((df.grp == "LUAD").sum()), "| LUSC", int((df.grp == "LUSC").sum()))

COVS = ["age", "male", "smk_fc", "stage4", "logtmb", "fga", "purity", "panel505"]
PAIRS = [("M_STK11", "ADRENAL_GLANDS"), ("M_STK11", "PLEURA"), ("M_STK11", "LIVER"),
         ("M_SMARCA4", "BONE"), ("M_SMARCA4", "ADRENAL_GLANDS"),
         ("M_EGFR", "ADRENAL_GLANDS"), ("A_FGFR1", "LIVER"), ("A_CDK4", "CNS_BRAIN")]
ORGAN_COL = {"BONE": "BONE_b", "LIVER": "LIVER_b", "CNS_BRAIN": "CNS_BRAIN_b",
             "PLEURA": "PLEURA_b", "ADRENAL_GLANDS": "ADRENAL_GLANDS_b",
             "INTRA_ABDOMINAL": "INTRA_ABDOMINAL_b", "REPRODUCTIVE_ORGANS": "REPRODUCTIVE_ORGANS_b",
             "OTHER": "OTHER_b"}


def fit_one(d, x, y, burden_col):
    cols = [x] + COVS + ([burden_col] if burden_col else [])
    d = d.dropna(subset=[y, x] + COVS + ([burden_col] if burden_col else []))
    X = sm.add_constant(d[cols].astype(float))
    try:
        r = sm.Logit(d[y].astype(float), X).fit(disp=0, maxiter=200)
    except Exception as e:
        return None, len(d), str(e)
    b, se = r.params[x], r.bse[x]
    return (np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), r.pvalues[x]), len(d), None


rows = []
for strata in ["LUAD", "LUSC"]:
    d0 = df[df.grp == strata].copy()
    for x, organ in PAIRS:
        y = ORGAN_COL[organ]
        d0["lgorg_loo"] = np.log1p(d0.n_distant - d0[y])
        d0["lgorg_orig"] = np.log(d0.n_distant)
        o_orig, n1, e1 = fit_one(d0, x, y, "lgorg_orig")
        o_loo, n2, e2 = fit_one(d0, x, y, "lgorg_loo")
        o_non, n3, e3 = fit_one(d0, x, y, None)
        rows.append({
            "stratum": strata, "exposure": x, "organ": organ,
            "OR_orig_logN": None if o_orig is None else round(o_orig[0], 2),
            "lo_orig": None if o_orig is None else round(o_orig[1], 2),
            "hi_orig": None if o_orig is None else round(o_orig[2], 2),
            "OR_loo": None if o_loo is None else round(o_loo[0], 2),
            "lo_loo": None if o_loo is None else round(o_loo[1], 2),
            "hi_loo": None if o_loo is None else round(o_loo[2], 2),
            "p_loo": None if o_loo is None else round(o_loo[3], 4),
            "OR_noburden": None if o_non is None else round(o_non[0], 2),
            "lo_nb": None if o_non is None else round(o_non[1], 2),
            "hi_nb": None if o_non is None else round(o_non[2], 2),
            "n": max(n1, n2), "err": e1 or e2 or e3,
        })

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "stage4b_loo_burden.csv"), index=False)
pd.set_option("display.width", 200)
print("\n=== LUAD only ===")
print(res[res.stratum == "LUAD"][["exposure", "organ", "OR_noburden", "OR_orig_logN", "lo_orig", "hi_orig",
                                  "OR_loo", "lo_loo", "hi_loo", "p_loo", "n"]].to_string(index=False))
print("\n=== LUSC only (carrier-limited) ===")
print(res[res.stratum == "LUSC"][["exposure", "organ", "OR_orig_logN", "OR_loo", "lo_loo", "hi_loo", "n"]].to_string(index=False))

# pooled interaction with LOO burden
print("\n=== pooled interaction, LOO burden (LUAD vs LUSC) ===")
pool = df.copy()
for x, organ in PAIRS:
    y = ORGAN_COL[organ]
    pool["lgorg_loo"] = np.log1p(pool.n_distant - pool[y])
    d = pool.dropna(subset=[y, x, "hist"] + COVS + ["lgorg_loo"])
    d["xhist"] = d[x] * d["hist"]
    X = sm.add_constant(d[[x, "hist", "xhist"] + COVS + ["lgorg_loo"]].astype(float))
    try:
        r = sm.Logit(d[y].astype(float), X).fit(disp=0, maxiter=200)
        ratio = np.exp(r.params[x] + r.params["xhist"])
        se = r.bse[x] + 0  # interaction SE for the ratio-difference test
        print(f"{x:10s} {organ:15s} ratio={ratio:5.2f}  p_interaction={r.pvalues['xhist']:.3f}  n={len(d)}")
    except Exception as e:
        print(f"{x:10s} {organ:15s} failed: {e}")
