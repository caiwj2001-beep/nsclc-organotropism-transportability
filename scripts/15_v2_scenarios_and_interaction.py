# -*- coding: utf-8 -*-
"""V2 rebuild: LOO-based sensitivity scenarios, carrier counts in the model population,
interaction with and without a histology-by-burden term, and refreshed figure data."""
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm

PROJ = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
D2 = os.path.join(PROJ, "02_data_probe")
BASE = r"F:\数据库汇总\msk_chord_2024"
FD = os.path.join(PROJ, "05_results", "figdata")
V2 = os.path.join(PROJ, "05_results", "v2")

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
numc = [c for c in df.columns if c.startswith(("M_", "A_", "D_"))]
df[numc] = df[numc].fillna(0).astype(int)
df = df[df.any_distant].copy()
df["n_distant"] = df.n_distant.astype(int)
df["panel505"] = (df.GENE_PANEL == "IMPACT505").astype(int)
df["stage4"] = df.stage4.astype(int)
smk = df.SMOKING_PREDICTIONS_3_CLASSES.astype(str)
df["smk_fc"] = np.where(smk.str.strip() == "Former/Current Smoker", 1,
                        np.where(smk.str.strip() == "Never", 0, np.nan))
df["logtmb"] = np.log1p(pd.to_numeric(df.tmb, errors="coerce"))
df["hist"] = (df.grp == "LUSC").astype(int)
COVS = ["age", "male", "smk_fc", "stage4", "logtmb", "fga", "purity", "panel505"]
OCOL = {"BONE": "BONE_b", "LIVER": "LIVER_b", "CNS_BRAIN": "CNS_BRAIN_b", "PLEURA": "PLEURA_b",
        "ADRENAL_GLANDS": "ADRENAL_GLANDS_b"}


def fit(d, x, y, covs):
    d = d.dropna(subset=[y, x] + covs)
    X = sm.add_constant(d[[x] + covs].astype(float))
    r = sm.Logit(d[y].astype(float), X).fit(disp=0, maxiter=200)
    b, se = r.params[x], r.bse[x]
    return dict(OR=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=r.pvalues[x], n=len(d))


# ---------- carriers in the model (complete-case) population ----------
print("=== carriers in complete-case model population ===")
car = {}
for g in ["LUAD", "LUSC"]:
    d = df[df.grp == g].dropna(subset=["age", "male", "smk_fc", "stage4", "logtmb", "fga", "purity", "panel505"])
    car[g] = {x: int(d[x].sum()) for x in ["M_STK11", "M_SMARCA4", "M_EGFR", "A_FGFR1", "A_CDK4"]}
    print(f"  {g}: n={len(d)}  " + "  ".join(f"{k}={v}" for k, v in car[g].items()))

# ---------- four sensitivity scenarios with LOO burden ----------
scen = []
for x, organ in [("M_STK11", "ADRENAL_GLANDS"), ("M_SMARCA4", "ADRENAL_GLANDS")]:
    y = OCOL[organ]
    d = df[df.grp == "LUAD"].copy()
    d["lg_loo"] = np.log1p(d.n_distant - d[y])
    base = COVS + ["lg_loo"]
    f1 = fit(d, x, y, base)
    scen.append(dict(label=x.replace("M_", "") + " \u2013 adrenal", scenario_en="Complete case",
                     OR=f1["OR"], lo=f1["lo"], hi=f1["hi"], n=f1["n"], carriers=car["LUAD"][x]))
    # missing indicator: median fill + smoking unknown flag
    d2 = d.copy()
    for c in ["age", "fga", "purity", "logtmb"]:
        d2[c] = d2[c].fillna(d2[c].median())
    d2["smk_unk"] = d2.smk_fc.isna().astype(int)
    d2["smk_fc"] = d2.smk_fc.fillna(0)
    f2_ = fit(d2, x, y, base + ["smk_unk"])
    scen.append(dict(label=x.replace("M_", "") + " \u2013 adrenal", scenario_en="Missing indicator",
                     OR=f2_["OR"], lo=f2_["lo"], hi=f2_["hi"], n=f2_["n"], carriers=car["LUAD"][x]))
    # primary samples only
    d3 = d[d.SAMPLE_TYPE.astype(str).str.strip().str.lower() == "primary"].copy()
    f3 = fit(d3, x, y, base)
    scen.append(dict(label=x.replace("M_", "") + " \u2013 adrenal", scenario_en="Primary sample",
                     OR=f3["OR"], lo=f3["lo"], hi=f3["hi"], n=f3["n"], carriers=int(d3[x].sum())))
    # single-organ patients: burden-free, no burden term
    d4 = d[d.n_distant == 1].copy()
    f4 = fit(d4, x, y, COVS)
    scen.append(dict(label=x.replace("M_", "") + " \u2013 adrenal", scenario_en="Single organ",
                     OR=f4["OR"], lo=f4["lo"], hi=f4["hi"], n=f4["n"], carriers=int(d4[x].sum())))
sc = pd.DataFrame(scen)
sc.round(3).to_csv(os.path.join(FD, "fig3a_scenarios.csv"), index=False)
print("\n=== four scenarios (LOO burden) ===")
print(sc.assign(OR=sc.OR.round(2), lo=sc.lo.round(2), hi=sc.hi.round(2)).to_string(index=False))

# ---------- interaction: with and without histology-by-burden ----------
print("\n=== interaction p, with vs without histology x burden ===")
rows = []
for x, organ in [("M_STK11", "ADRENAL_GLANDS"), ("M_STK11", "PLEURA"), ("M_STK11", "LIVER"),
                 ("M_SMARCA4", "BONE"), ("M_SMARCA4", "ADRENAL_GLANDS"),
                 ("M_EGFR", "ADRENAL_GLANDS"), ("A_FGFR1", "LIVER"), ("A_CDK4", "CNS_BRAIN")]:
    y = OCOL[organ]
    d = df.copy()
    d["lg_loo"] = np.log1p(d.n_distant - d[y])
    d["xhist"] = d[x] * d["hist"]
    base = [x, "hist", "xhist"] + COVS + ["lg_loo"]
    dd = d.dropna(subset=[y, x] + COVS + ["lg_loo"])
    try:
        X = sm.add_constant(dd[base].astype(float))
        rA = sm.Logit(dd[y].astype(float), X).fit(disp=0, maxiter=200)
        pA = rA.pvalues["xhist"]
    except Exception:
        pA = np.nan
    try:
        dd2 = dd.copy()
        dd2["hist_burd"] = dd2["hist"] * dd2["lg_loo"]
        X2 = sm.add_constant(dd2[base + ["hist_burd"]].astype(float))
        rB = sm.Logit(dd2[y].astype(float), X2).fit(disp=0, maxiter=200)
        pB = rB.pvalues["xhist"]
    except Exception:
        pB = np.nan
    rows.append(dict(exposure=x, organ=organ, p_without=pA, p_with=pB, n=len(dd)))
    print(f"  {x:10s} {organ:15s} p={pA:.3f}  p(+hist x burden)={pB:.3f}  n={len(dd)}")
pd.DataFrame(rows).round(4).to_csv(os.path.join(V2, "interaction_with_hist_burden.csv"), index=False)
print("\nmin p without:", round(min(r['p_without'] for r in rows), 3),
      "| min p with:", round(min(r['p_with'] for r in rows), 3))

# ---------- refresh figure data ----------
f2 = pd.read_csv(os.path.join(V2, "fig2_forest_loo.csv"))
f2o = f2.rename(columns={"OR_loo": "OR", "lo_loo": "lo", "hi_loo": "hi", "p_loo": "p"})[
    ["label", "stratum", "OR", "lo", "hi", "p", "carriers", "p_int", "mde"]]
f2o.to_csv(os.path.join(FD, "fig2_forest.csv"), index=False)
a = pd.read_csv(os.path.join(V2, "fig3c_attenuation_loo.csv"))
a.to_csv(os.path.join(FD, "fig3c_attenuation.csv"), index=False)
print("\nfigdata refreshed: fig2_forest.csv (LOO), fig3a_scenarios.csv (LOO), fig3c_attenuation.csv (LOO)")
