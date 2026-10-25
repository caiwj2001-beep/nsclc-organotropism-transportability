# -*- coding: utf-8 -*-
"""V2 numbers: LOO-based primary estimates, BH across the eight, attenuation series,
per-organ differences, missingness, and the full hit tables for the supplement."""
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2_contingency

PROJ = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
D2 = os.path.join(PROJ, "02_data_probe")
BASE = r"F:\数据库汇总\msk_chord_2024"
OUT = os.path.join(PROJ, "04_multivariable")
V2 = os.path.join(PROJ, "05_results", "v2")
os.makedirs(V2, exist_ok=True)

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

COVS = ["age", "male", "smk_fc", "stage4", "logtmb", "fga", "purity", "panel505"]
OCOL = {"BONE": "BONE_b", "LIVER": "LIVER_b", "CNS_BRAIN": "CNS_BRAIN_b", "PLEURA": "PLEURA_b",
        "ADRENAL_GLANDS": "ADRENAL_GLANDS_b", "INTRA_ABDOMINAL": "INTRA_ABDOMINAL_b",
        "REPRODUCTIVE_ORGANS": "REPRODUCTIVE_ORGANS_b", "OTHER": "OTHER_b"}
PAIRS = [("M_STK11", "ADRENAL_GLANDS"), ("M_STK11", "PLEURA"), ("M_STK11", "LIVER"),
         ("M_SMARCA4", "BONE"), ("M_SMARCA4", "ADRENAL_GLANDS"),
         ("M_EGFR", "ADRENAL_GLANDS"), ("A_FGFR1", "LIVER"), ("A_CDK4", "CNS_BRAIN")]


def fit(d, x, y, covs):
    cols = [x] + covs
    d = d.dropna(subset=[y] + cols)
    X = sm.add_constant(d[cols].astype(float))
    r = sm.Logit(d[y].astype(float), X).fit(disp=0, maxiter=200)
    b, se = r.params[x], r.bse[x]
    return dict(OR=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=r.pvalues[x], n=len(d))


# ---------- 1. primary LOO models + BH across the eight ----------
rows = []
for strata in ["LUAD", "LUSC"]:
    d0 = df[df.grp == strata].copy()
    for x, organ in PAIRS:
        y = OCOL[organ]
        d0["lg_loo"] = np.log1p(d0.n_distant - d0[y])
        d0["lg_orig"] = np.log(d0.n_distant)
        a = fit(d0, x, y, COVS + ["lg_loo"])
        b = fit(d0, x, y, COVS + ["lg_orig"])
        c = fit(d0, x, y, COVS)
        rows.append(dict(stratum=strata, exposure=x, organ=organ,
                         OR_loo=a["OR"], lo_loo=a["lo"], hi_loo=a["hi"], p_loo=a["p"], n=a["n"],
                         OR_logN=b["OR"], lo_logN=b["lo"], hi_logN=b["hi"],
                         OR_noburden=c["OR"], lo_nb=c["lo"], hi_nb=c["hi"],
                         carriers=int(d0[x].sum())))
res = pd.DataFrame(rows)
lu = res[res.stratum == "LUAD"].copy()
# BH across the eight pre-specified LUAD contrasts
o = np.argsort(lu.p_loo.values)
m = len(lu)
q = np.empty(m)
prev = 1.0
for rank in range(m, 0, -1):
    idx = o[rank - 1]
    val = min(prev, lu.p_loo.values[idx] * m / rank)
    q[idx] = val
    prev = val
lu["q_eight"] = q
res = res.merge(lu[["exposure", "organ", "q_eight"]], on=["exposure", "organ"], how="left")
res.round(3).to_csv(os.path.join(OUT, "stage4c_v2_primary.csv"), index=False)
print("=== V2 primary (LOO) LUAD, BH across the eight ===")
print(lu.assign(OR_loo=lu.OR_loo.round(2), lo=lu.lo_loo.round(2), hi=lu.hi_loo.round(2),
                OR_logN=lu.OR_logN.round(2), q=lu.q_eight.round(4))[
    ["exposure", "organ", "carriers", "OR_loo", "lo", "hi", "p_loo", "q", "OR_logN"]].to_string(index=False))

# ---------- 2. attenuation series for SMARCA4-bone with LOO burden ----------
d = df[df.grp == "LUAD"].copy()
y = "BONE_b"
d["lg_loo"] = np.log1p(d.n_distant - d[y])
c0 = [c for c in COVS if c != "logtmb"]
series = []
for lab, covs in [("Age, sex, smoking, stage,\nFGA, purity, panel\n(no burden, no TMB)", c0),
                  ("+ log(involved organs\nexcluding this organ)", c0 + ["lg_loo"]),
                  ("+ log(involved organs\nexcluding this organ)\n+ log(TMB)", c0 + ["lg_loo", "logtmb"])]:
    f = fit(d, "M_SMARCA4", y, covs)
    series.append(dict(model=lab, OR=f["OR"], lo=f["lo"], hi=f["hi"], n=f["n"]))
att = pd.DataFrame(series)
att.round(2).to_csv(os.path.join(V2, "fig3c_attenuation_loo.csv"), index=False)
print("\n=== SMARCA4-bone attenuation (LOO) ===")
print(att.assign(OR=att.OR.round(2), lo=att.lo.round(2), hi=att.hi.round(2)).to_string(index=False))

# ---------- 3. Fig 2 data (LOO-based) ----------
f2 = res[["stratum", "exposure", "organ", "OR_loo", "lo_loo", "hi_loo", "p_loo", "carriers"]].copy()
name = {("M_STK11", "ADRENAL_GLANDS"): "STK11 \u2013 adrenal", ("M_STK11", "PLEURA"): "STK11 \u2013 pleura",
        ("M_STK11", "LIVER"): "STK11 \u2013 liver", ("M_SMARCA4", "BONE"): "SMARCA4 \u2013 bone",
        ("M_SMARCA4", "ADRENAL_GLANDS"): "SMARCA4 \u2013 adrenal", ("M_EGFR", "ADRENAL_GLANDS"): "EGFR \u2013 adrenal",
        ("A_FGFR1", "LIVER"): "FGFR1 amp \u2013 liver", ("A_CDK4", "CNS_BRAIN"): "CDK4 amp \u2013 CNS"}
f2["label"] = [name[(e, o)] for e, o in zip(f2.exposure, f2.organ)]
it = pd.read_csv(os.path.join(OUT, "stage4b_loo_burden.csv"))  # not used; interaction recomputed below
pi = []
pool = df.copy()
for x, organ in PAIRS:
    yv = OCOL[organ]
    pool["lg_loo"] = np.log1p(pool.n_distant - pool[yv])
    dd = pool.dropna(subset=[yv, x, "hist"] + COVS + ["lg_loo"])
    dd["xhist"] = dd[x] * dd["hist"]
    X = sm.add_constant(dd[[x, "hist", "xhist"] + COVS + ["lg_loo"]].astype(float))
    r = sm.Logit(dd[yv].astype(float), X).fit(disp=0, maxiter=200)
    pi.append(dict(exposure=x, organ=organ, ratio=np.exp(r.params[x] + r.params["xhist"]),
                   p_int=r.pvalues["xhist"], n=len(dd)))
pit = pd.DataFrame(pi)
f2 = f2.merge(pit, on=["exposure", "organ"], how="left")
mde = pd.read_csv(os.path.join(OUT, "stage4_interaction_mde.csv"))
f2 = f2.merge(mde[["exposure", "organ", "MDE_ratio_alpha05"]], on=["exposure", "organ"], how="left")
f2.rename(columns={"OR_loo": "OR", "lo_loo": "lo", "hi_loo": "hi", "MDE_ratio_alpha05": "mde"}) \
  .to_csv(os.path.join(V2, "fig2_forest_loo.csv"), index=False)
print("\n=== Fig2 (LOO) + interaction + MDE ===")
print(f2.assign(OR=f2.OR_loo.round(2), mde=f2.MDE_ratio_alpha05.round(1), pint=pit.p_int.round(3))[
    ["label", "stratum", "OR_loo", "lo_loo", "hi_loo", "p_int", "MDE_ratio_alpha05"]].to_string(index=False))
print("\nMDE of the eight:", round(f2.MDE_ratio_alpha05.min(), 2), "to", round(f2.MDE_ratio_alpha05.max(), 2))

# ---------- 4. per-organ single-organ differences ----------
so = pd.read_csv(os.path.join(D2, "stage2_single_organ_pct.csv"), index_col=0)
d1 = df[(df.n_distant == 1)]
counts = {g: len(d1[d1.grp == g]) for g in ["LUAD", "LUSC"]}
print("\nsingle-organ counts:", counts)
keys7 = ["BONE", "LIVER", "CNS_BRAIN", "PLEURA", "ADRENAL_GLANDS", "INTRA_ABDOMINAL", "REPRODUCTIVE_ORGANS"]
tab = np.array([[so.loc["LUAD", k] / 100 * counts["LUAD"] for k in keys7],
                [so.loc["LUSC", k] / 100 * counts["LUSC"] for k in keys7]])
chi = chi2_contingency(tab)
print(f"7-organ chi2 = {chi[0]:.2f}, df={chi[2]}, p={chi[1]:.3f}")
tab8 = np.array([[so.loc["LUAD", k] / 100 * counts["LUAD"] for k in keys7 + ["OTHER"]],
                 [so.loc["LUSC", k] / 100 * counts["LUSC"] for k in keys7 + ["OTHER"]]])
chi8 = chi2_contingency(tab8)
print(f"8-category chi2 (incl OTHER) = {chi8[0]:.2f}, df={chi8[2]}, p={chi8[1]:.3f}")
per = []
for i, k in enumerate(keys7):
    a = tab[0, i] / counts["LUAD"]; b = tab[1, i] / counts["LUSC"]
    se = np.sqrt(a * (1 - a) / counts["LUAD"] + b * (1 - b) / counts["LUSC"])
    per.append(dict(organ=k, luad_pct=100 * a, lusc_pct=100 * b, diff_pp=100 * (a - b),
                    lo_pp=100 * (a - b - 1.96 * se), hi_pp=100 * (a - b + 1.96 * se),
                    luad_n=tab[0, i], lusc_n=tab[1, i]))
pd.DataFrame(per).round(2).to_csv(os.path.join(V2, "S4_single_organ_differences.csv"), index=False)
print(pd.DataFrame(per).round(1).to_string(index=False))

# ---------- 5. missingness by covariate and histology ----------
miss = []
for c in ["age", "male", "smk_fc", "stage4", "tmb", "fga", "purity", "GENE_PANEL"]:
    for g in ["LUAD", "LUSC"]:
        dd = df[df.grp == g]
        miss.append(dict(variable=c, histology=g,
                         missing_n=int(dd[c].isna().sum()), missing_pct=round(100 * dd[c].isna().mean(), 2),
                         n=len(dd)))
pd.DataFrame(miss).to_csv(os.path.join(V2, "S3_missingness.csv"), index=False)
print("\n=== complete-case retention ===")
for g in ["LUAD", "LUSC"]:
    dd = df[df.grp == g]
    cc = dd.dropna(subset=["age", "male", "smk_fc", "stage4", "logtmb", "fga", "purity", "panel505"])
    print(f"{g}: distant+ {len(dd)} -> complete case {len(cc)}  (lost {100*(1-len(cc)/len(dd)):.1f}%)")

# ---------- 6. full hit tables ----------
for f, lab in [("stage3_screen_mutation.csv", "S1_mutation_hits"),
               ("stage3_screen_cna_amp.csv", "S1_amp_hits"),
               ("stage3_screen_cna_deepdel.csv", "S1_del_hits")]:
    p = os.path.join(PROJ, "03_univariate", f)
    if os.path.exists(p):
        t = pd.read_csv(p)
        qc = [c for c in t.columns if c.startswith("q")]
        sig = t[t[qc[0]] < 0.05] if qc else t.head(0)
        sig.to_csv(os.path.join(V2, lab + ".csv"), index=False)
        print(f"{lab}: {len(t)} tests, {len(sig)} at first q < 0.05 (col {qc[0] if qc else 'n/a'})")
print("\nwritten to", V2)
