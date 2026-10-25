# -*- coding: utf-8 -*-
"""Prepare plotting tables for Figures 1-3 (all values traced to result CSVs)."""
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm

PROJ = r"E:\科研文件\肺癌器官特异性转移"
P2 = os.path.join(PROJ, "02_data_probe")
OUT = os.path.join(PROJ, "05_results", "figdata")
os.makedirs(OUT, exist_ok=True)

# ---------- Fig 1B: organ distribution (% of distant+ patients) ----------
org = pd.read_csv(os.path.join(P2, "stage2_organ_distribution.csv"), index_col=0)
ORGS = ["BONE", "CNS_BRAIN", "PLEURA", "LIVER", "ADRENAL_GLANDS"]
LAB = {"BONE": "Bone", "CNS_BRAIN": "CNS", "PLEURA": "Pleura",
       "LIVER": "Liver", "ADRENAL_GLANDS": "Adrenal"}
rows = []
for o in ORGS:
    rows.append({"organ": LAB[o], "LUAD": org.loc["LUAD", o + "_pct_distant"],
                 "LUSC": org.loc["LUSC", o + "_pct_distant"]})
pd.DataFrame(rows).to_csv(os.path.join(OUT, "fig1b_organ.csv"), index=False)

# ---------- Fig 1C: single-organ distribution ----------
so = pd.read_csv(os.path.join(P2, "stage2_single_organ_pct.csv"), index_col=0)
rows = []
for o in ORGS:
    rows.append({"organ": LAB[o], "LUAD": so.loc["LUAD", o], "LUSC": so.loc["LUSC", o]})
pd.DataFrame(rows).to_csv(os.path.join(OUT, "fig1c_single_organ.csv"), index=False)

# ---------- Fig 1D/1E: genomic landscape ----------
mut = pd.read_csv(os.path.join(P2, "stage2_mutation_landscape.csv"), index_col=0)
keep = ["TP53", "KRAS", "EGFR", "KMT2D", "CDKN2A", "STK11", "KEAP1", "RBM10",
        "SMARCA4", "PIK3CA", "FAT1", "ARID1A"]
mut.loc[[g for g in keep if g in mut.index], ["LUAD", "LUSC"]].reset_index() \
   .rename(columns={"index": "gene", "hugo": "gene"}) \
   .to_csv(os.path.join(OUT, "fig1d_mutation.csv"), index=False)

cna = pd.read_csv(os.path.join(P2, "stage2_cna_landscape_wide.csv"), index_col=0, header=[0, 1])
amp = cna["amp"] if "amp" in cna.columns.get_level_values(0) else None
if amp is not None:
    amp = amp.copy(); amp.columns = [c.replace("amp_", "") for c in amp.columns]
    amp.index.name = "gene"
    amp.reset_index().rename(columns=str).to_csv(os.path.join(OUT, "fig1e_cna.csv"), index=False)
    print("fig1e columns:", list(amp.columns))

# ---------- Fig 2: forest (LUAD vs LUSC + interaction) ----------
mv = pd.read_csv(os.path.join(PROJ, "04_multivariable", "stage4_multivariable.csv"))
it = pd.read_csv(os.path.join(PROJ, "04_multivariable", "stage4_interaction_mde.csv"))
lab = {("M_STK11", "ADRENAL_GLANDS"): "STK11 – adrenal",
       ("M_STK11", "PLEURA"): "STK11 – pleura",
       ("M_STK11", "LIVER"): "STK11 – liver",
       ("M_SMARCA4", "BONE"): "SMARCA4 – bone",
       ("M_SMARCA4", "ADRENAL_GLANDS"): "SMARCA4 – adrenal",
       ("M_EGFR", "ADRENAL_GLANDS"): "EGFR – adrenal",
       ("A_FGFR1", "LIVER"): "FGFR1 amp – liver",
       ("A_CDK4", "CNS_BRAIN"): "CDK4 amp – CNS"}
rows = []
for (exp, organ), name in lab.items():
    for stratum in ["LUAD", "LUSC"]:
        r = mv[(mv.exposure == exp) & (mv.organ == organ) & (mv.stratum == stratum)]
        if len(r):
            r = r.iloc[0]
            rows.append({"label": name, "stratum": stratum, "OR": r.OR,
                         "lo": r.lo, "hi": r.hi, "p": r.p, "carriers": int(r.carriers)})
f2 = pd.DataFrame(rows)
ii = it.set_index(["exposure", "organ"])
pint = {name: (float(ii.loc[(e, o), "p_interaction"]), float(ii.loc[(e, o), "MDE_ratio_alpha05"]))
        for (e, o), name in lab.items() if (e, o) in ii.index}
f2["p_int"] = f2["label"].map(lambda x: pint.get(x, (np.nan, np.nan))[0])
f2["mde"] = f2["label"].map(lambda x: pint.get(x, (np.nan, np.nan))[1])
f2.to_csv(os.path.join(OUT, "fig2_forest.csv"), index=False)
print(f2.round(3).to_string(index=False))

# ---------- Fig 3A: scenario ORs ----------
sc = pd.read_csv(os.path.join(PROJ, "06_validation", "stage6_sensitivity.csv"))
sc = sc[(sc.stratum == "LUAD") & sc.exposure.isin(["M_STK11", "M_SMARCA4"]) &
        sc.organ.isin(["ADRENAL_GLANDS"])]
sc["label"] = sc.exposure.map({"M_STK11": "STK11 – adrenal", "M_SMARCA4": "SMARCA4 – adrenal"})
sc["scenario_en"] = sc.scenario.map({"main_complete_case": "Complete case",
                                     "missing_indicator": "Missing indicator",
                                     "primary_only": "Primary sample",
                                     "single_organ": "Single organ"})
sc[["label", "scenario_en", "OR", "lo", "hi", "n", "carriers"]].to_csv(
    os.path.join(OUT, "fig3a_scenarios.csv"), index=False)

# ---------- Fig 3C: attenuation of SMARCA4 – bone (same covariate set, nested models) ----------
pat = pd.read_csv(os.path.join(P2, "stage2_patient_table.csv"), index_col=0, low_memory=False)
pat = pat[pat.grp == "LUAD"]
core = pd.read_csv(os.path.join(P2, "panel_core_genes.csv"))
core_e = set(core.entrezGeneId.astype(int))
mut = pd.read_csv(r"F:\数据库汇总\msk_chord_2024\data_mutations.csv", sep=",",
                  usecols=["sampleId", "patientId", "entrezGeneId"], low_memory=False)
mut = mut[mut.entrezGeneId.isin(core_e)].drop_duplicates(["sampleId", "entrezGeneId"])
smar = set(mut.loc[mut.entrezGeneId == 6597, "patientId"])
pat = pat[pat.any_distant].copy()
pat["SMARCA4"] = pat.index.isin(smar).astype(int)
pat["lgorg"] = np.log(pat.n_distant)
pat["smk_fc"] = np.where(pat.SMOKING_PREDICTIONS_3_CLASSES.astype(str).str.strip() == "Former/Current Smoker", 1,
                         np.where(pat.SMOKING_PREDICTIONS_3_CLASSES.astype(str).str.strip() == "Never", 0, np.nan))
pat["logtmb"] = np.log1p(pd.to_numeric(pat.tmb, errors="coerce"))
pat["panel505"] = (pat.GENE_PANEL == "IMPACT505").astype(int)
BASE = ["age", "male", "smk_fc", "stage4", "logtmb", "fga", "purity", "panel505"]
BASE0 = [c for c in BASE if c != "logtmb"]


def fit_logit(d, cols):
    d = d.dropna(subset=["BONE_b", "SMARCA4"] + cols)
    X = sm.add_constant(d[["SMARCA4"] + cols].astype(float))
    r = sm.Logit(d.BONE_b.astype(float), X).fit(disp=0)
    b, se = r.params["SMARCA4"], r.bse["SMARCA4"]
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), len(d)


models = [("Age, sex, smoking,\nstage, TMB, FGA, purity,\npanel (no burden)", [c for c in BASE]),
          ("+ log(involved organs)", BASE + ["lgorg"]),
          ("+ log(involved organs)\n+ log(TMB)", BASE + ["lgorg"])]
rows = []
for name, cols in models:
    o, lo, hi, n = fit_logit(pat, cols)
    rows.append({"model": name, "OR": o, "lo": lo, "hi": hi, "n": n})
att = pd.DataFrame(rows)
att.to_csv(os.path.join(OUT, "fig3c_attenuation.csv"), index=False)
print("\nattenuation (same covariate set):")
print(att.round(2).to_string(index=False))

# ---------- Fig 3B: co-occurrence ----------
def fisher_pair(a, b):
    from scipy.stats import fisher_exact
    tab = pd.crosstab(pat[a], pat[b]).reindex(index=[0, 1], columns=[0, 1]).fillna(0).values
    orr, p = fisher_exact(tab)
    return orr, p, int(tab[1, 1])


pat["KEAP1"] = pat.index.isin(set(mut.loc[mut.entrezGeneId == 9817, "patientId"])).astype(int)
rows = []
for g in ["KEAP1", "SMARCA4"]:
    orr, p, both = fisher_pair("STK11", g) if "STK11" in pat.columns else (np.nan, np.nan, 0)
    rows.append({"pair": f"STK11 – {g}", "OR": orr, "p": p, "both": both})
# add STK11 mutation column if missing
if "STK11" not in pat.columns:
    pat["STK11"] = pat.index.isin(set(mut.loc[mut.entrezGeneId == 6794, "patientId"])).astype(int)
rows = []
for g in ["KEAP1", "SMARCA4"]:
    orr, p, both = fisher_pair("STK11", g)
    rows.append({"pair": f"STK11 – {g}", "OR": orr, "p": p, "both": both})
pd.DataFrame(rows).to_csv(os.path.join(OUT, "fig3b_cooccurrence.csv"), index=False)
print("\nco-occurrence:")
print(pd.DataFrame(rows).round(2).to_string(index=False))
print("\nwrote figdata to", OUT)
