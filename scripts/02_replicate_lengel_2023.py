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
"""STEP B: LUAD replication of Lengel 2023 (Cancer Cell, PMID 37084736).
Estimand-matched part = cross-sectional logistic models (ever-metastatic vs non-metastatic;
site-specific odds within ever-metastatic). Time-to-metastasis HRs are NOT replicable
(no metastasis date in MSK-CHORD) and are reported as estimand-mismatched.

Published targets extracted verbatim from PMC10391526 full text:
  MV logistic, ever-metastatic (n=448): SMARCA4 -> BONE  OR 6.47 (1.70-42.60), p=0.017
                                       TP53    -> LYMPH OR 2.03 (1.28-3.23),  p=0.003
  ever-met (n=448) vs non-met (n=318): TP53, KEAP1, CDKN2A, MDM2, PIK3CA, NKX2-1, RB1,
                                       MYC, SMARCA4, FOXA1 higher; RBM10 lower (16% vs 7%)
"""
import json, os
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

BASE = _CHORD
OUT = _os.path.abspath(_OUT)
os.makedirs(OUT, exist_ok=True)

GENES = {7157: "TP53", 9817: "KEAP1", 1029: "CDKN2A", 4193: "MDM2", 5290: "PIK3CA",
         7080: "NKX2-1", 5925: "RB1", 4609: "MYC", 6597: "SMARCA4", 3169: "FOXA1",
         8241: "RBM10", 6794: "STK11", 4780: "NFE2L2", 3845: "KRAS", 1956: "EGFR"}
DISTANT = ["BONE", "LIVER", "CNS_BRAIN", "PLEURA", "ADRENAL_GLANDS",
           "INTRA_ABDOMINAL", "REPRODUCTIVE_ORGANS", "OTHER"]

pat = pd.read_csv(os.path.join(BASE, "data_clinical_patient.csv"), sep=",", low_memory=False)
smp = pd.read_csv(os.path.join(BASE, "data_clinical_sample.csv"), sep=",", low_memory=False)
mut = pd.read_csv(os.path.join(BASE, "data_mutations.csv"), sep=",",
                  usecols=["sampleId", "patientId", "entrezGeneId"], low_memory=False)
mp = mut.drop_duplicates("sampleId").set_index("sampleId")["patientId"]
smp["patientId"] = smp["sampleId"].map(mp)

# ---------- choose ONE sample per LUAD patient: Primary > Local Recurrence > Metastasis ----------
luad_ids = set(smp.loc[smp.CANCER_TYPE_DETAILED == "Lung Adenocarcinoma", "patientId"].dropna())
s = smp[smp.patientId.isin(luad_ids)].copy()
prio = {"Primary": 0, "Local Recurrence": 1, "Metastasis": 2, "Unknown": 3}
s["prio"] = s["SAMPLE_TYPE"].map(prio).fillna(4)
s = s.sort_values(["patientId", "prio", "sampleId"])
chosen = s.drop_duplicates("patientId", keep="first").set_index("patientId")
print("LUAD patients:", len(luad_ids), "| chosen samples:", len(chosen))
print("chosen SAMPLE_TYPE:", chosen["SAMPLE_TYPE"].value_counts().to_dict())
print("chosen GENE_PANEL:", chosen["GENE_PANEL"].value_counts().to_dict())

# ---------- patient-level covariate table ----------
df = pat[pat.patientId.isin(chosen.index)].set_index("patientId").copy()
df = df.join(chosen[["sampleId", "SAMPLE_TYPE", "GENE_PANEL", "TMB_NONSYNONYMOUS",
                     "FRACTION_GENOME_ALTERED", "TUMOR_PURITY"]])
for c in DISTANT:
    df[c + "_b"] = (df[c].astype(str).str.strip().str.lower() == "yes")
df["n_distant"] = df[[c + "_b" for c in DISTANT]].sum(axis=1)
df["ever_met"] = df["n_distant"] > 0
df["age"] = pd.to_numeric(df["CURRENT_AGE_DEID"], errors="coerce")
df["male"] = (df["GENDER"].astype(str).str.lower() == "male").astype(int)
df["smoke"] = df["SMOKING_PREDICTIONS_3_CLASSES"].astype(str).replace("nan", np.nan)
df["stage4"] = (df["STAGE_HIGHEST_RECORDED"].astype(str).str.contains("4")).astype(int)
df["tmb"] = pd.to_numeric(df["TMB_NONSYNONYMOUS"], errors="coerce")
df["fga"] = pd.to_numeric(df["FRACTION_GENOME_ALTERED"], errors="coerce")
df["purity"] = pd.to_numeric(df["TUMOR_PURITY"], errors="coerce")
print("\never-met:", int(df.ever_met.sum()), "| non-met:", int((~df.ever_met).sum()))
print("stage4 among ever-met:", df.loc[df.ever_met, "stage4"].mean().round(3))

# ---------- mutation matrix (chosen sample only) ----------
ids = set(GENES)
m = mut[mut.sampleId.isin(chosen["sampleId"]) & mut.entrezGeneId.isin(ids)]
m = m.drop_duplicates(["sampleId", "entrezGeneId"])
mm = m.assign(g=m.entrezGeneId.map(GENES)).pivot(index="sampleId", columns="g", values="entrezGeneId")
mm = mm.notna().astype(int).reindex(chosen["sampleId"].values).fillna(0).astype(int)
mm.index = mm.index.map(dict(zip(chosen["sampleId"], chosen.index)))
print("\nmut carrier count (LUAD, chosen sample):")
print(mm.sum().sort_values(ascending=False).to_string())

# ---------- CNA (chosen sample only, chunked) ----------
rows = []
reader = pd.read_csv(os.path.join(BASE, "data_cna_long.csv"), sep=",",
                     usecols=["sampleId", "entrezGeneId", "value"], chunksize=2_000_000, low_memory=False)
samp_set = set(chosen["sampleId"])
for ch in reader:
    ch = ch[ch.sampleId.isin(samp_set) & ch.entrezGeneId.isin(ids)]
    if len(ch):
        rows.append(ch)
cna = pd.concat(rows) if rows else pd.DataFrame(columns=["sampleId", "entrezGeneId", "value"])
cna["g"] = cna.entrezGeneId.map(GENES)
cna_deep = (cna.assign(v=(cna.value.astype(float) <= -1.5).astype(int))
            .pivot_table(index="sampleId", columns="g", values="v", aggfunc="max"))
cna_amp = (cna.assign(v=(cna.value.astype(float) == 2).astype(int))
           .pivot_table(index="sampleId", columns="g", values="v", aggfunc="max"))
cna_deep = cna_deep.reindex(chosen["sampleId"].values).fillna(0).astype(int)
cna_amp = cna_amp.reindex(chosen["sampleId"].values).fillna(0).astype(int)
cna_deep.index = cna_amp.index = chosen["sampleId"].values

# alteration = mutation OR deep deletion OR amplification
alt = pd.DataFrame(index=mm.index)
for g in mm.columns:
    d = cna_deep[g].values if g in cna_deep.columns else np.zeros(len(mm))
    a = cna_amp[g].values if g in cna_amp.columns else np.zeros(len(mm))
    alt[g] = np.maximum.reduce([mm[g].values, d, a])
print("\nany-alteration carrier counts (mutation OR deepdel OR amp):")
print(alt.sum().sort_values(ascending=False).to_string())

df = df.join(mm.add_prefix("mut_")).join(alt.add_prefix("alt_"))
df["cna_cdkn2a_deepdel"] = [cna_deep.loc[chosen.loc[p, "sampleId"], "CDKN2A"]
                            if "CDKN2A" in cna_deep.columns else 0 for p in df.index]


def fit(y, xgene, covs, data):
    d = data[[y, xgene] + covs].dropna()
    if d[y].nunique() < 2 or d[xgene].nunique() < 2:
        return None
    X = sm.add_constant(d[[xgene] + covs].astype(float))
    try:
        r = sm.Logit(d[y].astype(float), X).fit(disp=0, maxiter=200)
    except Exception as e:
        return None
    b, se = r.params[xgene], r.bse[xgene]
    return {"OR": float(np.exp(b)), "lo": float(np.exp(b - 1.96 * se)),
            "hi": float(np.exp(b + 1.96 * se)), "p": float(r.pvalues[xgene]),
            "n": int(len(d)), "events": int(d[y].sum()),
            "carrier_cases": int(d.loc[d[y] == 1, xgene].sum()),
            "carrier_ctrl": int(d.loc[d[y] == 0, xgene].sum())}


def fisher_row(y, xgene, data):
    d = data[[y, xgene]].dropna()
    a = int(((d[y] == 1) & (d[xgene] == 1)).sum()); b1 = int(((d[y] == 1) & (d[xgene] == 0)).sum())
    c = int(((d[y] == 0) & (d[xgene] == 1)).sum()); d1 = int(((d[y] == 0) & (d[xgene] == 0)).sum())
    from scipy.stats import fisher_exact
    orr, p = fisher_exact([[a, b1], [c, d1]])
    return {"OR_crude": orr, "p_crude": p, "n": len(d), "table": [a, b1, c, d1]}


COVS = ["age", "male", "tmb", "fga", "purity", "stage4"]
COVS_NOSTAGE = ["age", "male", "tmb", "fga", "purity"]

res = {}
# ---- (1) ever-met vs non-met (mirror Lengel NM vs ever-met; do NOT adjust stage here) ----
print("\n" + "=" * 100)
print("(1) EVER-METASTATIC vs NON-METASTATIC (LUAD, primary-preferred sample)")
print("=" * 100)
rows1 = []
for g in [v for v in GENES.values()]:
    col = "alt_" + g
    if col not in df.columns:
        continue
    f = fisher_row("ever_met", col, df)
    mv = fit("ever_met", col, COVS_NOSTAGE, df)
    if mv:
        rows1.append({"gene": g, "rate_evermet": f"{df.loc[df.ever_met,col].mean():.3f}",
                      "rate_nonmet": f"{df.loc[~df.ever_met,col].mean():.3f}",
                      "OR_crude": round(f["OR_crude"], 2), "p_crude": f["p_crude"],
                      "OR_adj": round(mv["OR"], 2), "lo": round(mv["lo"], 2), "hi": round(mv["hi"], 2),
                      "p_adj": mv["p"], "n": mv["n"]})
t1 = pd.DataFrame(rows1)
if len(t1):
    t1["q_adj"] = multipletests(t1.p_adj.fillna(1), method="fdr_bh")[1]
    t1 = t1.sort_values("p_adj")
    print(t1.to_string(index=False))
    t1.to_csv(os.path.join(OUT, "replication_step1_evermet_vs_nonmet.csv"), index=False)

# ---- (2) site-specific within ever-metastatic ----
print("\n" + "=" * 100)
print("(2) SITE-SPECIFIC ODDS WITHIN EVER-METASTATIC (mirror Lengel Fig3F)")
print("=" * 100)
em = df[df.ever_met].copy()
rows2 = []
for organ in ["BONE", "LIVER", "CNS_BRAIN", "PLEURA", "ADRENAL_GLANDS", "LYMPH_NODES"]:
    y = organ + "_b"
    if y not in em.columns:
        continue
    for g in ["SMARCA4", "TP53", "CDKN2A", "KEAP1", "STK11", "NFE2L2", "KRAS", "EGFR"]:
        col = "alt_" + g
        if col not in em.columns:
            continue
        f = fisher_row(y, col, em)
        mv = fit(y, col, COVS, em)
        if mv:
            rows2.append({"organ": organ, "gene": g,
                          "carrier_rate_in": round(em.loc[em[y], col].mean(), 3),
                          "carrier_rate_out": round(em.loc[~em[y], col].mean(), 3),
                          "OR_crude": round(f["OR_crude"], 2), "p_crude": f["p_crude"],
                          "OR_adj": round(mv["OR"], 2), "lo": round(mv["lo"], 2),
                          "hi": round(mv["hi"], 2), "p_adj": mv["p"],
                          "n": mv["n"], "cases": mv["events"]})
t2 = pd.DataFrame(rows2)
if len(t2):
    t2["q_adj"] = multipletests(t2.p_adj.fillna(1), method="fdr_bh")[1]
    t2 = t2.sort_values("p_adj")
    print(t2.head(30).to_string(index=False))
    t2.to_csv(os.path.join(OUT, "replication_step2_site_specific.csv"), index=False)

# ---- (3) CDKN2A deep deletion: metastatic lesions vs primary samples (Lengel 24% vs 14%) ----
print("\n" + "=" * 100)
print("(3) CDKN2A deep deletion rate: Primary vs Metastasis samples (LUAD)")
print("=" * 100)
luad_samples = smp[(smp.CANCER_TYPE_DETAILED == "Lung Adenocarcinoma") & smp.patientId.notna()]
cna_cd = cna[cna.g == "CDKN2A"] if "g" in cna.columns else pd.DataFrame()
if len(cna_cd):
    deep = cna_cd.assign(v=(cna_cd.value.astype(float) <= -1.5).astype(int))
    deep = deep.merge(luad_samples[["sampleId", "SAMPLE_TYPE"]], on="sampleId", how="inner")
    tab = deep.groupby("SAMPLE_TYPE").v.agg(["size", "sum", "mean"]).sort_values("size", ascending=False)
    print(tab.to_string())
    pr = tab.loc["Primary"] if "Primary" in tab.index else None
    me = tab.loc["Metastasis"] if "Metastasis" in tab.index else None
    if pr is not None and me is not None:
        from scipy.stats import fisher_exact
        a, n1 = int(me["sum"]), int(me["size"]); b1, n0 = int(pr["sum"]), int(pr["size"])
        orr, p = fisher_exact([[a, n1 - a], [b1, n0 - b1]])
        print(f"\nMetastasis {a}/{n1} ({100*a/n1:.1f}%) vs Primary {b1}/{n0} ({100*b1/n0:.1f}%) "
              f"OR={orr:.2f} p={p:.3g}")

df.to_csv(os.path.join(OUT, "luad_analysis_table.csv"), index_label="patientId")
print("\nwrote luad_analysis_table.csv", df.shape)
