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
"""STAGE 2A: analysis sets, baseline, organ distribution, mutation landscape,
METASTATIC_SITE vs patient-level organ flag concordance (kappa)."""
import os, json
import numpy as np
import pandas as pd

BASE = _CHORD
OUT = _os.path.abspath(_OUT)
os.makedirs(OUT, exist_ok=True)

pat = pd.read_csv(os.path.join(BASE, "data_clinical_patient.csv"), sep=",", low_memory=False)
smp = pd.read_csv(os.path.join(BASE, "data_clinical_sample.csv"), sep=",", low_memory=False)
mut = pd.read_csv(os.path.join(BASE, "data_mutations.csv"), sep=",",
                  usecols=["sampleId", "patientId", "entrezGeneId"], low_memory=False)
mp = mut.drop_duplicates("sampleId").set_index("sampleId")["patientId"]
smp["patientId"] = smp["sampleId"].map(mp)
core = pd.read_csv(os.path.join(OUT, "panel_core_genes.csv"))
core_e = set(core.entrezGeneId.astype(int))
gmap = pd.read_csv(_GENEMAP)
gmap.columns = [c.strip().lstrip("\ufeff") for c in gmap.columns]
gmap["entrezGeneId"] = pd.to_numeric(gmap.entrezGeneId, errors="coerce")
e2h = dict(zip(gmap.entrezGeneId.dropna().astype(int), gmap.hugoGeneSymbol))

# ---------- 1. analysis sets ----------
print("=" * 100)
print("1. ANALYSIS SETS (patient level, one histology label per patient)")
print("=" * 100)
_d = smp.dropna(subset=["patientId", "CANCER_TYPE_DETAILED"])
_d = _d.assign(_lung=_d.CANCER_TYPE_DETAILED.astype(str).str.contains("Lung|lung").astype(int))
_d = _d.sort_values(["patientId", "_lung"], ascending=[True, False])
det = _d.groupby("patientId")["CANCER_TYPE_DETAILED"].agg(lambda s: s.iloc[0])
pat = pat.set_index("patientId").join(det.rename("DETAILED"))
nsclc_ct = set(smp.loc[smp.CANCER_TYPE == "Non-Small Cell Lung Cancer", "patientId"].dropna())
pat["is_nsclc"] = pat.index.isin(nsclc_ct)
print("CHORD NSCLC patients:", int(pat.is_nsclc.sum()))
print("\nAll DETAILED labels inside CHORD-NSCLC (patient level):")
print(pat.loc[pat.is_nsclc, "DETAILED"].value_counts().to_string())

MAIN = {"Lung Adenocarcinoma": "LUAD", "Lung Squamous Cell Carcinoma": "LUSC"}
SECOND = {"Non-Small Cell Lung Cancer": "NSCLC-NOS",
          "Poorly Differentiated Non-Small Cell Lung Cancer": "NSCLC-PD"}
EXCLUDE_NE = ["Lung Carcinoid", "Atypical Lung Carcinoid", "Lung Neuroendocrine Tumor",
              "Large Cell Neuroendocrine Carcinoma"]
EXCLUDE_OTHER = ["Lung Adenosquamous Carcinoma", "Pleomorphic Carcinoma of the Lung",
                 "Sarcomatoid Carcinoma of the Lung", "Large Cell Lung Carcinoma",
                 "Adenoid Cystic Carcinoma of the Lung", "Mucoepidermoid Carcinoma of the Lung",
                 "NUT Carcinoma of the Lung", "Lymphoepithelioma-like Carcinoma of the Lung",
                 "Spindle Cell Carcinoma of the Lung", "Basaloid Large Cell Carcinoma of the Lung",
                 "Lung Adenocarcinoma In Situ"]
pat["grp"] = None
pat.loc[pat.DETAILED.isin(MAIN), "grp"] = pat.loc[pat.DETAILED.isin(MAIN), "DETAILED"].map(MAIN)
pat.loc[pat.DETAILED.isin(SECOND), "grp"] = pat.loc[pat.DETAILED.isin(SECOND), "DETAILED"].map(SECOND)
pat.loc[pat.DETAILED.isin(EXCLUDE_NE), "grp"] = "EXCL-neuroendocrine"
pat.loc[pat.DETAILED.isin(EXCLUDE_OTHER), "grp"] = "EXCL-other-histology"
pat.loc[pat.is_nsclc & pat.grp.isna(), "grp"] = "EXCL-unlabelled"
print("\ngroup assignment:")
print(pat.grp.value_counts(dropna=False).to_string())

# ---------- 2. phenotype + baseline ----------
DISTANT = ["BONE", "LIVER", "CNS_BRAIN", "PLEURA", "ADRENAL_GLANDS",
           "INTRA_ABDOMINAL", "REPRODUCTIVE_ORGANS", "OTHER"]
for c in DISTANT + ["LYMPH_NODES", "LUNG"]:
    pat[c + "_b"] = (pat[c].astype(str).str.strip().str.lower() == "yes")
pat["n_distant"] = pat[[c + "_b" for c in DISTANT]].sum(axis=1)
pat["any_distant"] = pat.n_distant > 0
pat["age"] = pd.to_numeric(pat.CURRENT_AGE_DEID, errors="coerce")
pat["male"] = (pat.GENDER.astype(str).str.lower() == "male").astype(int)

# chosen sample per patient (Primary preferred)
prio = {"Primary": 0, "Local Recurrence": 1, "Metastasis": 2, "Unknown": 3, "nan": 4}
s = smp[smp.patientId.notna()].copy()
s["prio"] = s.SAMPLE_TYPE.astype(str).map(prio).fillna(4)
s = s.sort_values(["patientId", "prio", "sampleId"])
ch = s.drop_duplicates("patientId", keep="first").set_index("patientId")
pat = pat.join(ch[["sampleId", "SAMPLE_TYPE", "GENE_PANEL", "TMB_NONSYNONYMOUS",
                   "FRACTION_GENOME_ALTERED", "TUMOR_PURITY", "MSI_TYPE"]])
pat["tmb"] = pd.to_numeric(pat.TMB_NONSYNONYMOUS, errors="coerce")
pat["fga"] = pd.to_numeric(pat.FRACTION_GENOME_ALTERED, errors="coerce")
pat["purity"] = pd.to_numeric(pat.TUMOR_PURITY, errors="coerce")
pat["stage4"] = pat.STAGE_HIGHEST_RECORDED.astype(str).str.contains("4").astype(int)

print("\n" + "=" * 100)
print("2. ORGAN DISTRIBUTION (patient level)")
print("=" * 100)
rows = []
for g in ["LUAD", "LUSC", "NSCLC-NOS", "NSCLC-PD"]:
    d = pat[pat.grp == g]
    if not len(d):
        continue
    dm = d[d.any_distant]
    r = {"group": g, "n": len(d), "any_distant": int(d.any_distant.sum()),
         "pct_distant": round(100 * d.any_distant.mean(), 1)}
    for c in DISTANT + ["LYMPH_NODES"]:
        r[c] = int(d[c + "_b"].sum())
        r[c + "_pct_all"] = round(100 * d[c + "_b"].mean(), 1)
        r[c + "_pct_distant"] = round(100 * dm[c + "_b"].mean(), 1) if len(dm) else None
    r["n_org_mean"] = round(dm.n_distant.mean(), 2) if len(dm) else None
    r["solitary_single_organ"] = int((dm.n_distant == 1).sum())
    rows.append(r)
organ = pd.DataFrame(rows).set_index("group")
print(organ[["n", "any_distant", "pct_distant", "BONE", "CNS_BRAIN", "PLEURA", "LIVER",
             "ADRENAL_GLANDS", "LYMPH_NODES", "n_org_mean", "solitary_single_organ"]].to_string())
print("\nBONE pct of distant+:", organ["BONE_pct_distant"].to_dict())
organ.to_csv(os.path.join(OUT, "stage2_organ_distribution.csv"))

print("\n" + "=" * 100)
print("3. BASELINE (LUAD vs LUSC vs NSCLC-NOS)")
print("=" * 100)
base = []
for g in ["LUAD", "LUSC", "NSCLC-NOS"]:
    d = pat[pat.grp == g]
    base.append({"group": g, "n": len(d),
                 "age_median": round(d.age.median(), 1),
                 "male_pct": round(100 * d.male.mean(), 1),
                 "stage4_pct": round(100 * d.stage4.mean(), 1),
                 "tmb_median": round(d.tmb.median(), 2),
                 "fga_median": round(d.fga.median(), 3),
                 "purity_median": round(d.purity.median(), 3),
                 "primary_sample_pct": round(100 * (d.SAMPLE_TYPE == "Primary").mean(), 1),
                 "panel_505_pct": round(100 * (d.GENE_PANEL == "IMPACT505").mean(), 1),
                 "smoke_never_pct": round(100 * (d.SMOKING_PREDICTIONS_3_CLASSES.astype(str)
                                                 .str.lower().str.startswith("never")).mean(), 1)})
bt = pd.DataFrame(base).set_index("group")
print(bt.to_string())
bt.to_csv(os.path.join(OUT, "stage2_baseline.csv"))

# ---------- 4. mutation landscape ----------
print("\n" + "=" * 100)
print("4. MUTATION LANDSCAPE (carrier %, chosen sample)")
print("=" * 100)
samp2grp = pat.dropna(subset=["sampleId"]).set_index("sampleId")["grp"].to_dict()
m = mut[mut.sampleId.isin(samp2grp)].copy()
m["grp"] = m.sampleId.map(samp2grp)
m = m.drop_duplicates(["sampleId", "entrezGeneId"])
m["hugo"] = m.entrezGeneId.map(e2h)
land = []
for g in ["LUAD", "LUSC", "NSCLC-NOS"]:
    d = m[(m.grp == g)]
    n = len(pat[pat.grp == g])
    vc = (d.groupby("hugo").sampleId.nunique() / n * 100).sort_values(ascending=False)
    land.append(vc.rename(g))
lt = pd.concat(land, axis=1).fillna(0).round(1)
core_h_set = set(core.hugoGeneSymbol.dropna())
lt["in_core335"] = [k in core_h_set for k in lt.index]
print("\nTOP 25 all genes (any panel):")
print(lt.sort_values("LUAD", ascending=False).head(25).to_string())
print("\nTOP 20 within 335-gene core set (comparable across panels):")
core_h = set(core.hugoGeneSymbol.dropna())
print(lt[lt.index.isin(core_h)].sort_values("LUAD", ascending=False).head(20).to_string())
lt.to_csv(os.path.join(OUT, "stage2_mutation_landscape.csv"))

# ---------- 5. METASTATIC_SITE vs organ flag concordance ----------
print("\n" + "=" * 100)
print("5. SAMPLE-LEVEL METASTATIC_SITE vs PATIENT-LEVEL ORGAN FLAG")
print("=" * 100)
MAP = {"liver": "LIVER", "bone": "BONE", "iliac bone": "BONE", "spine": "BONE",
       "brain": "CNS_BRAIN", "pleura": "PLEURA", "pleural fluid": "PLEURA",
       "adrenal gland": "ADRENAL_GLANDS", "lymph node": "LYMPH_NODES", "lung": "LUNG",
       "peritoneum": "INTRA_ABDOMINAL", "omentum": "INTRA_ABDOMINAL",
       "retroperitoneum": "INTRA_ABDOMINAL", "ovary": "REPRODUCTIVE_ORGANS"}
ss = smp[smp.patientId.notna()].copy()
ss["site"] = ss.METASTATIC_SITE.astype(str).str.strip().str.lower()
ss["organ"] = ss.site.map(MAP)
iso = ss[ss.organ.notna()].copy()
iso = iso.join(pat[[c + "_b" for c in set(MAP.values())]], on="patientId")
print("samples with mappable METASTATIC_SITE:", len(iso),
      "| distinct organs:", iso.organ.nunique())
rows = []
for org in sorted(iso.organ.unique()):
    d = iso[iso.organ == org]
    col = org + "_b"
    if col not in d.columns:
        continue
    flag = d[col].astype(int).values
    # 2x2 with the complementary class: samples annotated organ X vs all other annotated sites
    other = iso[iso.organ != org]
    yy = np.concatenate([flag, other[col].astype(int).values])
    xx = np.concatenate([np.ones(len(d), int), np.zeros(len(other), int)])
    rows.append({"organ": org, "samples_annotated": len(d),
                 "flag_yes_pct": round(100 * flag.mean(), 1),
                 "flag_yes_pct_in_other_sites": round(100 * other[col].mean(), 1),
                 "risk_ratio": round((flag.mean() + 1e-9) / (other[col].mean() + 1e-9), 2)})
kt = pd.DataFrame(rows).sort_values("samples_annotated", ascending=False)
print(kt.to_string(index=False))
kt.to_csv(os.path.join(OUT, "stage2_site_flag_concordance.csv"), index=False)

pat.to_csv(os.path.join(OUT, "stage2_patient_table.csv"), index_label="patientId")
print("\nwrote stage2_patient_table.csv", pat.shape)
