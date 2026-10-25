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
"""STEP A: empirical panel intersection gene set (341/410/468/505) + carrier counts.
Approximation documented: a gene's presence on a panel is inferred from being observed
mutated in >=1 sample run on that panel (observed-set approximation). Intersection across
the four panels is therefore a CONSERVATIVE core set (subset of the true common panel)."""
import json, os
import numpy as np
import pandas as pd

BASE = _CHORD
OUT = _os.path.abspath(_OUT)
os.makedirs(OUT, exist_ok=True)

smp = pd.read_csv(os.path.join(BASE, "data_clinical_sample.csv"), sep=",", low_memory=False,
                  usecols=["sampleId", "CANCER_TYPE", "CANCER_TYPE_DETAILED", "GENE_PANEL",
                           "SAMPLE_TYPE", "TMB_NONSYNONYMOUS", "FRACTION_GENOME_ALTERED",
                           "TUMOR_PURITY", "MSI_TYPE"])
mut = pd.read_csv(os.path.join(BASE, "data_mutations.csv"), sep=",",
                  usecols=["sampleId", "patientId", "entrezGeneId", "mutationStatus",
                           "variantType", "mutationType", "proteinChange"], low_memory=False)
print("mutation rows:", len(mut), "| unique samples:", mut.sampleId.nunique())

print("\nmutationStatus:", mut.mutationStatus.value_counts(dropna=False).to_dict())
print("variantType:", mut.variantType.value_counts(dropna=False).to_dict())
print("mutationType top:", mut.mutationType.value_counts(dropna=False).head(10).to_dict())

# attach panel
mut = mut.merge(smp[["sampleId", "GENE_PANEL", "CANCER_TYPE", "CANCER_TYPE_DETAILED", "SAMPLE_TYPE"]],
                on="sampleId", how="left")
print("\nrows without panel:", mut.GENE_PANEL.isna().sum())

# --- observed gene universe per panel ---
panels = ["IMPACT341", "IMPACT410", "IMPACT468", "IMPACT505"]
sets = {}
for p in panels:
    s = set(mut.loc[mut.GENE_PANEL == p, "entrezGeneId"].dropna().astype(int))
    sets[p] = s
    print(f"{p}: {len(s)} genes observed | samples with mutations: {mut.loc[mut.GENE_PANEL==p,'sampleId'].nunique()}")

core = set.intersection(*sets.values())
print(f"\nINTERSECTION (all four panels): {len(core)} genes")
print(f"union: {len(set.union(*sets.values()))} genes")

# sanity: canonical genes present?
gmap = pd.read_csv(_GENEMAP)
gmap.columns = [c.strip().lstrip("\ufeff") for c in gmap.columns]
gmap["entrezGeneId"] = pd.to_numeric(gmap["entrezGeneId"], errors="coerce")
e2h = dict(zip(gmap.entrezGeneId.dropna().astype(int), gmap.hugoGeneSymbol))
h2e = {v: k for k, v in e2h.items()}
check = ["TP53", "KRAS", "EGFR", "KEAP1", "STK11", "SMARCA4", "CDKN2A", "NFE2L2", "RB1",
         "NKX2-1", "MYC", "MDM2", "PIK3CA", "FOXA1", "RBM10", "ARID1A", "ERBB2", "BRAF",
         "MET", "ALK", "RET", "ROS1", "ATM", "STK11", "TERT", "APC", "NOTCH1", "FAT1", "PTPRD"]
print("\ncanonical gene check (in core set?):")
for g in check:
    e = h2e.get(g)
    print(f"  {g:10s} entrez={e} in_core={e in core if e else 'NA'} "
          f"per-panel={{ {', '.join(p+':'+str(e in sets[p]) for p in panels)} }}")

core_df = pd.DataFrame({
    "entrezGeneId": sorted(core),
    "hugoGeneSymbol": [e2h.get(e, "") for e in sorted(core)],
})
core_df.to_csv(os.path.join(OUT, "panel_core_genes.csv"), index=False)
print("\nwrote panel_core_genes.csv:", len(core_df))

# --- novel-cancer-type carriers within NSCLC ---
nsclc_det = ["Lung Adenocarcinoma", "Lung Squamous Cell Carcinoma",
             "Non-Small Cell Lung Cancer", "Poorly Differentiated Non-Small Cell Lung Cancer"]
sub = mut[mut.CANCER_TYPE_DETAILED.isin(nsclc_det) & mut.entrezGeneId.isin(core)]
print("\nNSCLC mutation rows in core set:", len(sub),
      "| patients:", sub.patientId.nunique())
# patients with >=1 core mutation, by panel
tmp = sub.drop_duplicates(["patientId", "entrezGeneId"]).groupby("GENE_PANEL").patientId.nunique()
print("patients with >=1 core mutation, by panel:\n", tmp.to_string())

print("\ndone.")
