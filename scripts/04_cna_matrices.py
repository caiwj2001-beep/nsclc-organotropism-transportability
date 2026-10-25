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
"""STAGE 2B: CNA gene universe, core-set intersection, patient-level CNA matrices,
CNA landscape per histology. Saves matrices for later stages."""
import os
import numpy as np
import pandas as pd

BASE = _CHORD
OUT = _os.path.abspath(_OUT)

pat = pd.read_csv(os.path.join(OUT, "stage2_patient_table.csv"), index_col=0, low_memory=False)
core = pd.read_csv(os.path.join(OUT, "panel_core_genes.csv"))
core_e = set(core.entrezGeneId.astype(int))
gmap = pd.read_csv(_GENEMAP)
gmap.columns = [c.strip().lstrip("\ufeff") for c in gmap.columns]
gmap["entrezGeneId"] = pd.to_numeric(gmap.entrezGeneId, errors="coerce")
e2h = dict(zip(gmap.entrezGeneId.dropna().astype(int), gmap.hugoGeneSymbol))

# ---- 1. universe ----
print("=" * 100)
print("1. CNA GENE UNIVERSE")
print("=" * 100)
uni = set()
reader = pd.read_csv(os.path.join(BASE, "data_cna_long.csv"), sep=",",
                     usecols=["entrezGeneId"], chunksize=3_000_000, low_memory=False)
for ch in reader:
    uni.update(ch.entrezGeneId.dropna().astype(int).unique())
print("CNA genes in file:", len(uni))
inter = sorted(uni & core_e)
print("CNA ∩ 335-gene core set:", len(inter))
missing = sorted(core_e - uni)
print("core genes WITHOUT CNA data:", len(missing),
      "| examples:", [e2h.get(e) for e in missing[:15]])
pd.DataFrame({"entrezGeneId": sorted(uni), "hugoGeneSymbol": [e2h.get(e, "") for e in sorted(uni)],
              "in_core335": [e in core_e for e in sorted(uni)]}
             ).to_csv(os.path.join(OUT, "cna_gene_universe.csv"), index=False)

# ---- 2. matrices for NSCLC patients, chosen sample ----
ns = pat[pat.grp.isin(["LUAD", "LUSC", "NSCLC-NOS", "NSCLC-PD"])].copy()
samp2pat = dict(zip(ns.sampleId, ns.index))
print("\nNSCLC patients with a chosen sample:", len(samp2pat))
rows = []
reader = pd.read_csv(os.path.join(BASE, "data_cna_long.csv"), sep=",",
                     usecols=["sampleId", "entrezGeneId", "value"], chunksize=3_000_000, low_memory=False)
for ch in reader:
    ch = ch[ch.sampleId.isin(samp2pat) & ch.entrezGeneId.isin(set(inter))]
    if len(ch):
        rows.append(ch)
cna = pd.concat(rows, ignore_index=True)
print("filtered CNA rows (NSCLC, core∩CNA):", len(cna))
cna["value"] = cna.value.astype(float)
cna["patientId"] = cna.sampleId.map(samp2pat)
cna["gene"] = cna.entrezGeneId.map(e2h)

deep = (cna.assign(v=(cna.value <= -1.5).astype(int))
        .pivot_table(index="patientId", columns="gene", values="v", aggfunc="max"))
amp = (cna.assign(v=(cna.value == 2).astype(int))
       .pivot_table(index="patientId", columns="gene", values="v", aggfunc="max"))
raw = cna.pivot_table(index="patientId", columns="gene", values="value", aggfunc="first")
print("matrix shapes: deep", deep.shape, "amp", amp.shape, "raw", raw.shape)

deep = deep.reindex(ns.index).fillna(0).astype(np.int8)
amp = amp.reindex(ns.index).fillna(0).astype(np.int8)
raw = raw.reindex(ns.index).astype(np.float32)
deep.index.name = amp.index.name = raw.index.name = "patientId"
deep.to_csv(os.path.join(OUT, "stage2_cna_deepdel.csv.gz"), compression="gzip")
amp.to_csv(os.path.join(OUT, "stage2_cna_amp.csv.gz"), compression="gzip")
raw.to_csv(os.path.join(OUT, "stage2_cna_raw.csv.gz"), compression="gzip")
print("saved cna matrices (.csv.gz)")

# ---- 3. landscape ----
print("\n" + "=" * 100)
print("3. CNA LANDSCAPE (% of patients, chosen sample)")
print("=" * 100)
land = []
for g in ["LUAD", "LUSC", "NSCLC-NOS"]:
    ids = ns.index[ns.grp == g]
    if len(ids) < 20:
        continue
    land.append(pd.DataFrame({"amp": amp.loc[ids].mean() * 100,
                              "del": deep.loc[ids].mean() * 100}).assign(grp=g))
L = pd.concat(land)
for g in ["LUAD", "LUSC"]:
    d = L[L.grp == g]
    print(f"\n--- {g}: TOP 12 amplifications")
    print(d.sort_values("amp", ascending=False).head(12)["amp"].round(1).to_string())
    print(f"--- {g}: TOP 12 deep deletions")
    print(d.sort_values("del", ascending=False).head(12)["del"].round(1).to_string())
L.to_csv(os.path.join(OUT, "stage2_cna_landscape_long.csv"))

wide = L.reset_index().pivot_table(index="gene", columns="grp", values=["amp", "del"])
wide.round(1).to_csv(os.path.join(OUT, "stage2_cna_landscape_wide.csv"))
print("\nwrote stage2_cna_landscape_wide.csv", wide.shape)
