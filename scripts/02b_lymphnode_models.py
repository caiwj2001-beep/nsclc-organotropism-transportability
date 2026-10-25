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
"""STEP B addendum: LYMPH_NODES site-specific models (missed in step_b because the
distant-organ list excludes LYMPH_NODES). Lengel's second MV target is TP53 -> nodal
metastasis OR 2.03 (1.28-3.23). Analyses reuse the saved LUAD table."""
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact

BASE = _CHORD
OUT = _os.path.abspath(_OUT)
df = pd.read_csv(os.path.join(OUT, "luad_analysis_table.csv"), low_memory=False, index_col=0)
df.index.name = "patientId"
df = df.reset_index()
pat = pd.read_csv(os.path.join(BASE, "data_clinical_patient.csv"), sep=",", low_memory=False,
                  usecols=["patientId", "LYMPH_NODES", "LUNG"])
pat = pat.set_index("patientId")
for c in ["LYMPH_NODES", "LUNG"]:
    df[c + "_b"] = (df.patientId.map(pat[c]).astype(str).str.strip().str.lower() == "yes")
print("ever_met:", int(df.ever_met.sum()), "| LN+ among ever_met:",
      int(df.loc[df.ever_met, "LYMPH_NODES_b"].sum()))

em = df[df.ever_met].copy()
COVS = ["age", "male", "tmb", "fga", "purity", "stage4"]


def fit(y, x, covs, d):
    d = d[[y, x] + covs].dropna()
    if d[y].nunique() < 2 or d[x].nunique() < 2:
        return None
    X = sm.add_constant(d[[x] + covs].astype(float))
    try:
        r = sm.Logit(d[y].astype(float), X).fit(disp=0, maxiter=200)
    except Exception:
        return None
    b, se = r.params[x], r.bse[x]
    return dict(OR=float(np.exp(b)), lo=float(np.exp(b - 1.96 * se)),
                hi=float(np.exp(b + 1.96 * se)), p=float(r.pvalues[x]),
                n=int(len(d)), cases=int(d[y].sum()))


rows = []
for g in ["TP53", "SMARCA4", "CDKN2A", "KEAP1", "STK11", "NFE2L2", "KRAS", "EGFR"]:
    col = "alt_" + g
    if col not in em.columns:
        continue
    d = em[[c for c in ["LYMPH_NODES_b", col] + COVS]].dropna()
    a = int(((d.LYMPH_NODES_b) & (d[col] == 1)).sum()); b1 = int(((d.LYMPH_NODES_b) & (d[col] == 0)).sum())
    c = int(((~d.LYMPH_NODES_b) & (d[col] == 1)).sum()); d1 = int(((~d.LYMPH_NODES_b) & (d[col] == 0)).sum())
    orr, pc = fisher_exact([[a, b1], [c, d1]])
    mv = fit("LYMPH_NODES_b", col, COVS, em)
    if mv:
        rows.append(dict(organ="LYMPH_NODES", gene=g,
                         rate_in=round(float(em.loc[em.LYMPH_NODES_b, col].mean()), 3),
                         rate_out=round(float(em.loc[~em.LYMPH_NODES_b, col].mean()), 3),
                         OR_crude=round(orr, 2), p_crude=pc, OR_adj=round(mv["OR"], 2),
                         lo=round(mv["lo"], 2), hi=round(mv["hi"], 2), p_adj=mv["p"],
                         n=mv["n"], cases=mv["cases"], tbl=[a, b1, c, d1]))
t = pd.DataFrame(rows)
t["q_adj"] = multipletests(t.p_adj, method="fdr_bh")[1]
print(t.sort_values("p_adj").to_string(index=False))
t.to_csv(os.path.join(OUT, "replication_step2b_lymphnode.csv"), index=False)

# merge into the main site-specific table for a complete record
main = pd.read_csv(os.path.join(OUT, "replication_step2_site_specific.csv"))
t2 = t.rename(columns={"rate_in": "carrier_rate_in", "rate_out": "carrier_rate_out"}).drop(columns=["tbl"])
full = pd.concat([main, t2], ignore_index=True)
full.to_csv(os.path.join(OUT, "replication_step2_site_specific_FULL.csv"), index=False)
print("\nfull site table rows:", len(full), "| organs:", sorted(full.organ.unique()))
