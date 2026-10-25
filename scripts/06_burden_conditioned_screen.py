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
"""STAGE 3: burden-conditioned single-marker screening.
Model per organ o, exposure x (mutation carrier / CNA amp / CNA deep deletion), stratum s:
    logit P(organ_o = 1) = b0 + b1*x + b2*log(n_distant)
Population: patients with >=1 distant metastasis (burden-conditioned), per histology.
Interaction: logit = b0 + b1*x + b2*hist + b3*x*hist + b4*log(n_distant) on LUAD+LUSC pooled.
BH-FDR within each (exposure type x organ x stratum) family; global BH reported separately.
"""
import os, sys, json
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.stats import norm

OUT = _os.path.abspath(_OUT)
BASE = _CHORD
ORGS = ["BONE", "LIVER", "CNS_BRAIN", "PLEURA", "ADRENAL_GLANDS", "INTRA_ABDOMINAL"]
MIN_CARRIERS = 10
MIN_CELL = 5

pat = pd.read_csv(os.path.join(OUT, "stage2_patient_table.csv"), index_col=0, low_memory=False)
pat = pat[pat.grp.isin(["LUAD", "LUSC", "NSCLC-NOS"])]
pat["hist"] = (pat.grp == "LUSC").astype(int)          # LUAD=0, LUSC=1
pat["logorg"] = np.log(pat.n_distant.clip(lower=1))

# --- mutation carrier matrix (mutation only, chosen sample) ---
core = pd.read_csv(os.path.join(OUT, "panel_core_genes.csv"))
core_e = set(core.entrezGeneId.astype(int))
mut = pd.read_csv(os.path.join(BASE, "data_mutations.csv"), sep=",",
                  usecols=["sampleId", "entrezGeneId"], low_memory=False)
samp2pat = dict(zip(pat.sampleId, pat.index))
mut = mut[mut.sampleId.isin(samp2pat) & mut.entrezGeneId.isin(core_e)].drop_duplicates()
mut["patientId"] = mut.sampleId.map(samp2pat)
MUT = (mut.assign(v=1).pivot_table(index="patientId", columns="entrezGeneId", values="v",
                                  aggfunc="max").reindex(pat.index).fillna(0).astype(np.int8))
MUT.columns = [f"M_{c}" for c in MUT.columns]
print("mutation matrix:", MUT.shape)

AMP = pd.read_csv(os.path.join(OUT, "stage2_cna_amp.csv.gz"), index_col=0)
DEL = pd.read_csv(os.path.join(OUT, "stage2_cna_deepdel.csv.gz"), index_col=0)
h2e = dict(zip(core.entrezGeneId.astype(int), core.hugoGeneSymbol))
AMP.columns = [f"A_{c}" for c in AMP.columns]
DEL.columns = [f"D_{c}" for c in DEL.columns]
print("cna matrices:", AMP.shape, DEL.shape)

X = pd.concat([pat[["grp", "hist", "logorg", "n_distant"] + [o + "_b" for o in ORGS]], MUT, AMP, DEL],
              axis=1, join="inner").copy()
print("analysis matrix after inner join:", X.shape)
numcols = list(MUT.columns) + list(AMP.columns) + list(DEL.columns)
X[numcols] = X[numcols].fillna(0)
for c in MUT.columns:
    X[c] = X[c].astype(int)
for c in list(AMP.columns) + list(DEL.columns):
    X[c] = X[c].astype(int)



def gname(c):
    k = c.split("_", 1)[1]
    try:
        return h2e.get(int(k), k)
    except ValueError:
        return k   # CNA matrices are keyed by HUGO symbol

def fit_logit(y, Xm, cols, max_iter=60, tol=1e-8):
    """Newton-IRLS; returns params, se, or None on failure/separation."""
    d = pd.DataFrame({"y": y})
    for c in cols:
        d[c] = Xm[c].values
    d = d.dropna()
    yv = d.y.values.astype(float)
    if yv.min() == yv.max():
        return None
    A = np.column_stack([np.ones(len(d))] + [d[c].values.astype(float) for c in cols])
    b = np.zeros(A.shape[1])
    for _ in range(max_iter):
        eta = A @ b
        p = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
        W = np.clip(p * (1 - p), 1e-8, None)
        z = eta + (yv - p) / W
        try:
            nb = np.linalg.solve((A * W[:, None]).T @ A, (A * W[:, None]).T @ z)
        except np.linalg.LinAlgError:
            return None
        if np.max(np.abs(nb - b)) < tol:
            b = nb
            break
        b = nb
    eta = A @ b
    p = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
    W = np.clip(p * (1 - p), 1e-8, None)
    try:
        cov = np.linalg.inv((A * W[:, None]).T @ A)
    except np.linalg.LinAlgError:
        return None
    se = np.sqrt(np.diag(cov))
    return b, se, len(d)


def screen(expo_prefix, label):
    rows = []
    for org in ORGS:
        y = X[org + "_b"].astype(int)
        for stratum, mask in [("LUAD", X.grp == "LUAD"), ("LUSC", X.grp == "LUSC")]:
            # restrict to distant+ (burden condition), same index space as X
            sub = X.loc[mask & (X.n_distant > 0)]
            yv = sub[org + "_b"].astype(int)
            for c in [c for c in X.columns if c.startswith(expo_prefix)]:
                x = sub[c].astype(int)
                nc = int(x.sum())
                if nc < MIN_CARRIERS:
                    continue
                a = int(((yv == 1) & (x == 1)).sum()); b1 = int(((yv == 1) & (x == 0)).sum())
                cc = int(((yv == 0) & (x == 1)).sum()); d1 = int(((yv == 0) & (x == 0)).sum())
                if min(a, b1, cc, d1) < MIN_CELL:
                    continue
                r = fit_logit(yv.values, sub, [c, "logorg"])
                if r is None:
                    continue
                bb, se, n = r
                rows.append({"exposure": label, "organ": org, "stratum": stratum,
                             "gene_id": c.split("_", 1)[1],
                             "gene": gname(c),
                             "n": n, "cases": int(yv.sum()), "carriers": nc,
                             "rate_in": round(100 * a / max(a + cc, 1), 1),
                             "rate_out": round(100 * b1 / max(b1 + d1, 1), 1),
                             "OR": float(np.exp(bb[1])),
                             "lo": float(np.exp(bb[1] - 1.96 * se[1])),
                             "hi": float(np.exp(bb[1] + 1.96 * se[1])),
                             "p": float(2 * (1 - norm.cdf(abs(bb[1] / se[1]))))})
    t = pd.DataFrame(rows)
    t["q_family"] = np.nan
    for (o, s), g in t.groupby(["organ", "stratum"]):
        t.loc[g.index, "q_family"] = multipletests(g.p, method="fdr_bh")[1]
    t["q_global"] = multipletests(t.p, method="fdr_bh")[1]
    return t.sort_values(["organ", "stratum", "p"])


def interaction(expo_prefix, label):
    rows = []
    pool = X[(X.n_distant > 0) & X.grp.isin(["LUAD", "LUSC"])].copy()
    pool["hist"] = (pool.grp == "LUSC").astype(int)
    for org in ORGS:
        yv = pool[org + "_b"].astype(int)
        for c in [c for c in pool.columns if c.startswith(expo_prefix)]:
            x = pool[c].astype(int)
            if x.sum() < 15:
                continue
            # need both histologies to have carriers
            hv = pool["hist"].values; xv = x.values
            n_lusc = int(xv[hv == 1].sum()); n_luad = int(xv[hv == 0].sum())
            if min(n_lusc, n_luad) < 3:
                continue
            d = pd.DataFrame({"y": yv.values, "x": x.values, "h": pool["hist"].values,
                              "lg": pool.logorg.values})
            if d.y.min() == d.y.max():
                continue
            A = np.column_stack([np.ones(len(d)), d.x, d.h, d.x * d.h, d.lg])
            r = fit_glm(yv.values, A)
            if r is None:
                continue
            bb, se = r
            z = bb[3] / se[3]
            rows.append({"exposure": label, "organ": org,
                         "gene_id": c.split("_", 1)[1],
                         "gene": gname(c),
                         "OR_LUAD": float(np.exp(bb[1])),
                         "OR_LUSC": float(np.exp(bb[1] + bb[3])),
                         "ratio_OR": float(np.exp(bb[3])),
                         "p_interaction": float(2 * (1 - norm.cdf(abs(z)))),
                         "carriers_LUAD": n_luad, "carriers_LUSC": n_lusc,
                         "n": len(d)})
    t = pd.DataFrame(rows)
    if len(t):
        t["q_family"] = np.nan
        for o, g in t.groupby("organ"):
            t.loc[g.index, "q_family"] = multipletests(g.p_interaction, method="fdr_bh")[1]
        t["q_global"] = multipletests(t.p_interaction, method="fdr_bh")[1]
        t = t.sort_values("p_interaction")
    return t


def fit_glm(yv, A, max_iter=80, tol=1e-9):
    b = np.zeros(A.shape[1])
    for _ in range(max_iter):
        eta = A @ b
        p = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
        W = np.clip(p * (1 - p), 1e-8, None)
        z = eta + (yv - p) / W
        try:
            nb = np.linalg.solve((A * W[:, None]).T @ A, (A * W[:, None]).T @ z)
        except np.linalg.LinAlgError:
            return None
        if np.max(np.abs(nb - b)) < tol:
            b = nb; break
        b = nb
    eta = A @ b
    p = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
    W = np.clip(p * (1 - p), 1e-8, None)
    try:
        cov = np.linalg.inv((A * W[:, None]).T @ A)
    except np.linalg.LinAlgError:
        return None
    return b, np.sqrt(np.diag(cov))


for prefix, label in [("M_", "mutation"), ("A_", "cna_amp"), ("D_", "cna_deepdel")]:
    t = screen(prefix, label)
    fp = os.path.join(OUT, f"stage3_screen_{label}.csv")
    t.to_csv(fp, index=False)
    print(f"\n=== {label}: {len(t)} tests | q_family<0.05: {int((t.q_family<0.05).sum())} "
          f"| q_global<0.05: {int((t.q_global<0.05).sum())} | wrote {os.path.basename(fp)}")
    tt = t[t.q_family < 0.05].sort_values("p")
    if len(tt):
        print(tt.head(15)[["organ", "stratum", "gene", "carriers", "rate_in", "rate_out",
                           "OR", "lo", "hi", "p", "q_family"]].to_string(index=False))
    it = interaction(prefix, label)
    ifp = os.path.join(OUT, f"stage3_interaction_{label}.csv")
    it.to_csv(ifp, index=False)
    print(f"--- interaction {label}: {len(it)} tests | q_family<0.05: {int((it.q_family<0.05).sum()) if len(it) else 0}")
    if len(it):
        print(it.head(10)[["organ", "gene", "OR_LUAD", "OR_LUSC", "ratio_OR",
                           "p_interaction", "q_family"]].to_string(index=False))
print("\nDONE")
