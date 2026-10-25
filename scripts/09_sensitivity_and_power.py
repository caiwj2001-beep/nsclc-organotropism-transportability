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
"""STAGE 6: internal validation & sensitivity
1) single-organ subset (burden-free tropism comparison)
2) primary-sample-only subset
3) missing-indicator / imputed covariates (recover the 21% lost to complete-case)
4) interaction power curve: MDE(ratio_OR) as a function of LUSC carrier count
Outputs a consolidated numbers registry for the manuscript.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

BASE = _CHORD
P2 = _os.path.abspath(_OUT)
OUT = _os.path.abspath(_OUT)
os.makedirs(OUT, exist_ok=True)
h2e = None

pat = pd.read_csv(os.path.join(P2, "stage2_patient_table.csv"), index_col=0, low_memory=False)
pat = pat[pat.grp.isin(["LUAD", "LUSC"])].copy()
core = pd.read_csv(os.path.join(P2, "panel_core_genes.csv"))
core_e = set(core.entrezGeneId.astype(int))
h2e = dict(zip(core.entrezGeneId.astype(int), core.hugoGeneSymbol))

mut = pd.read_csv(os.path.join(BASE, "data_mutations.csv"), sep=",",
                  usecols=["sampleId", "entrezGeneId"], low_memory=False)
s2p = dict(zip(pat.sampleId, pat.index))
mut = mut[mut.sampleId.isin(s2p) & mut.entrezGeneId.isin(core_e)].drop_duplicates()
mut["patientId"] = mut.sampleId.map(s2p)
M = mut.assign(v=1).pivot_table(index="patientId", columns="entrezGeneId", values="v", aggfunc="max")
M.columns = ["M_" + h2e.get(c, str(c)) for c in M.columns]
AMP = pd.read_csv(os.path.join(P2, "stage2_cna_amp.csv.gz"), index_col=0); AMP.columns = ["A_" + c for c in AMP.columns]
DEL = pd.read_csv(os.path.join(P2, "stage2_cna_deepdel.csv.gz"), index_col=0); DEL.columns = ["D_" + c for c in DEL.columns]

df = pat.join([M, AMP, DEL], how="inner").copy()
numcols = [c for c in df.columns if c.startswith(("M_", "A_", "D_"))]
df[numcols] = df[numcols].fillna(0).astype(int)
smk = df.SMOKING_PREDICTIONS_3_CLASSES.astype(str).str.strip()
df["smk_cat"] = np.where(smk == "Former/Current Smoker", 1, np.where(smk == "Never", 0, 2))
df["logtmb"] = np.log1p(pd.to_numeric(df.tmb, errors="coerce"))
df["panel505"] = (df.GENE_PANEL == "IMPACT505").astype(int)
df["stage4"] = df.stage4.astype(int)

TARGETS = [("M_STK11", "ADRENAL_GLANDS"), ("M_STK11", "PLEURA"), ("M_STK11", "LIVER"),
           ("M_SMARCA4", "ADRENAL_GLANDS"), ("M_SMARCA4", "BONE"),
           ("M_EGFR", "ADRENAL_GLANDS"), ("M_EGFR", "CNS_BRAIN"),
           ("A_CDK4", "CNS_BRAIN"), ("A_FGFR1", "LIVER"), ("D_CDKN2A", "CNS_BRAIN")]


def fit(y, A, max_iter=100, tol=1e-9):
    b = np.zeros(A.shape[1])
    for _ in range(max_iter):
        eta = A @ b
        p = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
        W = np.clip(p * (1 - p), 1e-8, None)
        z = eta + (y - p) / W
        H = (A * W[:, None]).T @ A
        try:
            nb = np.linalg.solve(H, (A * W[:, None]).T @ z)
        except np.linalg.LinAlgError:
            nb = np.linalg.lstsq(H, (A * W[:, None]).T @ z, rcond=None)[0]
        if np.max(np.abs(nb - b)) < tol:
            b = nb; break
        b = nb
    eta = A @ b
    p = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
    W = np.clip(p * (1 - p), 1e-8, None)
    H = (A * W[:, None]).T @ A
    cov = np.linalg.pinv(H)
    se = np.sqrt(np.abs(np.diag(cov)))
    return b, se


def model(d, x, org, covs, add_burden=True):
    covs = list(covs) + (["lgorg"] if add_burden else [])
    d = d.copy()
    d["lgorg"] = np.log(d.n_distant.clip(lower=1))
    d = d.dropna(subset=[x, org + "_b"] + covs)
    if len(d) < 40 or d[x].sum() < 8:
        return None
    y = d[org + "_b"].astype(float).values
    A = np.column_stack([np.ones(len(d)), d[x].astype(float)] +
                        [d[c].astype(float) for c in covs])
    r = fit(y, A)
    if r is None or abs(r[0][1]) > 10:
        return None
    b, se = r
    return dict(n=len(d), cases=int(y.sum()), carriers=int(d[x].sum()),
                OR=float(np.exp(b[1])), lo=float(np.exp(b[1] - 1.96 * se[1])),
                hi=float(np.exp(b[1] + 1.96 * se[1])),
                p=float(2 * (1 - norm.cdf(abs(b[1] / se[1])))))


BASE_COVS = ["age", "male", "smk_cat", "stage4", "logtmb", "fga", "purity", "panel505"]
rows = []

# ---- scenario 0: main (complete case, distant+) ----
EM = df[df.any_distant].copy()
EM["smk_cat"] = EM.smk_cat.where(EM.smk_cat != 2)      # Unknown -> NaN (complete case)
for x, org in TARGETS:
    for stratum in ["LUAD", "LUSC"]:
        r = model(EM[EM.grp == stratum], x, org, BASE_COVS)
        if r:
            rows.append(dict(scenario="main_complete_case", exposure=x, organ=org,
                             stratum=stratum, **r))

# ---- scenario 1: single-organ patients (burden-free) ----
SO = df[(df.any_distant) & (df.n_distant == 1)].copy()
SO["smk_cat"] = SO.smk_cat.where(SO.smk_cat != 2)
for x, org in TARGETS:
    for stratum in ["LUAD", "LUSC"]:
        r = model(SO[SO.grp == stratum], x, org, BASE_COVS, add_burden=False)
        if r:
            rows.append(dict(scenario="single_organ", exposure=x, organ=org,
                             stratum=stratum, **r))

# ---- scenario 2: primary-sample only ----
PO = df[(df.any_distant) & (df.SAMPLE_TYPE == "Primary")].copy()
PO["smk_cat"] = PO.smk_cat.where(PO.smk_cat != 2)
PO["lgorg"] = np.log(PO.n_distant)
for x, org in TARGETS:
    for stratum in ["LUAD", "LUSC"]:
        r = model(PO[PO.grp == stratum], x, org, BASE_COVS)
        if r:
            rows.append(dict(scenario="primary_only", exposure=x, organ=org,
                             stratum=stratum, **r))

# ---- scenario 3: missing-indicator (smoking Unknown as own level, purity imputed) ----
MI = df[df.any_distant].copy()
MI["purity"] = MI.purity.fillna(MI.purity.median())
MI["fga"] = MI.fga.fillna(MI.fga.median())
MI["logtmb"] = MI.logtmb.fillna(MI.logtmb.median())
MI["age"] = MI.age.fillna(MI.age.median())
MI["smk_unk"] = (MI.smk_cat == 2).astype(int)
MI["smk_cat"] = MI.smk_cat.where(MI.smk_cat != 2, 0)
MI_COVS = BASE_COVS + ["smk_unk"]
for x, org in TARGETS:
    for stratum in ["LUAD", "LUSC"]:
        r = model(MI[MI.grp == stratum], x, org, MI_COVS)
        if r:
            rows.append(dict(scenario="missing_indicator", exposure=x, organ=org,
                             stratum=stratum, **r))

sens = pd.DataFrame(rows)
sens.to_csv(os.path.join(OUT, "stage6_sensitivity.csv"), index=False)
piv = sens.pivot_table(index=["exposure", "organ", "stratum"], columns="scenario",
                       values=["OR", "p", "n"])
print("=== 敏感性：主分析 vs 单器官 vs 仅原发 vs 缺失指示（LUAD）")
lua = sens[(sens.stratum == "LUAD")].pivot_table(index=["exposure", "organ"], columns="scenario",
                                                values="OR")
print(lua.round(2).to_string())
print("\n=== 同时给出 p 与 n")
print(sens[(sens.stratum == "LUAD") & (sens.scenario.isin(["main_complete_case", "missing_indicator"]))]
      [["exposure", "organ", "scenario", "n", "carriers", "OR", "lo", "hi", "p"]].round(3).to_string(index=False))

# ---- 4) interaction power curve: simulation-based (observed-design, null main effect) ----
rng = np.random.default_rng(42)
pool = EM.copy()
pool["lgorg"] = np.log(pool.n_distant.clip(lower=1))
ORG = "ADRENAL_GLANDS"
base_cols = ["lgorg", "age", "male", "smk_cat", "stage4", "logtmb", "fga", "purity", "panel505"]
d0 = pool.dropna(subset=base_cols + [ORG + "_b"]).copy()
d0["hist"] = (d0.grp == "LUSC").astype(int)
y_obs = d0[ORG + "_b"].astype(float).values
# calibrate nuisance model: organ ~ hist + covariates (no exposure)
A0 = np.column_stack([np.ones(len(d0)), d0["hist"].astype(float),
                      d0[base_cols].astype(float)])
b0, _ = fit(y_obs, A0)
lin0 = A0 @ b0
print("calibration: observed adrenal rate {:.1f}%, model mean {:.1f}%".format(
    100 * y_obs.mean(), 100 * (1 / (1 + np.exp(-lin0))).mean()))

idx_lusc = np.where(d0.grp.values == "LUSC")[0]
idx_luad = np.where(d0.grp.values == "LUAD")[0]
# carriers must exist in BOTH histologies, otherwise x and x*hist are collinear
n_luad_carriers = int(round(0.153 * len(idx_luad)))     # ~ observed STK11 rate in LUAD
curve = []
for m in [20, 50, 100, 150, 250, 400]:
    for ratio in [1.5, 2.0, 2.5, 3.0, 4.0]:
        powers = []
        for _ in range(150):
            x = np.zeros(len(d0))
            x[rng.choice(idx_luad, n_luad_carriers, replace=False)] = 1
            x[rng.choice(idx_lusc, min(m, len(idx_lusc)), replace=False)] = 1
            lin = lin0 + np.log(ratio) * x * d0["hist"].values
            p = 1 / (1 + np.exp(-lin))
            ysim = rng.binomial(1, p)
            A = np.column_stack([np.ones(len(d0)), x, d0["hist"].astype(float),
                                 x * d0["hist"].astype(float), d0[base_cols].astype(float)])
            r = fit(ysim, A, max_iter=60)
            if r is None:
                continue
            b, se = r
            if se[3] <= 0 or not np.isfinite(se[3]):
                continue
            powers.append(2 * (1 - norm.cdf(abs(b[3] / se[3]))) < 0.05)
        if powers:
            curve.append({"LUSC_carriers": m, "true_ratio_OR": ratio,
                          "power_alpha05": float(np.mean(powers)), "n_sim": len(powers)})
cv = pd.DataFrame(curve)
cv.to_csv(os.path.join(OUT, "stage6_interaction_power_curve.csv"), index=False)
print("=== interaction power curve (adrenal, 150 sims/cell, alpha=0.05) ===")
print(cv.pivot(index="LUSC_carriers", columns="true_ratio_OR", values="power_alpha05").round(2).to_string())
print("LUSC distant-positive patients =", len(idx_lusc))
print("DONE")
