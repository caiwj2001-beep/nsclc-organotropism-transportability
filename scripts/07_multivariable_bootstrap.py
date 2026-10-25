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
"""STAGE 4: multivariable models + bootstrap(1000) + amplicon conditionalization +
minimum detectable interaction ratio (MDE).

Main model (per organ, per stratum):
  logit P(organ) = b0 + b1*exposure + age + male + smoking dummies + stage4
                   + TMB + FGA + purity + panel505 + log(n_distant)
Pooled interaction:
  ... + hist + exposure*hist   -> ratio_OR = exp(b_int); MDE = exp((1.96+0.84)*SE(b_int))
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = _CHORD
OUT = _os.path.abspath(_OUT)
os.makedirs(OUT, exist_ok=True)
D2 = _os.path.abspath(_OUT)

pat = pd.read_csv(os.path.join(D2, "stage2_patient_table.csv"), index_col=0, low_memory=False)
pat = pat[pat.grp.isin(["LUAD", "LUSC"])]
core = pd.read_csv(os.path.join(D2, "panel_core_genes.csv"))
core_e = set(core.entrezGeneId.astype(int))
h2e = dict(zip(core.entrezGeneId.astype(int), core.hugoGeneSymbol))

# exposure matrices
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
df["lgorg"] = np.log(df.n_distant)
df["panel505"] = (df.GENE_PANEL == "IMPACT505").astype(int)
df["stage4"] = df.stage4.astype(int)
sm = df.SMOKING_PREDICTIONS_3_CLASSES.astype(str)
# only two informative levels exist in MSK-CHORD: Never / Former-Current Smoker; Unknown -> missing
df["smk_fc"] = np.where(sm.str.strip() == "Former/Current Smoker", 1,
                        np.where(sm.str.strip() == "Never", 0, np.nan))
df["logtmb"] = np.log1p(pd.to_numeric(df.tmb, errors="coerce"))
df["hist"] = (df.grp == "LUSC").astype(int)
print("analysis frame (distant+):", df.shape, "| LUAD", int((df.grp == "LUAD").sum()),
      "| LUSC", int((df.grp == "LUSC").sum()))

COVS = ["age", "male", "smk_fc", "stage4", "logtmb", "fga", "purity", "panel505", "lgorg"]


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
            b = nb
            break
        b = nb
    eta = A @ b
    p = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
    W = np.clip(p * (1 - p), 1e-8, None)
    H = (A * W[:, None]).T @ A
    try:
        cov = np.linalg.inv(H)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(H)
    se = np.sqrt(np.abs(np.diag(cov)))
    if not np.all(np.isfinite(b)) or not np.all(np.isfinite(se)):
        return None
    return b, se


def design(d, xcols):
    d = d.dropna(subset=xcols + COVS)
    A = np.column_stack([np.ones(len(d))] + [d[c].astype(float).values for c in xcols + COVS])
    return d, A


TARGETS = [
    ("M_STK11", "ADRENAL_GLANDS"), ("M_STK11", "PLEURA"), ("M_STK11", "LIVER"),
    ("M_SMARCA4", "ADRENAL_GLANDS"), ("M_SMARCA4", "BONE"),
    ("M_KEAP1", "ADRENAL_GLANDS"), ("M_EGFR", "ADRENAL_GLANDS"), ("M_EGFR", "CNS_BRAIN"),
    ("M_KRAS", "ADRENAL_GLANDS"),
    ("A_EGFR", "CNS_BRAIN"), ("A_CDK4", "CNS_BRAIN"), ("A_SDHA", "CNS_BRAIN"),
    ("A_IL7R", "CNS_BRAIN"), ("D_CDKN2A", "CNS_BRAIN"), ("A_FGFR1", "LIVER"),
]

# ---------- 1. multivariable per stratum ----------
rows = []
for x, org in TARGETS:
    if x not in df.columns:
        continue
    for stratum in ["LUAD", "LUSC"]:
        d = df[(df.grp == stratum) & df[org + "_b"].notna()].copy()
        d, A = design(d, [x])
        if len(d) < 50 or A is None:
            continue
        y = d[org + "_b"].astype(float).values
        r = fit(y, A)
        if r is None:
            continue
        b, se = r
        if abs(b[1]) > 10 or se[1] <= 0:
            continue          # separation / degenerate
        z = b[1] / se[1]
        rows.append({"exposure": x, "organ": org, "stratum": stratum, "n": len(d),
                     "cases": int(y.sum()), "carriers": int(d[x].sum()),
                     "OR": np.exp(b[1]), "lo": np.exp(b[1] - 1.96 * se[1]),
                     "hi": np.exp(b[1] + 1.96 * se[1]),
                     "p": 2 * (1 - norm.cdf(abs(z)))})
mv = pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
mv["q"] = np.nan
for (o, s), g in mv.groupby(["organ", "stratum"]):
    mv.loc[g.index, "q"] = multipletests(g.p, method="fdr_bh")[1]
mv = mv.sort_values("p")
mv.to_csv(os.path.join(OUT, "stage4_multivariable.csv"), index=False)
print("\n=== 多因素（全协变量 + 负荷）")
print(mv[mv.stratum == "LUAD"].head(20)[["exposure", "organ", "n", "carriers", "OR", "lo", "hi", "p", "q"]].to_string(index=False))
print("\n--- LUSC 同靶点")
print(mv[mv.stratum == "LUSC"][["exposure", "organ", "n", "carriers", "OR", "lo", "hi", "p"]].to_string(index=False))

# ---------- 2. bootstrap (LUAD) ----------
rng = np.random.default_rng(20261025)
boot_rows = []
for x, org in [(t[0], t[1]) for t in TARGETS]:
    if x not in df.columns:
        continue
    d = df[(df.grp == "LUAD") & df[org + "_b"].notna()].copy()
    d, A = design(d, [x])
    if len(d) < 100:
        continue
    y = d[org + "_b"].astype(float).values
    ors = []
    for _ in range(1000):
        idx = rng.integers(0, len(d), len(d))
        Ab = A[idx]; yb = y[idx]
        if yb.min() == yb.max():
            continue
        r = fit(yb, Ab, max_iter=60)
        if r is not None and abs(r[0][1]) < 10:
            ors.append(float(np.exp(r[0][1])))
    if ors:
        boot_rows.append({"exposure": x, "organ": org, "n_boot": len(ors),
                          "boot_med": float(np.median(ors)),
                          "boot_2.5": float(np.percentile(ors, 2.5)),
                          "boot_97.5": float(np.percentile(ors, 97.5))})
bt = pd.DataFrame(boot_rows)
bt.to_csv(os.path.join(OUT, "stage4_bootstrap.csv"), index=False)
print("\n=== Bootstrap 1,000（LUAD，百分位 CI）")
print(bt.round(3).to_string(index=False))

# ---------- 3. amplicon conditionalization ----------
AMP_OF = {"A_CDK4": "A_MDM2", "A_SDHA": "A_MDM2", "A_IL7R": "A_MDM2",
          "D_CDKN2A": "D_CD274", "D_CD274": "D_CDKN2A", "A_FGFR1": None}
amp_rows = []
for x, org in [("A_CDK4", "CNS_BRAIN"), ("A_SDHA", "CNS_BRAIN"), ("A_IL7R", "CNS_BRAIN"),
               ("D_CDKN2A", "CNS_BRAIN"), ("A_FGFR1", "LIVER")]:
    partner = AMP_OF.get(x)
    d = df[(df.grp == "LUAD") & df[org + "_b"].notna()].copy()
    xc = [x] if not partner or partner not in d.columns else [x, partner]
    d, A = design(d, xc)
    if len(d) < 100:
        continue
    y = d[org + "_b"].astype(float).values
    r = fit(y, A)
    if r is None:
        continue
    b, se = r
    amp_rows.append({"exposure": x, "organ": org, "partner": partner or "none",
                     "n": len(d), "carriers": int(d[x].sum()),
                     "OR_adj_burden_only": np.nan,
                     "OR_conditional": np.exp(b[1]), "lo": np.exp(b[1] - 1.96 * se[1]),
                     "hi": np.exp(b[1] + 1.96 * se[1]),
                     "p": 2 * (1 - norm.cdf(abs(b[1] / se[1]))),
                     "partner_OR": np.exp(b[2]) if len(xc) == 2 else np.nan})
ac = pd.DataFrame(amp_rows)
ac.to_csv(os.path.join(OUT, "stage4_amplicon_conditional.csv"), index=False)
print("\n=== 扩增子条件化（LUAD）")
print(ac.round(3).to_string(index=False))

# ---------- 4. interaction + MDE ----------
int_rows = []
for x, org in TARGETS:
    if x not in df.columns:
        continue
    d = df[df[org + "_b"].notna()].copy()
    d = d.dropna(subset=[x] + COVS)
    A = np.column_stack([np.ones(len(d)), d[x].astype(float), d["hist"].astype(float),
                         (d[x].astype(float) * d["hist"].astype(float)),
                         d["lgorg"].astype(float)] +
                        [d[c].astype(float) for c in COVS if c != "lgorg"])
    y = d[org + "_b"].astype(float).values
    r = fit(y, A)
    if r is None:
        continue
    b, se = r
    se_int = se[3]
    int_rows.append({"exposure": x, "organ": org, "n": len(d),
                     "carriers_LUAD": int(d.loc[d["hist"] == 0, x].sum()),
                     "carriers_LUSC": int(d.loc[d["hist"] == 1, x].sum()),
                     "OR_LUAD": np.exp(b[1]), "OR_LUSC": np.exp(b[1] + b[3]),
                     "ratio_OR": np.exp(b[3]),
                     "p_interaction": 2 * (1 - norm.cdf(abs(b[3] / se_int))),
                     "MDE_ratio_alpha05": np.exp((1.96 + 0.84) * se_int),
                     "MDE_ratio_alpha001": np.exp((3.29 + 0.84) * se_int)})
it = pd.DataFrame(int_rows)
if len(it):
    it["q"] = multipletests(it.p_interaction, method="fdr_bh")[1]
    it = it.sort_values("p_interaction")
it.to_csv(os.path.join(OUT, "stage4_interaction_mde.csv"), index=False)
print("\n=== 交互 + 最小可检出 ratio_OR")
print(it.round(3).to_string(index=False))
print("\nDONE")
