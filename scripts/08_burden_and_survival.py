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
"""STAGE 5:
(A) burden phenotype: n_distant ~ exposure (Poisson GLM, robust SE) for 335 mutations + CNA
(B) CNA landscape vs burden: FGA / amp-burden / del-burden correlations
(C) overall survival: Cox PH (statsmodels PHReg) for organs, burden and key genomic exposures,
    cross-checked against a hand-coded log-rank for two-group comparisons.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr, chi2

BASE = _CHORD
P2 = _os.path.abspath(_OUT)
OUT = _os.path.abspath(_OUT)
os.makedirs(OUT, exist_ok=True)

import statsmodels.api as sm
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests

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
AMP = pd.read_csv(os.path.join(P2, "stage2_cna_amp.csv.gz"), index_col=0)
DEL = pd.read_csv(os.path.join(P2, "stage2_cna_deepdel.csv.gz"), index_col=0)
AMP.columns = ["A_" + c for c in AMP.columns]
DEL.columns = ["D_" + c for c in DEL.columns]

df = pat.join([M, AMP, DEL], how="inner").copy()
numcols = [c for c in df.columns if c.startswith(("M_", "A_", "D_"))]
df[numcols] = df[numcols].fillna(0).astype(int)
smk_raw = df.SMOKING_PREDICTIONS_3_CLASSES.astype(str)
df["smk_fc"] = np.where(smk_raw.str.strip() == "Former/Current Smoker", 1,
                        np.where(smk_raw.str.strip() == "Never", 0, np.nan))
df["logtmb"] = np.log1p(pd.to_numeric(df.tmb, errors="coerce"))
df["panel505"] = (df.GENE_PANEL == "IMPACT505").astype(int)
df["stage4"] = df.stage4.astype(int)
df["hist"] = (df.grp == "LUSC").astype(int)
df["os_m"] = pd.to_numeric(df.OS_MONTHS, errors="coerce")
df["os_event"] = (df.OS_STATUS.astype(str).str.strip().str.lower().isin(["1", "deceased", "1:deceased"])).astype(int)
print("frame:", df.shape, "| LUAD", int((df.grp == "LUAD").sum()), "| LUSC", int((df.grp == "LUSC").sum()))
print("OS 可用:", int(df.os_m.notna().sum()), "| 事件:", int(df.os_event.sum()),
      "| OS_MONTHS<=0:", int((df.os_m <= 0).sum()))

COVS = ["age", "male", "smk_fc", "stage4", "logtmb", "fga", "purity", "panel505"]
EM = df[df.any_distant & (df.n_distant >= 1)].copy()      # burden / organ models
print("burden 分析集(远处转移+):", EM.shape)

# ---------------- (A) burden: Poisson with robust SE ----------------
def burden_scan(cols, label):
    rows = []
    for c in cols:
        d = EM.dropna(subset=[c, "n_distant"] + COVS)
        if d[c].sum() < 30:
            continue
        X = sm.add_constant(d[[c] + COVS].astype(float))
        try:
            r = sm.GLM(d.n_distant.astype(float), X, family=sm.families.Poisson()).fit(
                cov_type="HC0", maxiter=100)
        except Exception:
            continue
        b, se = r.params[c], r.bse[c]
        rows.append({"exposure": label, "gene": c.split("_", 1)[1], "n": int(len(d)),
                     "carriers": int(d[c].sum()), "RR": float(np.exp(b)),
                     "lo": float(np.exp(b - 1.96 * se)), "hi": float(np.exp(b + 1.96 * se)),
                     "p": float(r.pvalues[c])})
    t = pd.DataFrame(rows)
    if len(t):
        t["q"] = multipletests(t.p, method="fdr_bh")[1]
        t = t.sort_values("p")
    return t


all_mut = [c for c in M.columns if c in df.columns]
all_amp = [c for c in AMP.columns if c in df.columns]
all_del = [c for c in DEL.columns if c in df.columns]
b1 = burden_scan(all_mut, "mutation")
b2 = burden_scan(all_amp, "cna_amp")
b3 = burden_scan(all_del, "cna_deepdel")
bur = pd.concat([b1, b2, b3], ignore_index=True)
bur.to_csv(os.path.join(OUT, "stage5_burden_all.csv"), index=False)
print("\n=== 器官负荷（Poisson RR，全协变量）TOP 15")
print(bur.head(15).round(3).to_string(index=False))

# per-histology for the top hits
top_hits = bur.head(12)[["exposure", "gene"]].values.tolist()
rows = []
for label, g in top_hits:
    col = f"{'M' if label=='mutation' else ('A' if label=='cna_amp' else 'D')}_{g}"
    if col not in EM.columns:
        continue
    for stratum in ["LUAD", "LUSC"]:
        d = EM[EM.grp == stratum].dropna(subset=[col, "n_distant"] + COVS)
        if d[col].sum() < 10:
            continue
        X = sm.add_constant(d[[col] + COVS].astype(float))
        try:
            r = sm.GLM(d.n_distant.astype(float), X, family=sm.families.Poisson()).fit(
                cov_type="HC0", maxiter=100)
        except Exception:
            continue
        rows.append({"gene": g, "exposure": label, "stratum": stratum, "n": len(d),
                     "carriers": int(d[col].sum()), "RR": np.exp(r.params[col]),
                     "lo": np.exp(r.params[col] - 1.96 * r.bse[col]),
                     "hi": np.exp(r.params[col] + 1.96 * r.bse[col]),
                     "p": r.pvalues[col]})
pd.DataFrame(rows).to_csv(os.path.join(OUT, "stage5_burden_by_histology.csv"), index=False)
print("\n--- 负荷关联的分层（TOP 命中）")
print(pd.DataFrame(rows).round(3).to_string(index=False))

# ---------------- (B) CNA / FGA vs burden ----------------
print("\n=== FGA 与负荷、负荷与 OS 的描述统计")
for stratum in ["LUAD", "LUSC"]:
    d = EM[EM.grp == stratum]
    rho, p = spearmanr(d.fga.dropna(), d.loc[d.fga.notna(), "n_distant"])
    a = (d[all_amp].sum(axis=1)); dd = (d[all_del].sum(axis=1))
    rho2, p2 = spearmanr(a, d.n_distant)
    rho3, p3 = spearmanr(dd, d.n_distant)
    print(f"  {stratum}: FGA~n_org rho={rho:.3f} p={p:.2e} | #amp~n_org rho={rho2:.3f} p={p2:.2e}"
          f" | #del~n_org rho={rho3:.3f} p={p3:.3e}")

# ---------------- (C) survival ----------------
OS = EM[EM.os_m.notna() & (EM.os_m > 0)].copy()
print("\nOS 分析集:", OS.shape, "| 事件", int(OS.os_event.sum()))
SCOV = ["age", "male", "smk_fc", "stage4", "logtmb", "fga", "purity", "panel505", "n_distant"]


def cox(d, cols, covs):
    d = d.dropna(subset=["os_m", "os_event"] + cols + covs)
    if len(d) < 50 or d.os_event.sum() < 10:
        return None
    X = sm.add_constant(d[covs + cols].astype(float), has_constant="add")
    try:
        r = PHReg(d.os_m.values, X.values, status=d.os_event.values, ties="efron").fit()
    except Exception:
        return None
    out = {}
    names = ["const"] + covs + cols
    for c in cols:
        i = names.index(c)
        b, se = r.params[i], r.bse[i]
        out[c] = dict(HR=float(np.exp(b)), lo=float(np.exp(b - 1.96 * se)),
                      hi=float(np.exp(b + 1.96 * se)),
                      p=float(2 * (1 - norm.cdf(abs(b / se)))), n=int(len(d)),
                      events=int(d.os_event.sum()))
    return out


def logrank(t1, e1, t0, e0):
    """hand-coded two-group log-rank (Mantel-Haenszel)."""
    t = np.concatenate([t1, t0]); e = np.concatenate([e1, e0])
    g = np.concatenate([np.ones(len(t1)), np.zeros(len(t0))])
    times = np.unique(t[e == 1])
    O1 = E1 = V = 0.0
    for tt in times:
        at = t >= tt
        n = at.sum(); n1 = (at & (g == 1)).sum()
        d = ((t == tt) & (e == 1)).sum(); d1 = ((t == tt) & (e == 1) & (g == 1)).sum()
        if n < 2:
            continue
        O1 += d1; E1 += d * n1 / n
        V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    if V <= 0:
        return np.nan, np.nan
    chi = (O1 - E1) ** 2 / V
    return chi, 1 - chi2.cdf(chi, 1)


os_rows, lr_rows = [], []
ORGS = ["BONE", "LIVER", "CNS_BRAIN", "PLEURA", "ADRENAL_GLANDS", "INTRA_ABDOMINAL"]
for stratum in ["LUAD", "LUSC"]:
    d = OS[OS.grp == stratum]
    # organs
    for o in ORGS:
        r = cox(d, [o + "_b"], [c for c in SCOV])
        if r:
            v = r[o + "_b"]
            os_rows.append(dict(model="organ", exposure=o, stratum=stratum, **v))
        dd = d.dropna(subset=[o + "_b", "os_m", "os_event"])
        if len(dd) > 50:
            c_, p_ = logrank(dd.loc[dd[o + "_b"], "os_m"].values, dd.loc[dd[o + "_b"], "os_event"].values,
                             dd.loc[~dd[o + "_b"], "os_m"].values, dd.loc[~dd[o + "_b"], "os_event"].values)
            lr_rows.append(dict(exposure=o, stratum=stratum, logrank_p=p_, chi2=c_, n=len(dd)))
    # burden
    for label, col in [("n_org_cont", "n_distant"), ("n_org_ge3", None)]:
        if col is None:
            d2 = d.copy(); d2["n_org_ge3"] = (d2.n_distant >= 3).astype(int); col = "n_org_ge3"
            r = cox(d2, [col], SCOV[:-1])
        else:
            r = cox(d, [col], SCOV[:-1])
        if r:
            os_rows.append(dict(model="burden", exposure=label, stratum=stratum, **r[col]))
    # key genomic exposures (adjusted for burden)
    for col in ["M_STK11", "M_SMARCA4", "M_KEAP1", "M_EGFR", "M_KRAS", "M_TP53", "D_CDKN2A"]:
        if col not in d.columns:
            continue
        r = cox(d, [col], SCOV)
        if r:
            os_rows.append(dict(model="genomic", exposure=col, stratum=stratum, **r[col]))
ost = pd.DataFrame(os_rows)
ost.to_csv(os.path.join(OUT, "stage5_os_cox.csv"), index=False)
print("\n=== Cox（校正年龄/性别/吸烟/分期/负荷 等）")
print(ost.round(3).to_string(index=False))
lr = pd.DataFrame(lr_rows)
lr.to_csv(os.path.join(OUT, "stage5_os_logrank_handcoded.csv"), index=False)
print("\n=== 手算 log-rank（两组合并比较，与 Cox 对照）")
print(lr.round(4).to_string(index=False))
print("\nDONE")
