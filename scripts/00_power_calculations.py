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
"""Minimal detectable OR (80% power, two-sided) for LUSC/LUAD organ phenotypes,
with alpha after BH over ~300 gene tests (alpha~1.7e-4) and nominal 0.05."""
import numpy as np
from scipy.stats import norm

def mdo(n1, n0, q, alpha):
    """smallest OR detectable for carrier freq q at given alpha, 80% power.
    Solve by scanning OR; carrier rate in cases = q*OR/(1+q*(OR-1)) approx via
    case-control sampling with prevalence ignored (balanced-ish)."""
    z_a = norm.ppf(1 - alpha / 2); z_b = norm.ppf(0.8)
    best = None
    for OR in np.arange(1.05, 8.0, 0.01):
        p1 = q * OR / (1 + q * (OR - 1))   # carriers enriched among cases
        p0 = q
        pbar = (n1 * p1 + n0 * p0) / (n1 + n0)
        se0 = np.sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n0))
        se1 = np.sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
        if (abs(p1 - p0) - z_a * se0) / se1 >= z_b:
            best = OR; break
    return best

for label, n1, n0 in [("LUAD BONE", 2851, 2836), ("LUAD CNS", 2154, 3533),
                      ("LUSC BONE", 286, 513), ("LUSC CNS", 164, 635),
                      ("LUSC PLEURA", 226, 573), ("LUSC LIVER", 169, 630)]:
    row = []
    for q in (0.05, 0.10, 0.20, 0.30):
        row.append(f"q={q:.2f}: OR>={mdo(n1, n0, q, 0.05/300):.2f} (nominal 0.05: OR>={mdo(n1,n0,q,0.05):.2f})")
    print(f"{label:12s} n1={n1:5d} n0={n0:5d} | " + " | ".join(row))
