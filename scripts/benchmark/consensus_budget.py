"""Budget-matched accuracy: float32 + more replicates vs float64 consensus.

The precision study showed float32's accuracy cost is noise in G that the
consensus barycenter averages out. float32+analytic is ~2x faster per replicate,
so at a fixed wall-clock budget we can afford many more replicates. Question:
does float32 with more replicates match or BEAT float64 consensus at equal (or
less) wall-clock? If so, the speed win converts directly into accuracy.

Reference: float64 n_iter=10 (from consensus_dtype): STARmap 0.5995/472s,
seqFISH 0.4757/384s, MOB 0.8366/565s.

Usage::
    uv run python scripts/benchmark/consensus_budget.py
"""
import contextlib
import os
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch

from benchmark_mosaic import match_and_score, REPO
from consensus_ablation import DATASETS, get_norm
from spotnmf.models import run_spotnmf_consensus

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
AGG = {
    "Dataset10_STARmap_li2022_sim_norm_mm": "wasserstein",
    "Dataset4_seqFISH_li2022_sim_norm_mm": "wasserstein",
    "MOB_dance_sim_norm_mm": "geomean",
}
# (dtype, n_iter) points: float64 baseline vs float32 at increasing replicates.
POINTS = [("float64", 10), ("float32", 10), ("float32", 20), ("float32", 30)]


def score(name, path, nhvf, w, dtype, n_iter):
    X, vi, oi, gt, k = get_norm(name, path, nhvf)
    sub = ad.AnnData(X=np.asarray(X, np.float64), obs=pd.DataFrame(index=oi),
                     var=pd.DataFrame(index=vi))
    sub.var["highly_variable"] = True
    agg = AGG.get(name, "geomean")
    t0 = time.time()
    with open(os.devnull, "w") as nul, contextlib.redirect_stdout(nul), \
            contextlib.redirect_stderr(nul):
        res, _ = run_spotnmf_consensus(
            sub, components=k, n_iter=n_iter, seed=0, device=DEVICE,
            consensus_method="match", aggregate=agg, refit="nnls",
            optim_name="adam", lr=0.01, max_iter=50, max_iter_inner=150,
            eps=0.02, h=1e-2, w=w, normalize_rows=True, cost="cosine", dtype=dtype)
    dt = time.time() - t0
    s = match_and_score(res["topics_per_spot"], gt)
    return s["PCC_mean"], s["RMSE"], dt, agg


def main():
    rows = []
    for name, path, nhvf, w in [d for d in DATASETS if not d[0].startswith("stereoseq")]:
        print(f"\n### {name}  (w={w})")
        for dtype, n_iter in POINTS:
            pcc, rmse, dt, agg = score(name, path, nhvf, w, dtype, n_iter)
            rows.append(dict(dataset=name, dtype=dtype, n_iter=n_iter, aggregate=agg,
                             PCC=pcc, RMSE=rmse, sec=dt))
            print(f"  {dtype:8s} n_iter={n_iter:3d} ({agg:11s})  "
                  f"PCC={pcc:.4f}  RMSE={rmse:.4f}  {dt:.0f}s")
    df = pd.DataFrame(rows)
    out = str(REPO / "results" / "consensus_budget.csv")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
