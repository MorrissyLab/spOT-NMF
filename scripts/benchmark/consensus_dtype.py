"""Consensus spOT-NMF in float32 vs float64.

Tests the hypothesis that the consensus barycenter washes out float32's extra
per-replicate rounding noise, so the (n_iter) replicate factorizations can run in
float32 -- 3-6x faster each on GPUs with weak FP64 -- at no consensus-accuracy
cost. Only the replicate factorizations are affected by ``dtype``; the Wasserstein
gene cost and the NNLS usage refit stay float64 regardless.

Usage::
    uv run python scripts/benchmark/consensus_dtype.py --n_iter 10
"""
import argparse
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

# Established best aggregation per dataset (see experiments doc).
AGG = {
    "Dataset10_STARmap_li2022_sim_norm_mm": "wasserstein",
    "Dataset4_seqFISH_li2022_sim_norm_mm": "wasserstein",
    "MOB_dance_sim_norm_mm": "geomean",
}


def score_consensus(name, adata_path, n_hvf, w, n_iter, dtype):
    X, var_idx, obs_idx, gt_prop, k = get_norm(name, adata_path, n_hvf)
    sub = ad.AnnData(X=np.asarray(X, np.float64), obs=pd.DataFrame(index=obs_idx),
                     var=pd.DataFrame(index=var_idx))
    sub.var["highly_variable"] = True
    agg = AGG.get(name, "geomean")
    if torch.cuda.is_available():
        torch.cuda.empty_cache(); torch.cuda.reset_peak_memory_stats()
    t0 = time.time()
    with open(os.devnull, "w") as nul, contextlib.redirect_stdout(nul), \
            contextlib.redirect_stderr(nul):
        res, _ = run_spotnmf_consensus(
            sub, components=k, n_iter=n_iter, seed=0, device=DEVICE,
            consensus_method="match", aggregate=agg, refit="nnls",
            optim_name="adam", lr=0.01, max_iter=50, max_iter_inner=150,
            eps=0.02, h=1e-2, w=w, normalize_rows=True, cost="cosine",
            dtype=dtype)
    dt = time.time() - t0
    peak = torch.cuda.max_memory_allocated() / 1e6 if torch.cuda.is_available() else float("nan")
    s = match_and_score(res["topics_per_spot"], gt_prop)
    return s["PCC_mean"], s["RMSE"], s["JSD"], dt, peak, agg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="*", default=None)
    ap.add_argument("--n_iter", type=int, default=10)
    ap.add_argument("--out", default=str(REPO / "results" / "consensus_dtype.csv"))
    args = ap.parse_args()

    selected = [d for d in DATASETS if not d[0].startswith("stereoseq")]
    if args.datasets:
        want = set(args.datasets)
        selected = [d for d in DATASETS if d[0] in want]

    rows = []
    for name, adata_path, n_hvf, w in selected:
        print(f"\n### {name}  (w={w}, n_iter={args.n_iter})")
        for dtype in ["float64", "float32"]:
            pcc, rmse, jsd, dt, peak, agg = score_consensus(
                name, adata_path, n_hvf, w, args.n_iter, dtype)
            rows.append(dict(dataset=name, dtype=dtype, aggregate=agg, PCC=pcc,
                             RMSE=rmse, JSD=jsd, sec=dt, peak_mb=peak))
            print(f"  {dtype:8s} ({agg:11s})  PCC={pcc:.4f}  RMSE={rmse:.4f}  "
                  f"{dt:.1f}s  peak={peak:.0f}MB")

    df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
