"""Consensus lever: does increasing the replicate count (n_iter) help?

cNMF typically uses 100+ replicates; spOT-NMF's consensus uses 12. This factorizes
``max_n`` replicates ONCE, then evaluates the shipped consensus (match + Wasserstein
barycenter + NNLS refit) on the first {12, 18, 24, ...} replicates -- so the cost is
one factorization pass, not one per n_iter. Focused on seqFISH (the dataset closest
to cNMF) but works on any.

Usage::  uv run python scripts/benchmark/consensus_niter.py --dataset Dataset4_seqFISH_li2022_sim_norm_mm --max_n 24
"""
import argparse
import contextlib
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.decomposition import non_negative_factorization

from benchmark_mosaic import match_and_score, REPO
from consensus_ablation import DATASETS, get_norm
from spotnmf.models import (run_spotnmf, _cluster_consensus_spectra,
                            _gene_ground_cost)

import torch
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="Dataset4_seqFISH_li2022_sim_norm_mm")
    ap.add_argument("--max_n", type=int, default=24)
    ap.add_argument("--steps", nargs="*", type=int, default=[12, 18, 24])
    ap.add_argument("--aggregate", default="wasserstein")
    ap.add_argument("--bary_reg", type=float, default=0.1)
    args = ap.parse_args()

    name, adata_path, n_hvf, w = {d[0]: d for d in DATASETS}[args.dataset]
    X, var_idx, obs_idx, gt_prop, k = get_norm(name, adata_path, n_hvf)
    print(f"### {name}  k={k}  w={w}  factorizing {args.max_n} replicates...")

    sub0 = ad.AnnData(X=np.asarray(X, np.float64), obs=pd.DataFrame(index=obs_idx),
                      var=pd.DataFrame(index=var_idx))
    sub0.var["highly_variable"] = True
    progs = []
    for i in range(args.max_n):
        with contextlib.redirect_stdout(open(os.devnull, "w")), \
                contextlib.redirect_stderr(open(os.devnull, "w")):
            res, _ = run_spotnmf(sub0.copy(), components=k, seed=i, device=DEVICE,
                                 optim_name="adam", lr=0.01, max_iter=50,
                                 max_iter_inner=150, eps=0.02, h=1e-2, w=w,
                                 normalize_rows=True, cost="cosine")
        progs.append(res["genes_per_topic"].T.to_numpy())
        print(f"  replicate {i + 1}/{args.max_n}")
    pool = np.stack(progs)

    M = _gene_ground_cost(sub0, var_idx, "cosine", True) if args.aggregate == "wasserstein" else None

    rows = []
    for n in args.steps:
        if n > args.max_n:
            continue
        sp = pd.DataFrame(pool[:n].reshape(n * k, -1), columns=var_idx,
                          index=[f"i{i}_t{t}" for i in range(n) for t in range(k)])
        H, agr = _cluster_consensus_spectra(sp, k, 0.5, 0.30, n_iter=n,
                                            method="match", aggregate=args.aggregate,
                                            cost_matrix=M, bary_reg=args.bary_reg)
        W, _, _ = non_negative_factorization(
            np.ascontiguousarray(X, np.float64), H=np.ascontiguousarray(H.to_numpy()),
            n_components=k, init="custom", update_H=False, solver="cd",
            max_iter=500, random_state=0)
        s = match_and_score(pd.DataFrame(W, index=obs_idx,
                            columns=[str(i + 1) for i in range(k)]), gt_prop)
        rows.append(dict(n_iter=n, PCC=s["PCC_mean"], RMSE=s["RMSE"],
                         JSD=s["JSD"], agreement=agr))
        print(f"  n_iter={n:3d}  PCC={s['PCC_mean']:.4f} RMSE={s['RMSE']:.4f} "
              f"JSD={s['JSD']:.4f}  agreement={agr:.3f}")

    out = REPO / "results" / f"consensus_niter_{name[:8]}.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
