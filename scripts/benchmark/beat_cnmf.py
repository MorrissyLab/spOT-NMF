"""Can spOT-NMF's consensus beat cNMF on the two datasets where cNMF leads
(STARmap, seqFISH)?

Since every consensus config here ends in the SAME NNLS usage refit, the only
thing that separates cNMF from spOT-NMF-consensus is the quality of the consensus
SPECTRA. So we factorize ``max_n`` OT replicates ONCE per dataset, then sweep the
consensus knobs cheaply on the cached pool:

  * n_iter        -- replicate count (prefixes of the pool)
  * (method, aggregate) -- clustering + distribution barycenter

each scored with the shipped NNLS refit against ground truth. Reports the best
config per dataset vs the cNMF reference from the full benchmark.

Usage::
    uv run python scripts/benchmark/beat_cnmf.py
    uv run python scripts/benchmark/beat_cnmf.py --datasets Dataset10_STARmap_li2022_sim_norm_mm --hreg 0.01 0.03
"""
import argparse
import contextlib
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
from sklearn.decomposition import non_negative_factorization

from benchmark_mosaic import match_and_score, REPO
from consensus_ablation import DATASETS, get_norm
from spotnmf.models import (run_spotnmf, _cluster_consensus_spectra,
                            _gene_ground_cost)

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# cNMF reference PCC from the full best-config benchmark (results/benchmark_all_bestcfg).
CNMF_REF = {
    "Dataset10_STARmap_li2022_sim_norm_mm": 0.611,
    "Dataset4_seqFISH_li2022_sim_norm_mm": 0.460,
}
CONFIGS = [("match", "wasserstein"), ("match", "mean"), ("match", "geomean"),
           ("kmeans", "mean"), ("kmeans", "median")]


def factorize_pool(X, var_idx, obs_idx, k, w, h, max_n):
    """Factorize ``max_n`` OT replicates once; return (n, k, genes) spectra pool."""
    sub = ad.AnnData(X=np.asarray(X, np.float64), obs=pd.DataFrame(index=obs_idx),
                     var=pd.DataFrame(index=var_idx))
    sub.var["highly_variable"] = True
    progs = []
    for i in range(max_n):
        with contextlib.redirect_stdout(open(os.devnull, "w")), \
                contextlib.redirect_stderr(open(os.devnull, "w")):
            res, _ = run_spotnmf(sub.copy(), components=k, seed=i, device=DEVICE,
                                 optim_name="adam", lr=0.01, max_iter=50,
                                 max_iter_inner=150, eps=0.02, h=h, w=w,
                                 normalize_rows=True, cost="cosine")
        progs.append(res["genes_per_topic"].T.to_numpy())
    return np.stack(progs)


def score_consensus(pool, n, k, X, var_idx, obs_idx, gt_prop, method, aggregate, M):
    """Aggregate the first n replicates and NNLS-refit; return (PCC, RMSE, agr)."""
    sp = pd.DataFrame(pool[:n].reshape(n * k, -1), columns=var_idx,
                      index=[f"i{i}_t{t}" for i in range(n) for t in range(k)])
    H, agr = _cluster_consensus_spectra(
        sp, k, 0.5, 0.30, n_iter=n, method=method, aggregate=aggregate,
        cost_matrix=M if aggregate == "wasserstein" else None, bary_reg=0.1)
    W, _, _ = non_negative_factorization(
        np.ascontiguousarray(X, np.float64),
        H=np.ascontiguousarray(H.to_numpy()), n_components=k, init="custom",
        update_H=False, solver="cd", max_iter=500, random_state=0)
    s = match_and_score(pd.DataFrame(W, index=obs_idx,
                        columns=[str(i + 1) for i in range(k)]), gt_prop)
    return s["PCC_mean"], s["RMSE"], agr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="*",
                    default=["Dataset10_STARmap_li2022_sim_norm_mm",
                             "Dataset4_seqFISH_li2022_sim_norm_mm"])
    ap.add_argument("--max_n", type=int, default=30)
    ap.add_argument("--steps", nargs="*", type=int, default=[12, 18, 24, 30])
    ap.add_argument("--hreg", nargs="*", type=float, default=[0.01])
    args = ap.parse_args()

    for name in args.datasets:
        _, adata_path, n_hvf, w = {d[0]: d for d in DATASETS}[name]
        X, var_idx, obs_idx, gt_prop, k = get_norm(name, adata_path, n_hvf)
        ref = CNMF_REF.get(name, float("nan"))
        sub_cost = ad.AnnData(X=np.asarray(X, np.float64),
                              obs=pd.DataFrame(index=obs_idx),
                              var=pd.DataFrame(index=var_idx))
        sub_cost.var["highly_variable"] = True
        M = _gene_ground_cost(sub_cost, var_idx, "cosine", True)
        rows = []
        for h in args.hreg:
            print(f"### {name}  k={k} w={w} h={h}  (cNMF ref PCC={ref:.3f})  "
                  f"factorizing {args.max_n} replicates...")
            pool = factorize_pool(X, var_idx, obs_idx, k, w, h, args.max_n)
            for method, aggregate in CONFIGS:
                for n in args.steps:
                    if n > args.max_n:
                        continue
                    try:
                        pcc, rmse, agr = score_consensus(
                            pool, n, k, X, var_idx, obs_idx, gt_prop,
                            method, aggregate, M)
                    except Exception as e:
                        print(f"  {method:6s}+{aggregate:11s} n={n:2d}  FAILED: {e}")
                        continue
                    flag = "  *BEATS cNMF*" if pcc > ref else ""
                    rows.append(dict(h=h, method=method, aggregate=aggregate,
                                     n_iter=n, PCC=pcc, RMSE=rmse, agreement=agr))
                    print(f"  h={h} {method:6s}+{aggregate:11s} n={n:2d}  "
                          f"PCC={pcc:.4f} RMSE={rmse:.4f}{flag}")
        df = pd.DataFrame(rows).sort_values("PCC", ascending=False)
        best = df.iloc[0]
        print(f"\n>>> {name}: best spOT-NMF-consensus PCC={best['PCC']:.4f} "
              f"(h={best['h']} {best['method']}+{best['aggregate']} n={best['n_iter']:.0f})"
              f"  vs cNMF {ref:.3f}  -> {'WIN' if best['PCC'] > ref else 'still behind'}\n")
        out = REPO / "results" / f"beat_cnmf_{name[:8]}.csv"
        df.to_csv(out, index=False)
        print(f"Wrote {out}\n")


if __name__ == "__main__":
    main()
