"""Phase 3 probe: can an UNSUPERVISED criterion pick the best consensus config?

The refit sweep showed the best ``(method, aggregate)`` differs by dataset
(STARmap->kmeans, seqFISH->match), so a fixed default forfeits ~0.02 PCC. A
selector must choose without ground truth. This probe tests two candidate
criteria against the (ground-truth) PCC ranking, on the cached replicate pools:

  * ``recon``  -- relative NNLS reconstruction error ||X - W H||_F / ||X||_F of
                  the refit (lower = spectra explain the data better).
  * ``self_pcc`` -- mean across-spot correlation between each consensus program's
                  NNLS usage and that program's *own* row in the pooled replicates
                  matched usage (a self-consistency proxy). [reconstruction only
                  here; kept minimal]

For each dataset we build every ``(method, aggregate)`` consensus, NNLS-refit,
and print PCC (oracle) alongside the recon error, then report whether argmin(recon)
== argmax(PCC).

Usage::  uv run python scripts/benchmark/consensus_select_probe.py
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import non_negative_factorization

from benchmark_mosaic import match_and_score, REPO
from consensus_ablation import DATASETS, replicate_pool
from spotnmf.models import _cluster_consensus_spectra

METHODS = [("kmeans", "median"), ("kmeans", "mean"), ("kmeans", "geomean"),
           ("match", "median"), ("match", "mean"), ("match", "geomean")]


def nnls_refit(H, X, k, seed=0):
    W, Hout, _ = non_negative_factorization(
        np.ascontiguousarray(X, np.float64), H=np.ascontiguousarray(H),
        n_components=k, init="custom", update_H=False, solver="cd",
        max_iter=500, random_state=seed)
    recon = np.linalg.norm(X - W @ Hout) / (np.linalg.norm(X) + 1e-12)
    return W, recon


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="*",
                    default=["Dataset10_STARmap_li2022_sim_norm_mm",
                             "Dataset4_seqFISH_li2022_sim_norm_mm",
                             "MOB_dance_sim_norm_mm"])
    ap.add_argument("--n_iter", type=int, default=12)
    args = ap.parse_args()

    want = set(args.datasets)
    selected = [d for d in DATASETS if d[0] in want]

    rows = []
    for name, adata_path, n_hvf, train_w in selected:
        if adata_path and not Path(adata_path).exists():
            continue
        pool, X, obs_idx, var_idx, gt_prop, k = replicate_pool(
            name, adata_path, n_hvf, train_w, args.n_iter)
        n_iter = pool.shape[0]
        spectra_pool = pd.DataFrame(
            pool.reshape(n_iter * k, -1), columns=var_idx,
            index=[f"iter{i}_topic{t + 1}" for i in range(n_iter) for t in range(k)])

        print(f"\n### {name}  (k={k})")
        recs = []
        for method, agg in METHODS:
            H, _ = _cluster_consensus_spectra(
                spectra_pool, k, 0.5, 0.30, n_iter=n_iter, method=method, aggregate=agg)
            W, recon = nnls_refit(H.to_numpy(), X, k)
            usage = pd.DataFrame(W, index=obs_idx, columns=[str(i + 1) for i in range(k)])
            s = match_and_score(usage, gt_prop)
            recs.append((method, agg, s["PCC_mean"], recon))
            rows.append(dict(dataset=name, method=method, aggregate=agg,
                             PCC=s["PCC_mean"], recon=recon))
            print(f"  {method:6s}+{agg:8s}  PCC={s['PCC_mean']:.4f}  recon={recon:.5f}")

        best_pcc = max(recs, key=lambda r: r[2])
        best_recon = min(recs, key=lambda r: r[3])
        match_ok = (best_pcc[0], best_pcc[1]) == (best_recon[0], best_recon[1])
        print(f"  -> oracle PCC : {best_pcc[0]}+{best_pcc[1]} ({best_pcc[2]:.4f})")
        print(f"  -> min recon  : {best_recon[0]}+{best_recon[1]} "
              f"(PCC there = {best_recon[2]:.4f})  {'MATCH' if match_ok else 'MISS'}")

    df = pd.DataFrame(rows)
    out = REPO / "results" / "consensus_select_probe.csv"
    df.to_csv(out, index=False)
    # Correlation of recon vs PCC within each dataset (want strongly negative).
    print("\nSpearman(recon, PCC) per dataset (want negative):")
    for name in df["dataset"].unique():
        sub = df[df["dataset"] == name]
        rho = sub["recon"].corr(sub["PCC"], method="spearman")
        print(f"  {name}: {rho:+.3f}")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
