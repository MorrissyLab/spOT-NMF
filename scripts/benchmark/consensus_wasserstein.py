"""Phase 2: entropic Wasserstein-barycenter aggregation for consensus spOT-NMF.

spOT-NMF spectra are probability distributions over genes with an optimal-transport
ground cost ``M`` between genes (cosine on the normalized reference, exactly as the
model builds it). The current consensus aggregates matched replicate programs with
a *mean* / *geomean* / *median* -- none of which use that gene geometry. The
OT-native consensus of gene-distributions is their **entropic Wasserstein
barycenter** under ``M`` (Agueh-Carlier; computed with Benamou et al.'s Iterative
Bregman Projections). This script tests whether it beats the barycenter-free
aggregates, especially on seqFISH (the one dataset where consensus spOT-NMF still
trails cNMF).

Protocol (reuses the cached replicate pools -- no re-factorization):
  1. Match every replicate's k programs to a medoid reference (same Hungarian
     matching as the package ``match`` consensus).
  2. For each cluster, aggregate the n_iter aligned gene-distributions with
     {mean, geomean, wasserstein(reg)} and NNLS-refit + score vs ground truth.

Writes ``results/consensus_wasserstein.csv``.

Usage::  uv run python scripts/benchmark/consensus_wasserstein.py
         uv run python scripts/benchmark/consensus_wasserstein.py --datasets Dataset4_seqFISH_li2022_sim_norm_mm
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from sklearn.decomposition import non_negative_factorization

from benchmark_mosaic import match_and_score, REPO
from consensus_ablation import DATASETS, replicate_pool
from spotnmf.models import _cross_replicate_agreement, _l1_normalize_rows


def gene_cost(X, normalize_rows=True):
    """Genes x genes cosine cost on the model's normalized reference (from X)."""
    A = np.asarray(X, np.float64).T.copy()          # genes x spots
    A += 1e-6
    if normalize_rows:
        mean_row_sum = A.sum(1).mean()
        A /= A.sum(1, keepdims=True) * mean_row_sum
    A /= A.sum(0, keepdims=True)
    features = 1e-6 + A
    M = cdist(features, features, metric="cosine")
    return M / (M.max() + 1e-12)


def geometric_bar(weights, alpha):
    return np.exp(np.log(np.clip(alpha, 1e-300, None)) @ weights)


def wass_barycenter(A, K, weights=None, numItermax=1000, stopThr=1e-5):
    """Entropic Wasserstein barycenter (IBP) of columns of ``A`` (n x S), kernel K."""
    n, S = A.shape
    if weights is None:
        weights = np.ones(S) / S
    u = np.ones((n, S)) / n
    UKv = u * (K @ (A / np.clip(K @ u, 1e-300, None)))
    bar = geometric_bar(weights, UKv)
    for cpt in range(numItermax):
        u = (bar[:, None] / np.clip(UKv, 1e-300, None)) * u
        UKv = u * (K @ (A / np.clip(K @ u, 1e-300, None)))
        bar_new = geometric_bar(weights, UKv)
        if cpt % 10 == 0 and np.abs(bar_new - bar).sum() < stopThr:
            bar = bar_new
            break
        bar = bar_new
    s = bar.sum()
    return bar / s if s > 0 else np.ones(n) / n


def aligned_clusters(pool):
    """Return ``(k, S, genes)`` matched distributions: [cluster][replicate]."""
    _, ref, units, dists = _cross_replicate_agreement(pool)
    k = pool.shape[1]
    aligned = [dists[ref]]
    for j in range(len(dists)):
        if j == ref:
            continue
        C = units[ref] @ units[j].T
        r, c = linear_sum_assignment(-C)
        order = np.empty(k, dtype=int)
        order[r] = c
        aligned.append(dists[j][order])
    stack = np.stack(aligned)          # (S, k, genes)
    return np.transpose(stack, (1, 0, 2))   # (k, S, genes)


def nnls_score(H, X, obs_idx, gt_prop, k, seed=0):
    W, _, _ = non_negative_factorization(
        np.ascontiguousarray(X, np.float64), H=np.ascontiguousarray(H),
        n_components=k, init="custom", update_H=False, solver="cd",
        max_iter=500, random_state=seed)
    usage = pd.DataFrame(W, index=obs_idx, columns=[str(i + 1) for i in range(k)])
    s = match_and_score(usage, gt_prop)
    return s["PCC_mean"], s["RMSE"], s["JSD"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="*",
                    default=["Dataset10_STARmap_li2022_sim_norm_mm",
                             "Dataset4_seqFISH_li2022_sim_norm_mm",
                             "MOB_dance_sim_norm_mm"])
    ap.add_argument("--n_iter", type=int, default=12)
    ap.add_argument("--regs", nargs="*", type=float, default=[0.01, 0.02, 0.05, 0.1])
    ap.add_argument("--out", default=str(REPO / "results" / "consensus_wasserstein.csv"))
    args = ap.parse_args()

    want = set(args.datasets)
    selected = [d for d in DATASETS if d[0] in want]

    rows = []
    for name, adata_path, n_hvf, train_w in selected:
        if adata_path and not Path(adata_path).exists():
            continue
        pool, X, obs_idx, var_idx, gt_prop, k = replicate_pool(
            name, adata_path, n_hvf, train_w, args.n_iter)
        print(f"\n### {name}  (k={k}, genes={pool.shape[2]})")
        clusters = aligned_clusters(pool)          # (k, S, genes)

        # Baselines: match + {mean, geomean} (barycenter-free).
        H_mean = _l1_normalize_rows(clusters.mean(axis=1))
        H_geo = _l1_normalize_rows(np.exp(np.log(np.clip(clusters, 1e-12, None)).mean(axis=1)))
        for tag, H in [("mean", H_mean), ("geomean", H_geo)]:
            p, r, j = nnls_score(H, X, obs_idx, gt_prop, k)
            rows.append(dict(dataset=name, agg=tag, reg=np.nan, PCC=p, RMSE=r, JSD=j))
            print(f"  match+{tag:10s}          PCC={p:.4f} RMSE={r:.4f} JSD={j:.4f}")

        # Wasserstein barycenter under the gene cost, swept reg.
        M = gene_cost(X)
        for reg in args.regs:
            K = np.exp(-M / reg)
            H_w = np.stack([wass_barycenter(clusters[c].T, K) for c in range(k)])
            H_w = _l1_normalize_rows(H_w)
            p, r, j = nnls_score(H_w, X, obs_idx, gt_prop, k)
            rows.append(dict(dataset=name, agg="wasserstein", reg=reg, PCC=p, RMSE=r, JSD=j))
            print(f"  match+wasserstein reg={reg:<5g}  PCC={p:.4f} RMSE={r:.4f} JSD={j:.4f}")

    df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\nWrote {args.out}")
    print("\nBest per dataset:")
    for name in df["dataset"].unique():
        sub = df[df["dataset"] == name]
        b = sub.loc[sub["PCC"].idxmax()]
        print(f"  {name}: {b['agg']} reg={b['reg']}  PCC={b['PCC']:.4f}")


if __name__ == "__main__":
    main()
