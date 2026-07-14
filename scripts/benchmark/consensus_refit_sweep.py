"""Phase 0-1: refit strategy for consensus spOT-NMF (NNLS vs OT vs OT+refit-w).

Audit finding #1: ``run_spotnmf_consensus`` can only refit usages with the OT
solve, but on mixed-spot data an NNLS refit of the *same* consensus spectra
deconvolves better (OT usages are too peaked). This script quantifies the effect
of the refit *choice* -- and of a **decoupled refit-w** for the OT path -- while
holding the consensus spectra fixed, so the refit is isolated from the
aggregation.

Protocol (reuses the cached replicate pools from ``consensus_ablation.py``, so no
re-factorization):

  1. For each dataset, load the ``(n_iter, k, genes)`` replicate pool + the
     mosaicMPI-normalized matrix ``X`` + ground-truth proportions.
  2. For every ``(method, aggregate)`` consensus, build the consensus spectra once.
  3. Score those spectra with:
       * NNLS refit  (scikit-learn non_negative_factorization, update_H=False)
       * OT refit    (spotnmf.transform) at each refit-w in ``--refit_w``.
  4. Also score the single-run (no-consensus) baseline under both refits.

Writes ``results/consensus_refit_sweep.csv`` with one row per
``(dataset, method, aggregate, refit, refit_w)``.

Usage::

    uv run python scripts/benchmark/consensus_refit_sweep.py
    uv run python scripts/benchmark/consensus_refit_sweep.py --datasets MOB_dance_sim_norm_mm
"""
import argparse
import contextlib
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch

from benchmark_mosaic import match_and_score, REPO
from consensus_ablation import DATASETS, replicate_pool, nnls_score
from spotnmf.models import _cluster_consensus_spectra, spotnmf

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# Consensus variants to test (clustering x aggregation).
METHODS = [("kmeans", "median"), ("kmeans", "mean"), ("kmeans", "geomean"),
           ("match", "median"), ("match", "mean"), ("match", "geomean")]


@contextlib.contextmanager
def _quiet():
    with open(os.devnull, "w") as nul, contextlib.redirect_stdout(nul), \
            contextlib.redirect_stderr(nul):
        yield


def ot_score(H, X, obs_idx, var_idx, gt_prop, k, w, train_w):
    """Score consensus spectra ``H`` (k x genes) via the OT usage refit.

    ``w`` is the *refit* usage-entropy (decoupled from the replicate training w).
    Everything else matches the tuned spOT-NMF defaults.
    """
    sub = ad.AnnData(X=np.asarray(X, np.float64),
                     obs=pd.DataFrame(index=obs_idx),
                     var=pd.DataFrame(index=var_idx))
    sub.var["highly_variable"] = True
    model = spotnmf(factors=k, h_regularization=1e-2, w_regularization=w,
                    eps=2e-2, cost="cosine")
    with _quiet():
        model.transform(sub, spectra=np.ascontiguousarray(H),
                        lr=1e-2, optim_name="adam", tol_inner=1e-12,
                        tol_outer=1e-5, normalize_rows=True,
                        max_iter=50, max_iter_inner=150, device=DEVICE)
    usage = pd.DataFrame(sub.obsm["W_OT"], index=obs_idx,
                         columns=[str(i + 1) for i in range(k)])
    s = match_and_score(usage, gt_prop)
    return s["PCC_mean"], s["RMSE"], s["JSD"], s["usage_dispersion"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="*",
                    default=["Dataset10_STARmap_li2022_sim_norm_mm",
                             "Dataset4_seqFISH_li2022_sim_norm_mm",
                             "MOB_dance_sim_norm_mm"])
    ap.add_argument("--n_iter", type=int, default=12)
    ap.add_argument("--refit_w", nargs="*", type=float,
                    default=[0.01, 0.02, 0.05, 0.1, 0.2])
    ap.add_argument("--out", default=str(REPO / "results" / "consensus_refit_sweep.csv"))
    args = ap.parse_args()

    want = set(args.datasets)
    selected = [d for d in DATASETS if d[0] in want]

    rows = []
    for name, adata_path, n_hvf, train_w in selected:
        if adata_path and not Path(adata_path).exists():
            print(f"SKIP {name}: {adata_path} not found")
            continue
        print(f"\n### {name}  (train_w={train_w})")
        pool, X, obs_idx, var_idx, gt_prop, k = replicate_pool(
            name, adata_path, n_hvf, train_w, args.n_iter)
        n_iter = pool.shape[0]
        spectra_pool = pd.DataFrame(
            pool.reshape(n_iter * k, -1), columns=var_idx,
            index=[f"iter{i}_topic{t + 1}" for i in range(n_iter) for t in range(k)])

        def record(method, agg, H, agreement):
            p, r, j = nnls_score(H, X, obs_idx, gt_prop, k)
            rows.append(dict(dataset=name, k=k, method=method, aggregate=agg,
                             refit="nnls", refit_w=np.nan, PCC=p, RMSE=r, JSD=j,
                             dispersion=np.nan, agreement=agreement))
            print(f"  {method:6s}+{agg:8s} nnls           PCC={p:.3f} RMSE={r:.3f} JSD={j:.3f}")
            for w in args.refit_w:
                p, r, j, disp = ot_score(H, X, obs_idx, var_idx, gt_prop, k, w, train_w)
                rows.append(dict(dataset=name, k=k, method=method, aggregate=agg,
                                 refit="ot", refit_w=w, PCC=p, RMSE=r, JSD=j,
                                 dispersion=disp, agreement=agreement))
                print(f"  {method:6s}+{agg:8s} ot   w={w:<5g}  PCC={p:.3f} RMSE={r:.3f} "
                      f"JSD={j:.3f} disp={disp:.3f}")

        # Single-run (no consensus) baseline, both refits.
        H0 = pool[0] / pool[0].sum(1, keepdims=True)
        record("single-run", "-", H0, np.nan)

        # Consensus variants.
        for method, agg in METHODS:
            H, agreement = _cluster_consensus_spectra(
                spectra_pool, k, 0.5, 0.30, n_iter=n_iter, method=method, aggregate=agg)
            record(method, agg, H.to_numpy(), agreement)

    df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\nWrote {args.out}")

    # Per-dataset best config (by PCC) vs the shipped consensus-OT baseline.
    print("\n" + "=" * 78)
    print("BEST CONFIG PER DATASET (by PCC)")
    print("=" * 78)
    for name in df["dataset"].unique():
        sub = df[df["dataset"] == name]
        best = sub.loc[sub["PCC"].idxmax()]
        print(f"{name}\n  best: {best['method']}+{best['aggregate']} "
              f"refit={best['refit']} w={best['refit_w']}  "
              f"PCC={best['PCC']:.3f} RMSE={best['RMSE']:.3f} JSD={best['JSD']:.3f}")


if __name__ == "__main__":
    main()
