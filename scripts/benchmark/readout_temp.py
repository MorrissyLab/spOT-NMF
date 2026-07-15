"""Readout-temperature calibration of the native OT usages.

The whole precision->accuracy chain runs through ONE line of the model::

    W = softmin( (log k / w) * H^T G )          # inverse-temperature coef = log(k)/w

``w`` here does double duty: it is the entropy regularizer *during* the dual
ascent AND it sets the peakiness of the FINAL usage readout. But the peakiness
that minimizes the deconvolution error (vs graded ground-truth proportions) need
not equal the ``w`` that is best for the optimization. This script decouples
them, essentially for free: train ONCE per dataset (at its shipped config) to get
G and H, then recompute the usages at a grid of temperature multipliers ``tau``::

    W(tau) = softmin( tau * coef * H^T G )

and score each against ground truth. tau<1 softens (more graded, better for mixed
spots); tau>1 sharpens (more peaked, better for crisp panels). tau=1 is the
current behaviour. This isolates exactly how much accuracy the fixed training
temperature leaves on the table -- no retraining per tau.

Writes ``results/readout_temp_<ds>.csv`` (one row per tau).

Usage::
    uv run python scripts/benchmark/readout_temp.py
    uv run python scripts/benchmark/readout_temp.py --dataset MOB_dance_sim_norm_mm
"""
import argparse
import contextlib
import os

import anndata as ad
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

from benchmark_mosaic import match_and_score, REPO
from consensus_ablation import DATASETS, get_norm
from spotnmf.models import spotnmf, seed_all, _resolve_dtype

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# Per-dataset row-normalization (matches benchmark_all OVERRIDES: MOB -> False).
NORMROWS = {
    "Dataset10_STARmap_li2022_sim_norm_mm": True,
    "Dataset4_seqFISH_li2022_sim_norm_mm": True,
    "MOB_dance_sim_norm_mm": False,
}
TAUS = [0.25, 0.4, 0.5, 0.65, 0.8, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0]


def train_model(X, var_idx, obs_idx, k, w, seed, normalize_rows):
    """Train one native spOT-NMF model on the pre-loaded matrix; return the model."""
    sub = ad.AnnData(X=np.asarray(X, np.float64), obs=pd.DataFrame(index=obs_idx),
                     var=pd.DataFrame(index=var_idx))
    sub.var["highly_variable"] = True
    seed_all(seed)
    model = spotnmf(factors=int(k), h_regularization=1e-2, w_regularization=w,
                    eps=2e-2, cost="cosine", pca_cost=False)
    with contextlib.redirect_stdout(open(os.devnull, "w")), \
            contextlib.redirect_stderr(open(os.devnull, "w")):
        model.train(sub, lr=1e-2, optim_name="adam", tol_inner=1e-12,
                    tol_outer=1e-5, normalize_rows=normalize_rows, max_iter=50,
                    max_iter_inner=150, dtype=_resolve_dtype("float64"),
                    device=DEVICE)
    return model


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default=None, help="one dataset name; default all")
    ap.add_argument("--seeds", type=int, default=3)
    args = ap.parse_args()

    ds_list = [d for d in DATASETS if d[0] in NORMROWS]
    if args.dataset:
        ds_list = [d for d in ds_list if d[0] == args.dataset]

    for name, adata_path, n_hvf, w in ds_list:
        nr = NORMROWS[name]
        print(f"### readout-temp  {name}  (w={w}, normalize_rows={nr}, device={DEVICE})")
        # Load + normalize the matrix ONCE (seed-independent); only training varies.
        X, var_idx, obs_idx, gt_prop, k = get_norm(name, adata_path, n_hvf)
        cols = [f"ot_{i+1}" for i in range(int(k))]
        acc = {tau: [] for tau in TAUS}  # tau -> list of (pcc, rmse, jsd) per seed
        for seed in range(args.seeds):
            model = train_model(X, var_idx, obs_idx, k, w, seed, nr)
            # H^T G is tau-independent: compute once, then softmin per tau (free).
            coef = np.log(model.factors) / model.w_regularization
            htg = (model.H.T @ model.G).detach()
            for tau in TAUS:
                W = F.softmin(tau * coef * htg, dim=0).T.cpu().numpy()
                s = match_and_score(pd.DataFrame(W, index=obs_idx, columns=cols),
                                    gt_prop)
                acc[tau].append((s["PCC_mean"], s["RMSE"], s["JSD"]))
        rows = []
        for tau in TAUS:
            a = np.array(acc[tau])  # (seeds, 3)
            rows.append(dict(tau=tau,
                             PCC_mean=a[:, 0].mean(), PCC_std=a[:, 0].std(),
                             RMSE=a[:, 1].mean(), JSD=a[:, 2].mean()))
            star = "  <-- tau=1 (current)" if tau == 1.0 else ""
            print(f"  tau={tau:<5g}  PCC={a[:,0].mean():.4f}±{a[:,0].std():.4f}"
                  f"  RMSE={a[:,1].mean():.4f}  JSD={a[:,2].mean():.4f}{star}")
        best = max(rows, key=lambda r: r["PCC_mean"])
        base = next(r for r in rows if r["tau"] == 1.0)
        print(f"  >>> best tau={best['tau']:g}  PCC={best['PCC_mean']:.4f}"
              f"  (tau=1: {base['PCC_mean']:.4f}, delta {best['PCC_mean']-base['PCC_mean']:+.4f})")
        out = REPO / "results" / f"readout_temp_{name[:8]}.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(out, index=False)
        print(f"  wrote {out}\n")


if __name__ == "__main__":
    main()
