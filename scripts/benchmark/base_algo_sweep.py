"""Deep-dive experiments on the BASE spOT-NMF optimizer (single run, native OT usages).

Answers, empirically, the audit questions about the core optimizer:
  * inner/outer iteration counts -- where does PCC plateau? are we over/under-iterating?
  * optimizer -- adam (lr 1e-2) vs lbfgs (lr 1.0, properly scaled) for the dual G.
  * convergence -- does the outer objective (model.losses) actually flatten?

All runs score the *native* OT per-spot usages (no consensus, no refit) against
ground-truth proportions, isolating raw algorithm quality. Uses the mosaicMPI-
normalized matrix (same input the benchmark feeds), on one dataset at a time
(default STARmap ds10: 189x882, fast, full panel).

Usage::
    uv run python scripts/benchmark/base_algo_sweep.py --exp iters
    uv run python scripts/benchmark/base_algo_sweep.py --exp optim
    uv run python scripts/benchmark/base_algo_sweep.py --dataset MOB_dance_sim_norm_mm --exp iters
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
from spotnmf.models import run_spotnmf

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


def _quiet():
    return contextlib.redirect_stdout(open(os.devnull, "w"))


def _apply_transform(X, kind):
    """Monotone, dynamic-range-compressing transforms applied before the model's
    own row/column normalization. All keep non-negativity."""
    X = np.asarray(X, np.float64)
    if kind in (None, "none"):
        return X
    if kind == "log1p":
        return np.log1p(X)
    if kind == "sqrt":
        return np.sqrt(np.clip(X, 0, None))
    if kind == "arcsinh":
        return np.arcsinh(X)
    if kind == "cbrt":
        return np.cbrt(np.clip(X, 0, None))
    raise ValueError(f"unknown transform {kind!r}")


def score_native(name, adata_path, n_hvf, w, seed, eps_override=0.02,
                 track_mem=False, x_transform=None, **run_kw):
    X, var_idx, obs_idx, gt_prop, k = get_norm(name, adata_path, n_hvf)
    X = _apply_transform(X, x_transform)
    sub = ad.AnnData(X=np.asarray(X, np.float64), obs=pd.DataFrame(index=obs_idx),
                     var=pd.DataFrame(index=var_idx))
    sub.var["highly_variable"] = True
    # Defaults; any key in run_kw overrides them (no duplicate-keyword clash).
    params = dict(device=DEVICE, eps=eps_override, h=1e-2, w=w,
                  normalize_rows=True, cost="cosine")
    params.update(run_kw)
    if track_mem and torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.reset_peak_memory_stats()
    t0 = time.time()
    with contextlib.redirect_stdout(open(os.devnull, "w")), \
            contextlib.redirect_stderr(open(os.devnull, "w")):
        res, losses = run_spotnmf(sub, components=k, seed=seed, **params)
    dt = time.time() - t0
    peak_mb = (torch.cuda.max_memory_allocated() / 1e6
               if track_mem and torch.cuda.is_available() else float("nan"))
    s = match_and_score(res["topics_per_spot"], gt_prop)
    return s["PCC_mean"], s["RMSE"], s["JSD"], dt, len(losses), float(losses[-1]), peak_mb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="Dataset10_STARmap_li2022_sim_norm_mm")
    ap.add_argument("--exp", choices=["iters", "optim", "seeds", "init", "damping",
                                      "eps", "normrows", "cost", "hreg", "lropt",
                                      "persist", "anneal", "dtype", "agrad",
                                      "transform", "rowscale", "batch", "gtail"], default="iters")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    ds = {d[0]: d for d in DATASETS}[args.dataset]
    name, adata_path, n_hvf, w = ds
    print(f"### base spOT-NMF  {name}  (w={w}, device={DEVICE})  exp={args.exp}")

    rows = []
    if args.exp == "iters":
        # Grid over inner/outer with adam. Multiple seeds to see variance.
        grid = [(10, 30), (30, 30), (60, 50), (150, 50), (150, 100), (400, 50)]
        for inner, outer in grid:
            pccs, times = [], []
            for seed in range(args.seeds):
                pcc, rmse, jsd, dt, nl, last, _ = score_native(
                    name, adata_path, n_hvf, w, seed, optim_name="adam", lr=1e-2,
                    max_iter=outer, max_iter_inner=inner)
                pccs.append(pcc); times.append(dt)
            rows.append(dict(inner=inner, outer=outer, optim="adam",
                             PCC_mean=np.mean(pccs), PCC_std=np.std(pccs),
                             sec=np.mean(times)))
            print(f"  inner={inner:4d} outer={outer:4d} adam   "
                  f"PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}  {np.mean(times):.0f}s")
    elif args.exp == "optim":
        for optim, lr, inner, outer in [
                ("adam", 1e-2, 150, 50), ("adam", 1e-2, 400, 100),
                ("lbfgs", 1.0, 50, 50), ("lbfgs", 1.0, 150, 50),
                ("lbfgs", 0.5, 150, 50)]:
            pccs, times, lasts = [], [], []
            for seed in range(args.seeds):
                pcc, rmse, jsd, dt, nl, last, _ = score_native(
                    name, adata_path, n_hvf, w, seed, optim_name=optim, lr=lr,
                    max_iter=outer, max_iter_inner=inner)
                pccs.append(pcc); times.append(dt); lasts.append(last)
            rows.append(dict(optim=optim, lr=lr, inner=inner, outer=outer,
                             PCC_mean=np.mean(pccs), PCC_std=np.std(pccs),
                             sec=np.mean(times), last_loss=np.mean(lasts)))
            print(f"  {optim:5s} lr={lr:<5g} inner={inner:4d} outer={outer:3d}  "
                  f"PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}  {np.mean(times):.0f}s  "
                  f"loss={np.mean(lasts):.4g}")
    elif args.exp == "normrows":
        for nr in [True, False]:
            pccs = []
            for seed in range(args.seeds):
                pcc, *_ = score_native(name, adata_path, n_hvf, w, seed,
                                       optim_name="adam", lr=1e-2, max_iter=50,
                                       max_iter_inner=150, normalize_rows=nr)
                pccs.append(pcc)
            rows.append(dict(normalize_rows=nr, PCC_mean=np.mean(pccs), PCC_std=np.std(pccs)))
            print(f"  normalize_rows={str(nr):5s}  PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}")
    elif args.exp == "cost":
        for cost in ["cosine", "correlation", "euclidean", "sqeuclidean", "cityblock"]:
            pccs = []
            for seed in range(args.seeds):
                pcc, *_ = score_native(name, adata_path, n_hvf, w, seed,
                                       optim_name="adam", lr=1e-2, max_iter=50,
                                       max_iter_inner=150, cost=cost)
                pccs.append(pcc)
            rows.append(dict(cost=cost, PCC_mean=np.mean(pccs), PCC_std=np.std(pccs)))
            print(f"  cost={cost:12s}  PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}")
    elif args.exp == "hreg":
        for h in [1e-3, 3e-3, 1e-2, 3e-2, 1e-1]:
            pccs = []
            for seed in range(args.seeds):
                pcc, *_ = score_native(name, adata_path, n_hvf, w, seed,
                                       optim_name="adam", lr=1e-2, max_iter=50,
                                       max_iter_inner=150, h=h)
                pccs.append(pcc)
            rows.append(dict(h=h, PCC_mean=np.mean(pccs), PCC_std=np.std(pccs)))
            print(f"  h={h:<6g}  PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}")
    elif args.exp == "lropt":
        for optim, lr in [("adam", 3e-3), ("adam", 1e-2), ("adam", 3e-2),
                          ("adamw", 1e-2), ("nadam", 1e-2), ("rmsprop", 1e-2),
                          ("rmsprop", 3e-3)]:
            pccs, times = [], []
            for seed in range(args.seeds):
                pcc, rmse, jsd, dt, *_ = score_native(
                    name, adata_path, n_hvf, w, seed, optim_name=optim, lr=lr,
                    max_iter=50, max_iter_inner=150)
                pccs.append(pcc); times.append(dt)
            rows.append(dict(optim=optim, lr=lr, PCC_mean=np.mean(pccs),
                             PCC_std=np.std(pccs), sec=np.mean(times)))
            print(f"  {optim:7s} lr={lr:<6g}  PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}  {np.mean(times):.0f}s")
    elif args.exp == "persist":
        for persist in [False, True]:
            for inner in [60, 150]:
                pccs, times = [], []
                for seed in range(args.seeds):
                    pcc, rmse, jsd, dt, *_ = score_native(
                        name, adata_path, n_hvf, w, seed, optim_name="adam", lr=1e-2,
                        max_iter=50, max_iter_inner=inner, persist_optimizer=persist)
                    pccs.append(pcc); times.append(dt)
                rows.append(dict(persist=persist, inner=inner, PCC_mean=np.mean(pccs),
                                 PCC_std=np.std(pccs), sec=np.mean(times)))
                print(f"  persist={str(persist):5s} inner={inner:4d}  "
                      f"PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}  {np.mean(times):.0f}s")
    elif args.exp == "anneal":
        # w annealing: start gentler (w*factor) and decay to target w. Also test
        # whether annealing lets a higher inner count (400) avoid the collapse.
        for inner in [150, 400]:
            for wa in [1.0, 2.0, 5.0, 10.0]:
                pccs = []
                for seed in range(args.seeds):
                    pcc, *_ = score_native(name, adata_path, n_hvf, w, seed,
                                           optim_name="adam", lr=1e-2, max_iter=50,
                                           max_iter_inner=inner, w_anneal=wa)
                    pccs.append(pcc)
                rows.append(dict(inner=inner, w_anneal=wa, PCC_mean=np.mean(pccs),
                                 PCC_std=np.std(pccs)))
                print(f"  inner={inner:4d} w_anneal={wa:<5g}  "
                      f"PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}")
    elif args.exp == "eps":
        # OT-geometry lever: entropic smoothing eps (default 0.02).
        for eps in [0.005, 0.01, 0.02, 0.05, 0.1]:
            pccs, times = [], []
            for seed in range(args.seeds):
                pcc, rmse, jsd, dt, nl, last, _ = score_native(
                    name, adata_path, n_hvf, w, seed, optim_name="adam", lr=1e-2,
                    max_iter=50, max_iter_inner=150, eps_override=eps)
                pccs.append(pcc); times.append(dt)
            rows.append(dict(eps=eps, PCC_mean=np.mean(pccs), PCC_std=np.std(pccs),
                             sec=np.mean(times)))
            print(f"  eps={eps:<6g}  PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}  {np.mean(times):.0f}s")
    elif args.exp == "damping":
        # Does damping stabilize the outer alternation -- beating 0.573 and/or
        # rescuing the over-converged (inner=400) collapse?
        for inner in [150, 400]:
            for damp in [1.0, 0.7, 0.5, 0.3]:
                pccs, times = [], []
                for seed in range(args.seeds):
                    pcc, rmse, jsd, dt, nl, last, _ = score_native(
                        name, adata_path, n_hvf, w, seed, optim_name="adam", lr=1e-2,
                        max_iter=50, max_iter_inner=inner, damping=damp)
                    pccs.append(pcc); times.append(dt)
                rows.append(dict(inner=inner, damping=damp, PCC_mean=np.mean(pccs),
                                 PCC_std=np.std(pccs), sec=np.mean(times)))
                print(f"  inner={inner:4d} damp={damp:<4g}  "
                      f"PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}  {np.mean(times):.0f}s")
    elif args.exp == "init":
        for init in ["random", "nndsvd"]:
            pccs, times = [], []
            for seed in range(args.seeds):
                pcc, rmse, jsd, dt, nl, last, _ = score_native(
                    name, adata_path, n_hvf, w, seed, optim_name="adam", lr=1e-2,
                    max_iter=50, max_iter_inner=150, init=init)
                pccs.append(pcc); times.append(dt)
            rows.append(dict(init=init, PCC_mean=np.mean(pccs), PCC_std=np.std(pccs),
                             PCC_min=np.min(pccs), PCC_max=np.max(pccs), sec=np.mean(times)))
            print(f"  init={init:8s}  PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f} "
                  f"[{np.min(pccs):.4f},{np.max(pccs):.4f}]  {np.mean(times):.0f}s")
    elif args.exp == "batch":
        # Stochastic spot mini-batching in the inner dual ascent (float32 so the
        # analytic path -- which supports batching -- is active). Hypothesis: SGD
        # noise regularizes like the inexact inner loop and finds a better optimum.
        for bs in [None, 32, 64, 128]:
            pccs, times = [], []
            for seed in range(args.seeds):
                pcc, rmse, jsd, dt, *_ = score_native(
                    name, adata_path, n_hvf, w, seed, optim_name="adam", lr=1e-2,
                    max_iter=50, max_iter_inner=150, dtype="float32", batch_size=bs)
                pccs.append(pcc); times.append(dt)
            rows.append(dict(batch_size=bs, PCC_mean=np.mean(pccs), PCC_std=np.std(pccs),
                             sec=np.mean(times)))
            print(f"  batch_size={str(bs):5s}  PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}"
                  f"  {np.mean(times):.1f}s")
    elif args.exp == "gtail":
        # Polyak tail-averaging of the dual G before the peaked-softmin block
        # update. Hypothesis: averaging the last g_tail near-optimum iterates
        # cancels the zero-mean dual jitter the softmin amplifies -> higher PCC
        # and/or lower seed variance (a within-run analogue of consensus).
        for gt in [0, 5, 10, 20, 40]:
            pccs, times = [], []
            for seed in range(args.seeds):
                pcc, rmse, jsd, dt, *_ = score_native(
                    name, adata_path, n_hvf, w, seed, optim_name="adam", lr=1e-2,
                    max_iter=50, max_iter_inner=150, g_tail=gt)
                pccs.append(pcc); times.append(dt)
            rows.append(dict(g_tail=gt, PCC_mean=np.mean(pccs), PCC_std=np.std(pccs),
                             sec=np.mean(times)))
            print(f"  g_tail={gt:4d}  PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}"
                  f"  {np.mean(times):.0f}s")
    elif args.exp == "transform":
        # Dynamic-range-compressing input transforms applied BEFORE the model's
        # own row/column normalization. Hypothesis: compressing the range makes
        # the OT solve better-conditioned (less noise-sensitive) -> higher PCC.
        for tf in ["none", "log1p", "sqrt", "cbrt", "arcsinh"]:
            pccs = []
            for seed in range(args.seeds):
                pcc, *_ = score_native(name, adata_path, n_hvf, w, seed,
                                       optim_name="adam", lr=1e-2, max_iter=50,
                                       max_iter_inner=150, x_transform=tf)
                pccs.append(pcc)
            rows.append(dict(transform=tf, PCC_mean=np.mean(pccs), PCC_std=np.std(pccs)))
            print(f"  transform={tf:8s}  PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}")
    elif args.exp == "rowscale":
        # Partial gene-mass equalization: divide rows by rowsum**alpha before the
        # column normalization. alpha=1 is the current default, 0 is no row norm.
        for a in [0.0, 0.25, 0.5, 0.75, 1.0]:
            pccs = []
            for seed in range(args.seeds):
                pcc, *_ = score_native(name, adata_path, n_hvf, w, seed,
                                       optim_name="adam", lr=1e-2, max_iter=50,
                                       max_iter_inner=150, normalize_rows=(a > 0),
                                       row_power=a)
                pccs.append(pcc)
            rows.append(dict(row_power=a, PCC_mean=np.mean(pccs), PCC_std=np.std(pccs)))
            print(f"  row_power={a:<4g}  PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}")
    elif args.exp == "agrad":
        # Closed-form dual gradient vs autograd backward. Same gradient to
        # ~1e-14, so PCC must match within seed noise; the win is wall-time.
        for dt_name in ["float64", "float32"]:
            for ag in [False, True]:
                pccs, times = [], []
                for seed in range(args.seeds):
                    pcc, rmse, jsd, dt, nl, last, _ = score_native(
                        name, adata_path, n_hvf, w, seed, optim_name="adam", lr=1e-2,
                        max_iter=50, max_iter_inner=150, dtype=dt_name, analytic_grad=ag)
                    pccs.append(pcc); times.append(dt)
                rows.append(dict(dtype=dt_name, analytic_grad=ag, PCC_mean=np.mean(pccs),
                                 PCC_std=np.std(pccs), sec=np.mean(times)))
                print(f"  {dt_name:8s} analytic_grad={str(ag):5s}  "
                      f"PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}  {np.mean(times):.1f}s")
    elif args.exp == "dtype":
        # Precision lever: float64 (current default) vs float32. Measures the
        # accuracy/speed/memory trade. The OT dual is log-sum-exp stabilized, so
        # float32 should hold PCC while halving G/A/K + Adam buffers and running
        # much faster on GPUs with weak FP64 throughput (e.g. RTX 3080 ~1:32).
        for dt_name in ["float64", "float32"]:
            pccs, rmses, times, mems = [], [], [], []
            for seed in range(args.seeds):
                pcc, rmse, jsd, dt, nl, last, pk = score_native(
                    name, adata_path, n_hvf, w, seed, optim_name="adam", lr=1e-2,
                    max_iter=50, max_iter_inner=150, dtype=dt_name, track_mem=True)
                pccs.append(pcc); rmses.append(rmse); times.append(dt); mems.append(pk)
            rows.append(dict(dtype=dt_name, PCC_mean=np.mean(pccs), PCC_std=np.std(pccs),
                             RMSE=np.mean(rmses), sec=np.mean(times),
                             peak_mb=np.nanmean(mems)))
            print(f"  dtype={dt_name:8s}  PCC={np.mean(pccs):.4f}±{np.std(pccs):.4f}  "
                  f"RMSE={np.mean(rmses):.4f}  {np.mean(times):.1f}s  "
                  f"peak={np.nanmean(mems):.0f}MB")
    elif args.exp == "seeds":
        for seed in range(args.seeds):
            pcc, rmse, jsd, dt, nl, last, _ = score_native(
                name, adata_path, n_hvf, w, seed, optim_name="adam", lr=1e-2,
                max_iter=50, max_iter_inner=150)
            rows.append(dict(seed=seed, PCC_mean=pcc, RMSE=rmse, JSD=jsd, sec=dt))
            print(f"  seed={seed}  PCC={pcc:.4f} RMSE={rmse:.4f}  {dt:.0f}s")

    df = pd.DataFrame(rows)
    out = args.out or str(REPO / "results" / f"base_algo_{args.exp}_{name[:8]}.csv")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
