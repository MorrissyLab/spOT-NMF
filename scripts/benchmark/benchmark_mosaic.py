"""Benchmark NMF vs. cNMF vs. spOT-NMF as factorization backends inside mosaicMPI.

The mosaicMPI paper (Verhey & Morrissy, *Nucleic Acids Research* 2024,
doi:10.1093/nar/gkae421) is an **integration** framework paper: it fixes cNMF as
its deconvolution engine and does not benchmark factorization algorithms for
accuracy against alternatives. So there is no ready-made benchmark to reuse; this
script builds one on the simulated STARmap slide that ships with spOT-NMF
(``data/test_data/dataset10_adata_spatial.h5ad``), which carries a ground-truth
cell-type-count matrix per spot (``adata.uns["ground_truth"]``).

All three algorithms run through the *identical* mosaicMPI pipeline (HVF
selection -> per-replicate factorization -> consensus -> refit). Only the
per-replicate factorizer changes:

  * ``NMF``      -- scikit-learn NMF, single replicate (n_iter=1, no consensus)
  * ``cNMF``     -- scikit-learn NMF, consensus over ``n_iter`` replicates
  * ``spOT-NMF`` -- optimal-transport NMF, consensus over ``n_iter`` replicates

We factorize at rank k = number of ground-truth cell types, match each learned
program to a cell type (Hungarian assignment on across-spot correlation), and
score the recovered per-spot usage against the ground-truth proportions with
metrics standard in spatial-deconvolution benchmarking:

  * PCC   -- mean across-spot Pearson correlation per cell type (higher better)
  * SCC   -- same, Spearman (higher better)
  * RMSE  -- root mean squared error of per-spot proportions (lower better)
  * JSD   -- mean per-spot Jensen-Shannon divergence (lower better)

mosaicMPI's own rank-selection stats (stability, reconstruction error) are also
reported for context.

Usage::

    python scripts/benchmark/benchmark_mosaic.py            # full run
    python scripts/benchmark/benchmark_mosaic.py --quick    # fast smoke test
"""

import argparse
import contextlib
import os
import shutil
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr, spearmanr

from sklearn.decomposition import non_negative_factorization

import mosaicmpi
import spotnmf.mosaic  # noqa: F401  (registers the "spotnmf" backend on import)


REPO = Path(__file__).resolve().parents[2]
DEFAULT_DATA = REPO / "data" / "test_data" / "dataset10_adata_spatial.h5ad"


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #
def _row_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Convert a non-negative matrix to per-row proportions (rows summing to 1)."""
    s = df.sum(axis=1).replace(0, np.nan)
    return df.div(s, axis=0).fillna(0.0)


def match_and_score(usage: pd.DataFrame, gt_prop: pd.DataFrame) -> dict:
    """Match learned programs to ground-truth cell types and score the fit.

    :param usage: spots x programs, non-negative (learned per-spot usage).
    :param gt_prop: spots x cell-types, per-spot proportions (ground truth).
    :return: dict of aggregate metrics plus the per-cell-type correlation table.
    """
    # Align spots.
    usage = usage.loc[gt_prop.index]
    pred_prop = _row_normalize(usage)

    P = pred_prop.to_numpy()            # spots x programs
    G = gt_prop.to_numpy()              # spots x celltypes
    n_prog, n_ct = P.shape[1], G.shape[1]

    # Cost = negative Pearson correlation between each (program, celltype) pair
    # across spots. Hungarian assignment maximises total correlation.
    corr = np.zeros((n_prog, n_ct))
    for i in range(n_prog):
        for j in range(n_ct):
            if P[:, i].std() == 0 or G[:, j].std() == 0:
                corr[i, j] = 0.0
            else:
                corr[i, j] = pearsonr(P[:, i], G[:, j])[0]
    corr = np.nan_to_num(corr)
    row_ind, col_ind = linear_sum_assignment(-corr)  # program i -> celltype col_ind

    # Per-cell-type scores on the matched pairs.
    pccs, sccs = [], []
    matched = {}
    for i, j in zip(row_ind, col_ind):
        ct = gt_prop.columns[j]
        prog = pred_prop.columns[i]
        matched[ct] = prog
        p = P[:, i]
        g = G[:, j]
        pcc = 0.0 if p.std() == 0 or g.std() == 0 else pearsonr(p, g)[0]
        scc = 0.0 if p.std() == 0 or g.std() == 0 else spearmanr(p, g)[0]
        pccs.append((ct, pcc))
        sccs.append((ct, scc))

    per_ct = pd.DataFrame({
        "cell_type": [c for c, _ in pccs],
        "matched_program": [matched[c] for c, _ in pccs],
        "PCC": [v for _, v in pccs],
        "SCC": [v for _, v in sccs],
    }).set_index("cell_type")

    # Reorder predicted proportions to the matched cell-type order for RMSE/JSD.
    pred_aligned = pred_prop[[matched[ct] for ct in gt_prop.columns]].to_numpy()
    pred_aligned = np.nan_to_num(pred_aligned)
    rmse = float(np.sqrt(np.mean((pred_aligned - G) ** 2)))
    jsds = []
    for s in range(G.shape[0]):
        if G[s].sum() > 0 and pred_aligned[s].sum() > 0:
            d = jensenshannon(pred_aligned[s], G[s], base=2)
            if np.isfinite(d):
                jsds.append(d ** 2)
    jsd = float(np.mean(jsds)) if jsds else float("nan")

    # Mean per-spot spread of the predicted mixture. ~0 means every spot gets an
    # (near-)identical program mixture, i.e. no spatial deconvolution signal.
    dispersion = float(pred_prop.std(axis=1).mean())

    return {
        "PCC_mean": float(per_ct["PCC"].mean()),
        "PCC_median": float(per_ct["PCC"].median()),
        "SCC_mean": float(per_ct["SCC"].mean()),
        "RMSE": rmse,
        "JSD": jsd,
        "usage_dispersion": dispersion,
        "per_celltype": per_ct,
    }


# --------------------------------------------------------------------------- #
# Running one algorithm through mosaicMPI
# --------------------------------------------------------------------------- #
def run_algorithm(adata, name, algorithm, n_iter, k, out_root,
                  factorizer_params=None, consensus=True):
    """Run one backend through mosaicMPI and return usages + rank stats.

    When ``consensus`` is False (the plain-NMF baseline) we run a *single*
    scikit-learn factorization on the exact ``norm_counts`` matrix mosaicMPI would
    feed its consensus replicates, using mosaicMPI's own KL solver settings. This
    isolates the effect of consensus: same input, same solver, no replicate
    clustering / refit.
    """
    run_name = f"bench_{name}"
    out_dir = str(out_root / name)
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)

    ds = mosaicmpi.Dataset(adata.copy())
    # Use the full targeted gene panel as features so all backends see the same
    # inputs (the STARmap panel is small and already curated).
    ds.select_hvf(feature_list=list(adata.var_names), score_type="odscore")

    t0 = time.time()
    run = ds.initialize_cnmf(
        cnmf_output_dir=out_dir,
        cnmf_name=run_name,
        kvals=[k],
        n_iter=n_iter,
        algorithm=algorithm,
        factorizer_params=factorizer_params,
    )

    if not consensus:
        # Single factorization on mosaicMPI's normalized input matrix, using each
        # method's *own* usage matrix (no consensus, no NNLS refit). This isolates
        # raw algorithm quality: for spOT-NMF it keeps the OT-computed usages that
        # the mosaicMPI consensus refit would otherwise replace.
        norm = ad.read_h5ad(run.paths["normalized_counts"])
        X = norm.X.toarray() if hasattr(norm.X, "toarray") else np.asarray(norm.X)
        if algorithm == "cnmf":
            W, _H, _n = non_negative_factorization(
                X, n_components=k, init="random", solver="mu",
                beta_loss="kullback-leibler", max_iter=1000, random_state=0)
        else:
            factorizer = mosaicmpi.factorization.get_factorizer(algorithm)
            with open(os.devnull, "w") as _null, \
                    contextlib.redirect_stdout(_null), contextlib.redirect_stderr(_null):
                _spectra, W = factorizer(X, k, 0, var_names=list(norm.var.index),
                                         params=factorizer_params)
        usage = pd.DataFrame(np.asarray(W), index=norm.obs.index,
                             columns=[str(i + 1) for i in range(k)])
        return {"usage": usage, "seconds": time.time() - t0,
                "stability": np.nan, "recon_err": np.nan}

    with open(os.devnull, "w") as _null, \
            contextlib.redirect_stdout(_null), contextlib.redirect_stderr(_null):
        run.factorize(verbose=False)
        run.postprocess()
    ds.add_cnmf_results(cnmf_output_dir=out_dir, cnmf_name=run_name)
    elapsed = time.time() - t0

    usage = ds.get_usages(k=k)
    if isinstance(usage.columns, pd.MultiIndex):
        usage = usage.xs(k, axis=1, level=0) if k in usage.columns.get_level_values(0) else usage
        usage.columns = [str(c[-1]) if isinstance(c, tuple) else str(c) for c in usage.columns]
    usage.columns = [str(c) for c in usage.columns]

    # mosaicMPI rank-selection stats for context.
    kvals = ds.adata.uns.get("kvals")
    stability = recon_err = np.nan
    if isinstance(kvals, pd.DataFrame) and k in kvals.index:
        for col in kvals.columns:
            lc = col.lower()
            if "stability" in lc or "silhouette" in lc:
                stability = float(kvals.loc[k, col])
            if "error" in lc or "recon" in lc:
                recon_err = float(kvals.loc[k, col])

    return {"usage": usage, "seconds": elapsed,
            "stability": stability, "recon_err": recon_err}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=str(DEFAULT_DATA))
    ap.add_argument("--out", default=str(REPO / "results" / "benchmark_mosaic"))
    ap.add_argument("--quick", action="store_true",
                    help="fast smoke test (few replicates / iterations)")
    args = ap.parse_args()

    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.data)
    # mosaicMPI round-trips usages through TSV files; a purely numeric obs index
    # gets reparsed as int and then fails to align, so use non-numeric spot names.
    adata.obs_names = ["spot" + str(i) for i in range(adata.n_obs)]
    adata.var_names = adata.var_names.astype(str)
    gt = adata.uns["ground_truth"]
    gt.index = adata.obs_names[: len(gt)] if len(gt) == adata.n_obs else gt.index
    gt.index = adata.obs_names
    gt_prop = _row_normalize(gt.astype(float))
    k = gt.shape[1]
    print(f"Dataset: {adata.uns.get('dataset_name')}  spots={adata.n_obs} "
          f"genes={adata.n_vars}  ground-truth cell types (k)={k}")

    # n_iter must be >= 4 or mosaicMPI's consensus neighbourhood collapses to 0.
    n_iter = 6 if args.quick else 12
    # IMPORTANT: use the adam optimizer. spOT-NMF's LBFGS path leaves the OT dual
    # variable G ~unchanged, so W = softmin(coef * H.T @ G) collapses to a uniform
    # per-spot mixture (usage_dispersion == 0) and carries no deconvolution signal
    # -- regardless of lr. adam moves G and recovers spot-varying usages.
    import torch
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"spOT-NMF device: {device}")
    # Tuned via scripts/benchmark/sweep_spotnmf.py (72-config GPU sweep):
    # eps=0.02, h=0.01, w=0.01, normalize_rows=True, cosine cost was best.
    spot_params = {"device": device, "optim_name": "adam", "lr": 1e-2,
                   "max_iter": 20 if args.quick else 50,
                   "max_iter_inner": 60 if args.quick else 150,
                   "eps": 0.02, "h": 1e-2, "w": 1e-2,
                   "normalize_rows": True, "cost": "cosine"}

    configs = [
        # (label, backend, n_iter, factorizer_params, consensus)
        ("NMF",             "cnmf",    1,      None,        False),  # single scikit-learn NMF (own usages)
        ("spOT-NMF-native", "spotnmf", 1,      spot_params, False),  # single OT-NMF (own OT usages)
        ("cNMF",            "cnmf",    n_iter, None,        True),   # consensus scikit-learn NMF
        ("spOT-NMF",        "spotnmf", n_iter, spot_params, True),   # consensus OT-NMF (mosaicMPI refit)
    ]

    rows = []
    per_ct_all = {}
    for label, backend, niter, fp, consensus in configs:
        print(f"\n=== {label}  (backend={backend}, n_iter={niter}, consensus={consensus}) ===")
        res = run_algorithm(adata, label, backend, niter, k, out_root, fp, consensus)
        scores = match_and_score(res["usage"], gt_prop)
        per_ct_all[label] = scores.pop("per_celltype")
        rows.append({
            "algorithm": label, "backend": backend, "n_iter": niter,
            "PCC_mean": scores["PCC_mean"], "PCC_median": scores["PCC_median"],
            "SCC_mean": scores["SCC_mean"], "RMSE": scores["RMSE"],
            "JSD": scores["JSD"], "usage_dispersion": scores["usage_dispersion"],
            "stability": res["stability"],
            "recon_err": res["recon_err"], "seconds": res["seconds"],
        })
        print(f"    PCC_mean={scores['PCC_mean']:.3f}  SCC_mean={scores['SCC_mean']:.3f}  "
              f"RMSE={scores['RMSE']:.3f}  JSD={scores['JSD']:.3f}  "
              f"({res['seconds']:.1f}s)")

    summary = pd.DataFrame(rows).set_index("algorithm")
    summary.to_csv(out_root / "summary.csv")
    for label, df in per_ct_all.items():
        df.to_csv(out_root / f"per_celltype_{label}.csv")

    print("\n" + "=" * 70)
    print("BENCHMARK SUMMARY  (higher PCC/SCC better; lower RMSE/JSD better)")
    print("=" * 70)
    with pd.option_context("display.width", 200, "display.float_format", "{:.3f}".format):
        print(summary[["PCC_mean", "PCC_median", "SCC_mean", "RMSE", "JSD",
                       "usage_dispersion", "stability", "recon_err", "seconds"]])

    winner = summary["PCC_mean"].idxmax()
    print(f"\nBest per-spot deconvolution accuracy (PCC_mean): {winner}")

    try:
        make_plots(summary, per_ct_all, out_root)
        print(f"Plots written to: {out_root / 'benchmark.png'}")
    except Exception as e:  # plotting is best-effort
        print(f"(plotting skipped: {e})")
    print(f"Results written to: {out_root}")


def make_plots(summary, per_ct_all, out_root):
    """Bar charts of the aggregate metrics and per-cell-type PCC distributions."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    algos = list(summary.index)
    colors = {"NMF": "#9e9e9e", "spOT-NMF-native": "#f58518",
              "cNMF": "#4c78a8", "spOT-NMF": "#e45756"}
    palette = [colors.get(a, "#59a14f") for a in algos]

    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    for ax, metric, better in [
        (axes[0], "PCC_mean", "higher better"),
        (axes[1], "SCC_mean", "higher better"),
        (axes[2], "RMSE", "lower better"),
        (axes[3], "JSD", "lower better"),
    ]:
        ax.bar(algos, summary[metric].values, color=palette)
        ax.set_title(f"{metric}\n({better})")
        ax.set_ylabel(metric)
        for i, v in enumerate(summary[metric].values):
            ax.text(i, v, f"{v:.3f}", ha="center", va="bottom", fontsize=9)
        ax.tick_params(axis="x", rotation=15)
    fig.suptitle("mosaicMPI factorization backends — deconvolution accuracy "
                 "(simulated STARmap, ground-truth cell types)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_root / "benchmark.png", dpi=150)
    plt.close(fig)

    # Per-cell-type PCC boxplot.
    fig, ax = plt.subplots(figsize=(6, 4))
    data = [per_ct_all[a]["PCC"].values for a in algos]
    bp = ax.boxplot(data, labels=algos, patch_artist=True, showmeans=True)
    for patch, c in zip(bp["boxes"], palette):
        patch.set_facecolor(c)
        patch.set_alpha(0.6)
    ax.set_ylabel("per-cell-type PCC (across spots)")
    ax.set_title("Per-cell-type recovery")
    ax.tick_params(axis="x", rotation=15)
    fig.tight_layout()
    fig.savefig(out_root / "benchmark_per_celltype.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
