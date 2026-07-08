"""Ablation: how the consensus clustering *method* and *aggregation* affect
spOT-NMF's consensus quality.

spOT-NMF's consensus was originally a direct port of cNMF (L2-normalize -> KMeans
-> per-cluster **median**). But spOT-NMF spectra are probability distributions
over genes, so cNMF's Euclidean/median choices are geometrically mismatched. This
script quantifies the effect of two knobs added in ``spotnmf.models``:

  * ``method``    -- ``"kmeans"`` (cNMF-style, with local-density outlier filter)
                     vs ``"match"`` (balanced cross-replicate Hungarian matching)
  * ``aggregate`` -- ``"median"`` (cNMF) vs ``"mean"`` (mixture barycenter) vs
                     ``"geomean"`` (KL barycenter -- the OT/info-geometry native)

Protocol: factorize ``n_iter`` replicates per dataset ONCE (cached), then for
every (method, aggregate) pair build the consensus spectra with the *package*
:func:`spotnmf.models._cluster_consensus_spectra` and score them with a *common
NNLS usage refit* against ground truth. Using one fixed refit isolates the effect
of the consensus algorithm on the SPECTRA (independent of OT-vs-NNLS usage
choices). Also reports cross-replicate program agreement (reproducibility) and
the single-run (no-consensus) baseline.

Writes ``results/consensus_ablation.csv``.

Usage::

    python scripts/benchmark/consensus_ablation.py
    python scripts/benchmark/consensus_ablation.py --datasets MOB_dance_sim_norm_mm
"""
import argparse
import contextlib
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.decomposition import non_negative_factorization

import mosaicmpi
import spotnmf.mosaic  # noqa: F401  (registers the "spotnmf" backend)
from benchmark_mosaic import _row_normalize, match_and_score, REPO
from spotnmf.models import run_spotnmf, _cluster_consensus_spectra

import torch

BASE = Path(
    r"Z:\University of Calgary\Morrissy Lab - Documents\Visium_profiling"
    r"\spotnmf_benchmark\benchmark_data\data"
)
CACHE = REPO / "results" / "_consensus_diag"          # reuse the diagnostic pools
STEREO_CACHE = REPO / "results" / "benchmark_all" / "_cache" / \
    "stereoseq_mouse_brain_li2023_sim_norm_mm_adata_spatial.h5ad"
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# (name, adata_path override, n_hvf, tuned w)  -- w matches benchmark_all OVERRIDES
DATASETS = [
    ("Dataset10_STARmap_li2022_sim_norm_mm", None, None, 0.01),
    ("Dataset4_seqFISH_li2022_sim_norm_mm", None, 2000, 0.1),
    ("MOB_dance_sim_norm_mm", None, 2000, 0.05),
    ("stereoseq_mouse_brain_li2023_sim_norm_mm", str(STEREO_CACHE), 2000, 0.01),
]
METHODS = [("kmeans", "median"), ("kmeans", "mean"), ("kmeans", "geomean"),
           ("match", "median"), ("match", "mean"), ("match", "geomean")]


def get_norm(name, adata_path, n_hvf):
    """mosaicMPI-normalized (HVF-selected) matrix + ground-truth proportions."""
    path = adata_path or str(BASE / name / "adata_spatial.h5ad")
    adata = ad.read_h5ad(path)
    adata.obs_names = ["spot" + str(i) for i in range(adata.n_obs)]
    adata.var_names = adata.var_names.astype(str)
    gt = adata.uns["ground_truth"]; gt.index = adata.obs_names
    gt_prop = _row_normalize(gt.astype(float)); k = gt.shape[1]
    ds = mosaicmpi.Dataset(adata.copy())
    with open(os.devnull, "w") as nul, contextlib.redirect_stdout(nul), \
            contextlib.redirect_stderr(nul):
        if n_hvf is None:
            ds.select_hvf(feature_list=list(adata.var_names), score_type="odscore")
        else:
            ds.select_hvf(score_type="odscore", top_n=n_hvf)
        run = ds.initialize_cnmf(cnmf_output_dir=str(CACHE / f"mm_{name[:8]}"),
                                 cnmf_name="d", kvals=[k], n_iter=1, algorithm="spotnmf")
    norm = ad.read_h5ad(run.paths["normalized_counts"])
    X = norm.X.toarray() if hasattr(norm.X, "toarray") else np.asarray(norm.X)
    return np.asarray(X, np.float64), norm.var.index, norm.obs.index, gt_prop, k


def replicate_pool(name, adata_path, n_hvf, w, n_iter):
    """Return ``(pool (n_iter,k,genes), X, obs_idx, var_idx, gt_prop, k)``, cached."""
    X, var_idx, obs_idx, gt_prop, k = get_norm(name, adata_path, n_hvf)
    cache = CACHE / f"pool_{name[:10]}.npz"
    if cache.exists():
        pool = np.load(cache, allow_pickle=True)["pool"]
        return pool, X, obs_idx, var_idx, gt_prop, k
    CACHE.mkdir(parents=True, exist_ok=True)
    sub0 = ad.AnnData(X=X, obs=pd.DataFrame(index=obs_idx),
                      var=pd.DataFrame(index=var_idx))
    sub0.var["highly_variable"] = True
    max_iter, inner = (40, 80) if name.startswith("stereoseq") else (50, 150)
    progs = []
    for i in range(n_iter):
        sub = sub0.copy()
        with open(os.devnull, "w") as nul, contextlib.redirect_stdout(nul), \
                contextlib.redirect_stderr(nul):
            res, _ = run_spotnmf(sub, components=k, seed=i, device=DEVICE,
                                 optim_name="adam", lr=0.01, max_iter=max_iter,
                                 max_iter_inner=inner, eps=0.02, h=1e-2, w=w,
                                 normalize_rows=True, cost="cosine")
        progs.append(res["genes_per_topic"].T.to_numpy())
        print(f"  [{name[:10]}] replicate {i + 1}/{n_iter}")
    pool = np.stack(progs)
    np.savez_compressed(cache, pool=pool)
    return pool, X, obs_idx, var_idx, gt_prop, k


def nnls_score(H, X, obs_idx, gt_prop, k):
    W, _, _ = non_negative_factorization(
        X, H=np.ascontiguousarray(H), n_components=k, init="custom",
        update_H=False, solver="cd", max_iter=500, random_state=0)
    usage = pd.DataFrame(W, index=obs_idx, columns=[str(i + 1) for i in range(k)])
    s = match_and_score(usage, gt_prop)
    return s["PCC_mean"], s["RMSE"], s["JSD"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="*", default=None)
    ap.add_argument("--n_iter", type=int, default=12)
    ap.add_argument("--out", default=str(REPO / "results" / "consensus_ablation.csv"))
    args = ap.parse_args()

    selected = DATASETS
    if args.datasets:
        want = set(args.datasets)
        selected = [d for d in DATASETS if d[0] in want]

    rows = []
    for name, adata_path, n_hvf, w in selected:
        if adata_path and not Path(adata_path).exists():
            print(f"SKIP {name}: {adata_path} not found (build it via benchmark_all first)")
            continue
        print(f"\n### {name}  (w={w})")
        pool, X, obs_idx, var_idx, gt_prop, k = replicate_pool(
            name, adata_path, n_hvf, w, args.n_iter)
        n_iter = pool.shape[0]
        spectra_pool = pd.DataFrame(
            pool.reshape(n_iter * k, -1), columns=var_idx,
            index=[f"iter{i}_topic{t + 1}" for i in range(n_iter) for t in range(k)])

        pcc0, rmse0, jsd0 = nnls_score(pool[0] / pool[0].sum(1, keepdims=True),
                                       X, obs_idx, gt_prop, k)
        rows.append(dict(dataset=name, k=k, method="single-run", aggregate="-",
                         PCC=pcc0, RMSE=rmse0, JSD=jsd0, agreement=np.nan))
        print(f"  single-run              PCC={pcc0:.3f} RMSE={rmse0:.3f} JSD={jsd0:.3f}")
        for method, agg in METHODS:
            H, agreement = _cluster_consensus_spectra(
                spectra_pool, k, 0.5, 0.30, n_iter=n_iter, method=method, aggregate=agg)
            pcc, rmse, jsd = nnls_score(H.to_numpy(), X, obs_idx, gt_prop, k)
            rows.append(dict(dataset=name, k=k, method=method, aggregate=agg,
                             PCC=pcc, RMSE=rmse, JSD=jsd, agreement=agreement))
            print(f"  {method:6s}+{agg:8s}      PCC={pcc:.3f} RMSE={rmse:.3f} "
                  f"JSD={jsd:.3f}  agreement={agreement:.3f}")

    df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
