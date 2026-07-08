"""Run the mosaicMPI-backend deconvolution benchmark across all manuscript datasets.

This drives :func:`benchmark_mosaic.run_benchmark` over the simulated
spatial-transcriptomics datasets used in the spOT-NMF manuscript's *normal*
benchmark (``scripts/manuscript_notebooks/spotnmf3_normal_benchmark.ipynb`` on
the ``manuscript`` branch). For every dataset it compares the five factorization
configurations (single NMF, single OT-NMF, consensus cNMF, consensus OT + NNLS
refit, and consensus OT + OT usage refit) against the ground-truth cell-type
proportions, writes a per-dataset ``summary.csv`` + plots, and finally a combined
``summary_all.csv`` / per-metric win table across datasets.

The datasets live in a shared data root (Dropbox). Each expected as
``<root>/<dataset>/adata_spatial.h5ad`` carrying ``uns['ground_truth']``. Small,
curated panels use the full gene set; transcriptome-wide datasets select the
top-N overdispersed genes (``--n_hvf``) to stay tractable, matching the
manuscript's overdispersed-gene selection.

Usage::

    python scripts/benchmark/benchmark_all.py                 # full run, default root
    python scripts/benchmark/benchmark_all.py --quick         # fast smoke test
    python scripts/benchmark/benchmark_all.py --data_root PATH --n_hvf 2000
"""

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd

from benchmark_mosaic import run_benchmark, REPO


# Default shared data root (Morrissy Lab Dropbox). Override with --data_root.
DEFAULT_ROOT = Path(
    r"Z:\University of Calgary\Morrissy Lab - Documents\Visium_profiling"
    r"\spotnmf_benchmark\benchmark_data\data"
)

# Spots to subsample from the (large) stereo-seq slide. A single spOT-NMF
# replicate over 3000 spots x 2000 HVGs x k=29 takes ~15 min on an RTX 3080, so a
# full 12-replicate consensus is impractical; we subsample to 1000 spots (fixed
# seed) and use lighter iteration caps (see OVERRIDES) to keep the run tractable
# while preserving the k=29 multi-cell-per-spot deconvolution structure.
STEREO_N_SPOTS = 1000
STEREO_SEED = 0

# Datasets from the manuscript normal benchmark.
# ``n_hvf``: None keeps the full curated panel; an int selects top-N
# overdispersed genes for transcriptome-wide datasets. Datasets whose name is in
# LOADERS (below) are built on the fly from raw files instead of a shipped
# ``adata_spatial.h5ad``.
DATASETS = [
    ("Dataset10_STARmap_li2022_sim_norm_mm", None),   # 189 x 882 curated panel, k=15
    ("Dataset4_seqFISH_li2022_sim_norm_mm", 2000),    # 72 x 9684, k=13
    ("MOB_dance_sim_norm_mm", 2000),                  # 260 x 18263, k=6
    ("stereoseq_mouse_brain_li2023_sim_norm_mm", 2000),  # 8793 x 25879, k=29 (subsampled)
]
# Excluded:
#   * Synthetic_Spotlight -- adata_spatial.h5ad has 0 genes in X (counts stored
#                            separately in mix_count.csv)


def build_stereoseq_adata(dataset_dir, cache_path,
                          n_spots=STEREO_N_SPOTS, seed=STEREO_SEED):
    """Build a spOT-NMF-ready AnnData for the stereo-seq mouse-brain slide.

    The slide ships only raw simulation outputs (no ``adata_spatial.h5ad`` and no
    in-repo loader -- the manuscript used an external ``get_xeno_spatial_data``).
    We reconstruct it from ``simulated_data/``:

      * ``combined_spatial_count.txt``  -> expression matrix (spots x genes)
      * ``combined_spot_clusters.txt``  -> per-spot cell-type counts (ground truth)

    Spots are subsampled to ``n_spots`` (fixed seed) for tractability; all-zero
    genes and cell-type columns are dropped. The result is cached to disk so
    repeat runs skip the (slow) text parse.
    """
    import anndata as ad

    if os.path.exists(cache_path):
        print(f"  using cached stereo-seq AnnData: {cache_path}")
        return cache_path

    sim = os.path.join(dataset_dir, "simulated_data")
    print("  building stereo-seq AnnData from raw simulation files (slow parse)...")
    counts = pd.read_csv(os.path.join(sim, "combined_spatial_count.txt"),
                         sep="\t", index_col=0)
    gt = pd.read_csv(os.path.join(sim, "combined_spot_clusters.txt"),
                     sep="\t", index_col=0)
    common = counts.index.intersection(gt.index)
    counts, gt = counts.loc[common], gt.loc[common]
    if n_spots and len(common) > n_spots:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(len(common), size=n_spots, replace=False))
        counts, gt = counts.iloc[idx], gt.iloc[idx]
    counts = counts.loc[:, counts.sum(axis=0) > 0]
    gt = gt.loc[:, gt.sum(axis=0) > 0]
    a = ad.AnnData(X=counts.to_numpy(dtype=np.float32),
                   obs=pd.DataFrame(index=counts.index.astype(str)),
                   var=pd.DataFrame(index=counts.columns.astype(str)))
    a.var["highly_variable"] = True
    a.uns["dataset_name"] = "stereoseq_mouse_brain_li2023_sim_norm_mm"
    a.uns["ground_truth"] = gt.astype(float)
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    a.write_h5ad(cache_path)
    print(f"  cached stereo-seq AnnData ({a.shape}, k={gt.shape[1]}) -> {cache_path}")
    return cache_path


# Datasets built on the fly: name -> callable(dataset_dir, cache_path) -> data_path.
LOADERS = {
    "stereoseq_mouse_brain_li2023_sim_norm_mm": build_stereoseq_adata,
}

# Per-dataset overrides.
#   * spot_overrides["w"]: spOT-NMF usage-entropy, TUNED PER DATASET (w-sweep in
#     scripts/benchmark/). It is the key deconvolution knob -- small w gives
#     peaked (near one-hot) usages suited to crisp single-cell-type panels
#     (STARmap), larger w gives softer graded mixtures suited to genuinely mixed
#     spots (MOB, seqFISH). A single w does NOT transfer across datasets: w=0.01
#     is best for STARmap but collapses MOB/seqFISH accuracy, while w=0.05-0.1
#     recovers/beats cNMF on the mixed slides. This mirrors the manuscript's
#     per-dataset hyperparameter tuning.
#   * n_iter/max_iter caps: only the large stereo-seq slide needs lighter
#     settings (applied to ALL methods, so comparisons stay fair).
OVERRIDES = {
    "Dataset10_STARmap_li2022_sim_norm_mm": dict(spot_overrides={"w": 0.01}),
    "Dataset4_seqFISH_li2022_sim_norm_mm": dict(spot_overrides={"w": 0.1}),
    "MOB_dance_sim_norm_mm": dict(spot_overrides={"w": 0.05}),
    "stereoseq_mouse_brain_li2023_sim_norm_mm": dict(
        n_iter_override=6, spot_max_iter=40, spot_max_iter_inner=80,
        spot_overrides={"w": 0.01}),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_root", default=str(DEFAULT_ROOT),
                    help="directory holding <dataset>/adata_spatial.h5ad")
    ap.add_argument("--out", default=str(REPO / "results" / "benchmark_all"))
    ap.add_argument("--quick", action="store_true",
                    help="fast smoke test (few replicates / iterations)")
    ap.add_argument("--n_hvf", type=int, default=None,
                    help="override the per-dataset top-N overdispersed gene count")
    ap.add_argument("--datasets", nargs="*", default=None,
                    help="subset of dataset folder names to run (default: all)")
    args = ap.parse_args()

    root = Path(args.data_root)
    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    selected = DATASETS
    if args.datasets:
        wanted = set(args.datasets)
        selected = [(d, n) for d, n in DATASETS if d in wanted]

    cache_dir = out_root / "_cache"
    for dataset, n_hvf in selected:
        print("\n" + "#" * 78)
        print(f"# {dataset}")
        print("#" * 78)
        if dataset in LOADERS:
            dataset_dir = root / dataset
            if not dataset_dir.is_dir():
                print(f"  SKIP: {dataset_dir} not found")
                continue
            try:
                data_path = Path(LOADERS[dataset](
                    str(dataset_dir),
                    str(cache_dir / f"{dataset}_adata_spatial.h5ad")))
            except Exception as e:
                print(f"  ERROR building {dataset}: {e}")
                import traceback
                traceback.print_exc()
                continue
        else:
            data_path = root / dataset / "adata_spatial.h5ad"
            if not data_path.exists():
                print(f"  SKIP: {data_path} not found")
                continue
        n = args.n_hvf if args.n_hvf is not None else n_hvf
        try:
            run_benchmark(
                str(data_path), str(out_root / dataset),
                quick=args.quick, n_hvf=n, dataset_label=dataset,
                **OVERRIDES.get(dataset, {}),
            )
        except Exception as e:  # keep going across datasets
            print(f"  ERROR on {dataset}: {e}")
            import traceback
            traceback.print_exc()

    # Build the combined view from every per-dataset summary.csv present under
    # out_root -- this accumulates results across incremental invocations (e.g.
    # running one extra dataset later merges into the same combined table).
    summary_files = sorted(glob.glob(str(out_root / "*" / "summary.csv")))
    frames = []
    for f in summary_files:
        df = pd.read_csv(f)
        if "dataset" not in df.columns:
            df.insert(0, "dataset", Path(f).parent.name)
        frames.append(df)
    if not frames:
        print("\nNo datasets ran successfully.")
        return

    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(out_root / "summary_all.csv", index=False)

    # Per-dataset winner (best PCC_mean) and how often each algorithm wins.
    print("\n" + "=" * 78)
    print("COMBINED SUMMARY  (PCC_mean per algorithm x dataset; higher better)")
    print("=" * 78)
    pivot = combined.pivot(index="algorithm", columns="dataset", values="PCC_mean")
    with pd.option_context("display.width", 220, "display.float_format", "{:.3f}".format):
        print(pivot)
        print("\nMean PCC_mean across datasets:")
        print(pivot.mean(axis=1).sort_values(ascending=False).to_string(float_format="{:.3f}".format))

    winners = combined.loc[combined.groupby("dataset")["PCC_mean"].idxmax(),
                           ["dataset", "algorithm", "PCC_mean"]]
    print("\nPer-dataset winner (PCC_mean):")
    with pd.option_context("display.float_format", "{:.3f}".format):
        print(winners.to_string(index=False))
    print(f"\nWin counts:\n{winners['algorithm'].value_counts().to_string()}")

    try:
        make_combined_plot(pivot, out_root)
        print(f"\nCombined plot: {out_root / 'summary_all.png'}")
    except Exception as e:
        print(f"(combined plot skipped: {e})")
    print(f"\nAll results written to: {out_root}")


def make_combined_plot(pivot, out_root):
    """Grouped bar chart of PCC_mean per algorithm across datasets."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    algos = list(pivot.index)
    datasets = list(pivot.columns)
    x = np.arange(len(datasets))
    width = 0.8 / max(len(algos), 1)

    fig, ax = plt.subplots(figsize=(2.5 * len(datasets) + 4, 5))
    for i, algo in enumerate(algos):
        ax.bar(x + i * width, pivot.loc[algo].values, width, label=algo)
    ax.set_xticks(x + width * (len(algos) - 1) / 2)
    ax.set_xticklabels([d.replace("_sim_norm_mm", "") for d in datasets],
                       rotation=15, ha="right", fontsize=8)
    ax.set_ylabel("PCC_mean (higher better)")
    ax.set_title("Deconvolution accuracy across manuscript datasets")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(out_root / "summary_all.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
