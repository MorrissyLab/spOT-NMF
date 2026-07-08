"""Hyperparameter sweep for the spOT-NMF backend on the mosaicMPI benchmark input.

Runs a single spOT-NMF factorization per config (no consensus) on the exact
normalized matrix mosaicMPI feeds its factorizer, and scores per-spot
deconvolution accuracy against the ground-truth cell types. Prints configs
ranked by PCC. Uses CUDA automatically when available.

    python scripts/benchmark/sweep_spotnmf.py
"""
import warnings, os, contextlib, itertools, sys, time
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, torch
import mosaicmpi, spotnmf.mosaic
sys.path.insert(0, os.path.dirname(__file__))
from benchmark_mosaic import match_and_score, _row_normalize

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
print(f"device: {DEVICE}  (torch {torch.__version__})")

adata = ad.read_h5ad("data/test_data/dataset10_adata_spatial.h5ad")
adata.obs_names = ["spot" + str(i) for i in range(adata.n_obs)]
adata.var_names = adata.var_names.astype(str)
gt = adata.uns["ground_truth"]; gt.index = adata.obs_names
gtp = _row_normalize(gt.astype(float)); k = gt.shape[1]

ds = mosaicmpi.Dataset(adata.copy())
ds.select_hvf(feature_list=list(adata.var_names), score_type="odscore")
run = ds.initialize_cnmf(cnmf_output_dir="results/benchmark_mosaic/_sweep",
                         cnmf_name="s", kvals=[k], n_iter=1, algorithm="spotnmf")
norm = ad.read_h5ad(run.paths["normalized_counts"])
X = np.asarray(norm.X.todense()) if hasattr(norm.X, "todense") else np.asarray(norm.X)
fac = mosaicmpi.factorization.get_factorizer("spotnmf")


def evaluate(params):
    t = time.time()
    with open(os.devnull, "w") as n, contextlib.redirect_stdout(n), contextlib.redirect_stderr(n):
        _sp, W = fac(X, k, 0, var_names=list(norm.var.index), params=params)
    W = np.asarray(W)
    usage = pd.DataFrame(W, index=norm.obs.index, columns=[str(i + 1) for i in range(k)])
    sc = match_and_score(usage, gtp)
    return sc, time.time() - t


# Grid. Fixed optimizer/iterations; vary the OT-relevant knobs plus the
# README-recommended values (w=0.01, normalize_rows=True).
base = dict(device=DEVICE, optim_name="adam", lr=1e-2,
            max_iter=50, max_iter_inner=150)
grid = {
    "eps": [0.02, 0.05, 0.1],
    "h": [1e-3, 1e-2, 1e-1],
    "w": [1e-3, 1e-2],                # 1e-2 = README value
    "normalize_rows": [False, True],  # True = README value
    "cost": ["cosine", "euclidean"],
}
keys = list(grid)
combos = list(itertools.product(*[grid[k_] for k_ in keys]))
print(f"{len(combos)} configs to evaluate\n")

rows = []
for i, vals in enumerate(combos):
    p = dict(base); p.update(dict(zip(keys, vals)))
    sc, sec = evaluate(p)
    rows.append({**dict(zip(keys, vals)),
                 "PCC": sc["PCC_mean"], "SCC": sc["SCC_mean"],
                 "RMSE": sc["RMSE"], "JSD": sc["JSD"],
                 "disp": sc["usage_dispersion"], "sec": sec})
    print(f"[{i+1}/{len(combos)}] PCC={sc['PCC_mean']:.3f} RMSE={sc['RMSE']:.3f} "
          f"disp={sc['usage_dispersion']:.3f} ({sec:.1f}s)  {dict(zip(keys, vals))}")

df = pd.DataFrame(rows).sort_values("PCC", ascending=False)
df.to_csv("results/benchmark_mosaic/sweep_results.csv", index=False)
print("\n=== TOP 10 by PCC ===")
with pd.option_context("display.width", 200, "display.float_format", "{:.3f}".format):
    print(df.head(10).to_string(index=False))
print("\nsaved: results/benchmark_mosaic/sweep_results.csv")
