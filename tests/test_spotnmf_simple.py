# tests/test_spotnmf_simple.py
from __future__ import annotations

from pathlib import Path

def test_spotnmf_pipeline(tmp_path):
    import spotnmf as spot  # noqa: F401

    # Paths
    root = Path(__file__).resolve().parents[1]
    data_file = root / "data" / "test_data" / "dataset10_adata_spatial.h5ad"
    sample_name = "SAMPLE1_k5"
    results_dir = tmp_path
    # run_experiment writes its outputs under <results_dir>/<sample_name>.
    results_path = results_dir / sample_name

    # Load
    assert data_file.exists(), f"Missing test data at: {data_file}"

    # === Read Data === #
    adata = spot.io.read_adata(
        data_path=data_file,
        data_mode="h5ad"
    )


    # === Model Parameters === #
    model_params = {
        "lr": 0.001,         # Learning rate
        "h": 0.01,           # H regularization
        "w": 0.01,           # W regularization
        "eps": 0.05,         # Epsilon
        "normalize_rows": True,
    }

    # === Run Factorization === #
    res = spot.cli.run_experiment(
        adata_spatial=adata,
        k=5,                        # Number of ranks
        sample_name=sample_name,
        results_dir=str(results_dir),
        genome="mm10",
        annotate=False,
        plot=False,
        network=False,
        is_visium=True,
        model_params=model_params,
    )

    # Basic checks
    assert isinstance(res, dict), "run_experiment should return a dict."
    for key in ("topics_per_spot", "genes_per_topic", "adata"):
        assert key in res, f"run_experiment result missing '{key}'."

    out_adata = res["adata"]
    topics = res["topics_per_spot"]

    # Usage matrix: one row per spot, one column per topic (k=5).
    assert topics.shape == (out_adata.n_obs, 5), f"Unexpected usage shape: {topics.shape}"

    # Deconvolution outputs should be written to disk.
    assert (results_path / f"topics_per_spot_{sample_name}.csv").exists()
    assert (results_path / f"genescores_per_topic_{sample_name}.csv").exists()

    # Spatial plotting should produce a figure without error.
    spot.pl.plot_spatial_all_topics(
        out_adata,
        rf_usages=topics,
        results_dir_path=str(results_path),
        title_name=sample_name,
        is_show=False,
    )
    assert (results_path / f"topics_plot_{sample_name}.pdf").exists()
