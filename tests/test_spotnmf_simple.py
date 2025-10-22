# tests/test_spotnmf_simple.py
from __future__ import annotations

from pathlib import Path
import numpy as np

def test_spotnmf_pipeline():
    import spotnmf as spot  # noqa: F401

    # Paths
    root = Path(__file__).resolve().parents[1]
    data_file = root / "data" / "test_data" / "dataset10_adata_spatial.h5ad"
    sample_name = "SAMPLE1_k5"
    results_dir = root / "data" / "test_results" / sample_name
    results_dir.mkdir(parents=True, exist_ok=True)

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
    assert isinstance(res, dict) and "adata" in res, "run_spotnmf should return dict with 'adata'."
    out_adata = res["adata"]

    # Plot usage
    spot.pl.plot_usage(adata=out_adata, results_dir=str(results_dir), sample_name="SAMPLE1")

    # Ensure something was written
    assert any(results_dir.iterdir()), f"No files written under {results_dir}"
