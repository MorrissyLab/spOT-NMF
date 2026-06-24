# spOT-NMF

**Optimal Transport-Based Matrix Factorization for Accurate Deconvolution of Spatial Transcriptomics**
*Abdelkareem, A.O. et al.(2025)*

`spOT-NMF` is a Python package for unsupervised deconvolution and discovery of gene programs in spatial transcriptomics. It integrates **Optimal Transport (OT)** into a non-negative matrix factorization (NMF) framework, enabling robust topic modeling, high-resolution spatial deconvolution, and rich biological annotation.

This package supports the analyses in:
**spOT-NMF: Optimal Transport-Based Matrix Factorization for Accurate Deconvolution of Spatial Transcriptomics** — bioRxiv (2025). DOI: **10.1101/2025.08.02.668292**

---

## ✨ Key Features

* **OT-NMF Deconvolution**: Reference-free topic modeling with OT-regularized NMF.
* **HVG Selection**: Flexible, batch-aware highly variable gene selection.
* **Biological Annotation**: Automated enrichment and gene-set overlap of inferred programs.
* **Spatial Visualization**: Publication-quality spatial plots for topic/program usage.
* **Scalable & Modular**: Built for large datasets and multi-sample workflows.
* **CLI & Python API**: Run from the command line or import in notebooks.

---

## 📦 Installation

`spOT-NMF` requires **Python ≥ 3.12**. We recommend [**uv**](https://docs.astral.sh/uv/) for a fast,
reproducible setup. PyTorch is installed separately so you can pick the build (CPU or CUDA) for your platform.

### Recommended: uv

```bash
# 1. Create and activate an isolated environment (uv fetches Python 3.12 if needed)
uv venv --python 3.12
# Linux/macOS:  source .venv/bin/activate
# Windows:      .venv\Scripts\activate

# 2. Install PyTorch for your platform (see pytorch.org)
#    CPU-only:
uv pip install torch --index-url https://download.pytorch.org/whl/cpu
#    CUDA 12.x (Linux/Windows with NVIDIA GPUs):
#    uv pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# 3. Install spOT-NMF
uv pip install spot-nmf
```

### Alternative: pip

```bash
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install spot-nmf
```

### From source (development)

```bash
git clone https://github.com/MorrissyLab/spOT-NMF.git
cd spOT-NMF
uv venv --python 3.12
uv pip install torch --index-url https://download.pytorch.org/whl/cpu
uv pip install -e ".[dev]"     # editable install with test dependencies
uv run pytest -q               # run the test suite
```

### Verify the install

```bash
spotnmf --help
```

> If no GPU is available, spOT-NMF automatically runs on CPU.

---

## 🚀 Quick Start

### Command Line

Full pipeline (deconvolution → annotation → spatial plots → networks):

```bash
spotnmf spotnmf \
  --sample_name SAMPLE1 \
  --adata_path ./data/sample1.h5ad \
  --data_mode h5ad \
  --results_dir ./results \
  --k 5 \
  --genome GRCh38
```

> **`--data_mode`** selects how the input is read: `h5ad` for a single AnnData
> `.h5ad` file, `visium` (the default) for a Space Ranger output directory, or
> `visium_hd` for Visium HD. Pass `--data_mode h5ad` whenever `--adata_path`
> points to a `.h5ad` file.

Other commands:

```bash
spotnmf deconvolve --sample_name SAMPLE1 --adata_path ./data/sample1.h5ad --data_mode h5ad --results_dir ./results --k 5
spotnmf plot       --sample_name SAMPLE1 --adata_path ./data/sample1.h5ad --data_mode h5ad --results_dir ./results
spotnmf annotate   --sample_name SAMPLE1 --results_dir ./results --genome GRCh38
spotnmf network    --sample_name SAMPLE1 --results_dir ./results --usage_threshold 0 --n_bins 1000 --edge_threshold 0.199
```

> The `network` command reuses the per-spot usages written by `deconvolve`. On
> small datasets no topic pairs may pass `--n_bins` / `--edge_threshold`; in that
> case it prints a notice and skips plotting — lower the thresholds to force a graph.

### Python / Notebooks

```python
from pathlib import Path
import spotnmf as spot

# === Configuration === #
DATA_PATH = Path("data/test_data/dataset10_adata_spatial.h5ad")
RESULTS_DIR = Path(r"/data/test_results/")
SAMPLE_NAME = "TestSample"
GENOME = "mm10"

# === Read Data === #
adata = spot.io.read_adata(
    data_path=DATA_PATH,
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
results = spot.cli.run_experiment(
    adata_spatial=adata,
    k=5,                        # Number of ranks
    sample_name=SAMPLE_NAME,
    results_dir=str(RESULTS_DIR),
    genome=GENOME,
    annotate=False,
    plot=False,
    network=False,
    is_visium=True,
    model_params=model_params,
)

# === Annotate Programs === #
spot.cli.annotate_programs(
    results_dir=str(RESULTS_DIR),
    sample_name=SAMPLE_NAME,
    genome=GENOME,
)

```

---

## ⚙️ CLI Overview

| Command      | Description                                                  |
| ------------ | ------------------------------------------------------------ |
| `spotnmf`    | Full pipeline: deconvolution → annotation → spatial plotting |
| `deconvolve` | Run OT-NMF and save results                                  |
| `plot`       | Visualize spatial topic/program usage                        |
| `annotate`   | Enrich and annotate gene programs                            |
| `network`    | Visualize niche networks based on topic interactions         |

Run `spotnmf <command> --help` for per-command options.

---

## 📁 Outputs

* `topics_per_spot_{sample}.csv` — topic/program usage per spot
* `genescores_per_topic_{sample}.csv` — gene scores per topic
* `ranked_genescores_{sample}.csv` — ranked marker genes per topic
* Pathway enrichment and gene-set overlap tables
* Spatial plots & QC visualizations
* Network plots of topic–topic interactions

---

## 🔬 Reproducibility (Manuscript Notebooks)

The **main** branch provides the reusable software package.
The original Jupyter notebooks used to reproduce manuscript figures are maintained in the **`manuscript`** branch:

```bash
git fetch origin
git checkout manuscript
```

Notebooks are in:

```
scripts/manuscript_notebooks/
```

Use **`manuscript`** to regenerate paper figures; use **`main`** for running the package on your data.

---

## 🧾 Citation

Please cite:

> Abdelkareem, A.O., Gill, G.S., Manoharan, V.T., Verhey, T.B., & Morrissy, A.S.
> **spOT-NMF: Optimal Transport-Based Matrix Factorization for Accurate Deconvolution of Spatial Transcriptomics.**
> *bioRxiv* (2025). [https://doi.org/10.1101/2025.08.02.668292](https://doi.org/10.1101/2025.08.02.668292)

```bibtex
@article{abdelkareem2025spotnmf,
  title   = {spOT-NMF: Optimal Transport-Based Matrix Factorization for Accurate Deconvolution of Spatial Transcriptomics},
  author  = {Abdelkareem, Aly O. and Gill, Gurveer S. and Manoharan, Varsha Thoppey and Verhey, Theodore B. and Morrissy, A. Sorana},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.08.02.668292},
  url     = {https://www.biorxiv.org/content/10.1101/2025.08.02.668292v1},
  note    = {Preprint}
}
```


---

## 🤝 Contributing

We welcome ideas, bug reports, and feature requests—**please open a GitHub Issue**:
[https://github.com/MorrissyLab/spOT-NMF/issues](https://github.com/MorrissyLab/spOT-NMF/issues)

---

## 📜 License

GPL-3.0. See **LICENSE** for details.

---

## 💬 Support

Questions or need help? Open an Issue:
[https://github.com/MorrissyLab/spOT-NMF/issues](https://github.com/MorrissyLab/spOT-NMF/issues)
