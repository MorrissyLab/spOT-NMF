# spOT-NMF

[![PyPI version](https://img.shields.io/pypi/v/spot-nmf.svg)](https://pypi.org/project/spot-nmf/)
[![Python versions](https://img.shields.io/pypi/pyversions/spot-nmf.svg)](https://pypi.org/project/spot-nmf/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Publish to PyPI](https://github.com/MorrissyLab/spOT-NMF/actions/workflows/publish.yml/badge.svg)](https://github.com/MorrissyLab/spOT-NMF/actions/workflows/publish.yml)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2025.08.02.668292-red.svg)](https://doi.org/10.1101/2025.08.02.668292)
[![Documentation Status](https://readthedocs.org/projects/spot-nmf/badge/?version=latest)](https://spot-nmf.readthedocs.io/en/latest/)

📖 **Documentation:** <https://spot-nmf.readthedocs.io>

**Optimal Transport-Based Matrix Factorization for Accurate Deconvolution of Spatial Transcriptomics**
*Abdelkareem, A.O. et al.(2025)*

`spOT-NMF` is a Python package for unsupervised deconvolution and discovery of gene programs in spatial transcriptomics. It integrates **Optimal Transport (OT)** into a non-negative matrix factorization (NMF) framework, enabling robust topic modeling, high-resolution spatial deconvolution, and rich biological annotation.

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
> small datasets no topic pairs may pass `--n_bins` / `--edge_threshold`; lower
> the thresholds to force a graph.

### Python / Notebooks

```python
import spotnmf as spot

adata = spot.io.read_adata("data/sample1.h5ad", data_mode="h5ad")

results = spot.cli.run_experiment(
    adata_spatial=adata,
    k=5,                       # number of programs
    sample_name="SAMPLE1",
    results_dir="./results",
    genome="GRCh38",
)
```

Model hyperparameters (`lr`, `h`, `w`, `eps`, …) default to tuned values; pass a
`model_params` dict to override them. See the **[full pipeline tutorial](docs/source/tutorials/full_pipeline.ipynb)**
for a complete walkthrough (HVG selection, annotation, spatial plots, and validation).

### Consensus mode (more robust, more accurate)

By default spOT-NMF runs a **single** factorization. In **consensus** mode it fits
several replicates, clusters them into consensus programs (cNMF-style), and refits
the per-spot usages with a fixed-spectra **optimal-transport** solve — keeping
spOT-NMF's OT usage geometry rather than reverting to least squares. On the packaged
simulated STARmap benchmark it was the top performer across accuracy metrics
(see `scripts/benchmark/`).

Enable it with `--consensus` on the CLI or `consensus=True` in Python; `--n_iter`
(default 10) controls the replicate count.

```bash
spotnmf deconvolve --sample_name SAMPLE1 --adata_path ./data/sample1.h5ad \
  --data_mode h5ad --results_dir ./results --k 5 --consensus --n_iter 10
```

> Consensus runs `n_iter` factorizations, so it takes roughly `n_iter`× longer
> than a single run. A GPU is recommended for larger datasets.

---

## 📓 Tutorials

Fully worked, **well-commented notebooks** run on the packaged example dataset with
figures pre-rendered, so GitHub displays them directly in the browser.

* **[Full pipeline tutorial →](docs/source/tutorials/full_pipeline.ipynb)** — data
  loading, HVG selection, OT-NMF deconvolution, spatial mapping, marker extraction, and
  validation against ground-truth cell types, end-to-end.
* **[Proximity (niche) networks →](docs/source/tutorials/proximity_networks.ipynb)** —
  build a program–program co-occurrence network from the usage matrix, detect niches with
  Infomap, and render the network + connection heatmap (step-by-step and one-call).

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

## 🤝 Contributing & Support

Questions, ideas, bug reports, and feature requests—**please open a GitHub Issue**:
[https://github.com/MorrissyLab/spOT-NMF/issues](https://github.com/MorrissyLab/spOT-NMF/issues)

---

## 📜 License

GPL-3.0. See **LICENSE** for details.
