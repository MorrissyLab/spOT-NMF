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

1. **Install PyTorch** (CPU or CUDA) for your platform (see [pytorch.org](https://pytorch.org)). Examples:

```bash
# CPU-only
pip install torch --index-url https://download.pytorch.org/whl/cpu
# CUDA 11.8 (Linux/Windows with NVIDIA GPUs)
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
```

2. **Install spOT-NMF**:

```bash
pip install spot-nmf
```

3. **Verify the CLI**:

```bash
spotnmf --help
```

> **Conda users:**
>
> ```bash
> conda create -n spotnmf python=3.12
> conda activate spotnmf
> # install torch as above, then:
> pip install spot-nmf
> ```

---

## 🚀 Quick Start

### Command Line

Full pipeline (deconvolution → annotation → spatial plots):

```bash
spotnmf spotnmf \
  --sample_name SAMPLE1 \
  --adata_path ./data/sample1.h5ad \
  --results_dir ./results \
  --k 5
```

Other commands:

```bash
spotnmf deconvolve --sample_name SAMPLE1 --adata_path ./data/sample1.h5ad --results_dir ./results --k 5
spotnmf plot       --sample_name SAMPLE1 --adata_path ./data/sample1.h5ad --results_dir ./results
spotnmf annotate   --sample_name SAMPLE1 --results_dir ./results --genome GRCh38
spotnmf network    --sample_name SAMPLE1 --results_dir ./results --usage_threshold 0 --n_bins 1000 --edge_threshold 0.199
```

### Python / Notebooks

```python
import spotnmf
from spotnmf.io import read_adata
from spotnmf.models import run_spotnmf
from spotnmf.pl import plot_usage

adata = read_adata("./data/sample1.h5ad")
res = run_spotnmf(adata, k=5, sample_name="SAMPLE1", results_dir="./results")
plot_usage(adata=res["adata"], results_dir="./results", sample_name="SAMPLE1")
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
