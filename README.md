# spOT-NMF

**Optimal Transport-Based Matrix Factorization for Accurate Deconvolution in Spatial Transcriptomics**  
*Abdelkareem, et al. (2025)*

---

`spOT-NMF` is a Python package for unsupervised deconvolution and discovery of gene programs in spatial transcriptomics data. It integrates **Optimal Transport (OT)** into a non-negative matrix factorization (NMF) framework, enabling robust topic modeling, high-resolution spatial deconvolution, and rich biological annotation.

This package powers the analyses in the paper:  
**"spOT-NMF: Optimal Transport-Based Matrix Factorization for Accurate Deconvolution of Spatial Transcriptomics"**  
[**bioRxiv** (2025)](https://www.biorxiv.org/content/10.1101/2025.08.02.668292v1) — [DOI: 10.1101/2025.08.02.668292](https://doi.org/10.1101/2025.08.02.668292)

---

## 🧬 Key Features

- **OT-NMF Deconvolution**: Unsupervised topic modeling of spatial transcriptomics data via OT-regularized NMF.
- **HVG Selection**: Flexible options for selecting highly variable genes, with batch-awareness.
- **Biological Annotation**: Automated enrichment and gene set overlap of inferred gene programs.
- **Spatial Visualization**: High-quality spatial plots of usage and programs.
- **Scalable and Modular**: Designed for large datasets and multi-sample workflows.

---

## 📦 Installation

Clone and install:

```bash
git clone https://github.com/MorrissyLab/spOT-NMF.git
cd spOT-NMF
````

Install dependencies (recommended via Conda):

```bash
conda create -n spotnmf python=3.12.3
conda activate spotnmf

conda install -c conda-forge scanpy python-igraph leidenalg scikit-learn statsmodels pandas pygam scipy=1.12 adjusttext
pip install rbo distinctipy gprofiler-official==1.0.0 fastcluster==1.2.6

pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
```

---

## 🚀 Quick Start

Run full deconvolution and annotation:

```bash
python cli.py spotnmf --sample_name SAMPLE1 --adata_path ./data/sample1.h5ad --results_dir ./results --k 5
```

Basic deconvolution only:

```bash
python cli.py deconvolve --sample_name SAMPLE1 --adata_path ./data/sample1.h5ad --results_dir ./results --k 5
```

Plot topic usage:

```bash
python cli.py plot --sample_name SAMPLE1 --adata_path ./data/sample1.h5ad --results_dir ./results
```

Annotate gene programs:

```bash
python cli.py annotate --sample_name SAMPLE1 --results_dir ./results --genome GRCh38
```

Plot niche networks:

```bash
python cli.py network --sample_name SAMPLE1 --results_dir ./results --usage_threshold 0 --n_bins 1000 --edge_threshold 0.199
```

---

## ⚙️ CLI Overview

| Command      | Description                                                    |
| ------------ | -------------------------------------------------------------- |
| `spotnmf`    | Full pipeline: deconvolution, annotation, and spatial plotting |
| `deconvolve` | Run OT-NMF and save results                                    |
| `plot`       | Visualize spatial topic/program usage                          |
| `annotate`   | Enrich and annotate gene programs                              |
| `network`    | Visualize niche networks based on topic interactions           |

Use `python cli.py --help` for all available options and parameters.

---

## 📁 Outputs

* `topics_per_spot_{sample}.csv`: Usage per topic/program per spot
* `genescores_per_topic_{sample}.csv`: Gene scores per topic
* `ranked_genescores_{sample}.csv`: Sorted marker genes per topic
* Pathway enrichment results and gene set overlap tables
* Spatial plots and QC visualizations
* Network plots of topic-topic interactions

---

## 🔬 Example Run

```bash
python cli.py spotnmf \
    --sample_name Brain1 \
    --adata_path ./data/Brain1.h5ad \
    --results_dir ./results \
    --k 8 \
    --genome GRCh38 \
    --data_mode visium \
    --is_aggr
```

---

## 📖 How to Cite

If you use **spOT-NMF** in your work, please cite:

> Abdelkareem, A.O., Manoharan, V.T., Gill, G.S., Verhey, T.B., & Morrissy, A.S.
> *spOT-NMF: Optimal Transport-Based Matrix Factorization for Accurate Deconvolution of Spatial Transcriptomics.*
> bioRxiv 2025. [https://doi.org/10.1101/2025.08.02.668292](https://doi.org/10.1101/2025.08.02.668292)

```bibtex
@article{abdelkareem2025spotnmf,
  title={spOT-NMF: Optimal Transport-Based Matrix Factorization for Accurate Deconvolution of Spatial Transcriptomics},
  author={Abdelkareem, Aly O and Manoharan, Varsha T and Gill, Gurveer S and Verhey, Theodore B and Morrissy, A Sorana},
  journal={bioRxiv},
  year={2025},
  doi={10.1101/2025.08.02.668292},
  url={https://www.biorxiv.org/content/10.1101/2025.08.02.668292v1}
}
```

---

## 🤝 Contributing

Contributions are welcome! See [CONTRIBUTING.md](./CONTRIBUTING.md) or open an issue.

---

## 📜 License

GPL-3.0 License. See [LICENSE](./LICENSE) for details.

---

For questions or help, open a [GitHub issue](https://github.com/MorrissyLab/spOT-NMF/issues).

