# spOT-NMF

**Optimal Transport-Based Matrix Factorization for Accurate Deconvolution in Spatial Transcriptomics**
*Aly O. Abdelkareem et al., 2025*

---

spOT-NMF is a Python package for unsupervised deconvolution and discovery of gene programs in spatial transcriptomics data. By integrating **Optimal Transport (OT)** into a non-negative matrix factorization (NMF) framework, spOT-NMF enables robust topic modeling, high-resolution spatial deconvolution, and comprehensive pathway and gene set annotation.

This package powers the analyses in the paper:
**Optimal Transport-Based Matrix Factorization for Accurate Deconvolution in Spatial Transcriptomics**
*Abdelkareem, A.O. et al., 2025*

## Key Features

* **OT-NMF Deconvolution**: Unsupervised topic modeling of spatial transcriptomics data via OT-regularized NMF, capturing gene programs and their spatial usage.
* **Highly Variable Gene (HVG) Selection**: Flexible strategies for gene selection and batch-aware analysis.
* **Program Annotation**: Automated enrichment and annotation of inferred topics with biological pathways and custom gene sets.
* **Spatial Visualization**: High-quality spatial mapping of deconvolved topics/programs.
* **Scalability**: Efficient for large spatial datasets and multi-sample (aggregated) analysis.

---

## Installation

Clone the repository and set up the environment:

```bash
git clone https://github.com/MorrissyLab/spOT-NMF.git
cd spOT-NMF
```

Install dependencies (recommended via Conda):

```bash
# Create environment
conda create -n spotnmf python=3.12.3
conda activate spotnmf

# Install major dependencies
conda install -c conda-forge scanpy python-igraph leidenalg scikit-learn statsmodels pandas pygam scipy=1.12
pip install rbo distinctipy gprofiler-official==1.0.0 fastcluster==1.2.6

# Install PyTorch (adjust CUDA version as needed)
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
```

---

## Quick Start

**Run full deconvolution and annotation:**

```bash
python cli.py spotnmf --sample_name SAMPLE1 --adata_path ./data/sample1.h5ad --results_dir ./results --k 5
```

* This runs topic modeling (OT-NMF), saves spatial usage, gene scores, performs annotation, and generates spatial plots.

**Basic deconvolution only:**

```bash
python cli.py deconvolve --sample_name SAMPLE1 --adata_path ./data/sample1.h5ad --results_dir ./results --k 5
```

**Plot spatial topic usage:**

```bash
python cli.py plot --sample_name SAMPLE1 --adata_path ./data/sample1.h5ad --results_dir ./results
```

**Annotate gene programs (enrichment and gene set overlap):**

```bash
python cli.py annotate --sample_name SAMPLE1 --results_dir ./results --genome GRCh38
```

**Plot niche network plots:**

```bash
python cli.py network --sample_name SAMPLE1 --results_dir ./results --usage_threshold 0 --n_bins 1000 --edge_threshold 0.199
```

---

## Command-Line Interface (CLI)

The main CLI (`cli.py`) supports the following commands:

| Command      | Description                                                    |
| ------------ | -------------------------------------------------------------- |
| spotnmf      | Full pipeline: deconvolution, annotation, and spatial plotting |
| deconvolve   | Only run OT-NMF and save results                               |
| plot         | Visualize spatial topic/program usage                          |
| annotate     | Annotate gene programs with pathway and gene set enrichment    |
| network      | Visualize niche networks                                       |

### Key Arguments

* `--sample_name`: Name for this analysis run (required)
* `--adata_path`: Path to input AnnData (`.h5ad`) file (required for deconvolve/spotnmf/plot)
* `--results_dir`: Output directory (required)
* `--k`: Number of topics/components (required for deconvolve/spotnmf)
* `--genome`: Reference genome label (default: mm10)
* `--data_mode`: Data type (`visium`, `visium_hd`, etc.)
* `--is_xeno`: Flag for xenograft data
* `--is_aggr`: Flag for aggregated libraries
* `--hvg_file`: Precomputed highly variable genes (optional)
* `--usage_threshold`: Usage threshold
* `--n_bins`: Number of bins
* `--edge_threshold`: Edge threshold
* `--annot_file`: Annotation file

**Model parameters** (customizable):

* `--lr`: Learning rate
* `--h`, `--w`, `--eps`: OT-NMF regularization parameters
* `--normalize_rows`: Normalize input matrix rows

*See* `python cli.py --help` *for all options.*

---

## Outputs

* **Deconvolution**:

  * `topics_per_spot_{sample}.csv`: Topic (program) usage per spot
  * `genescores_per_topic_{sample}.csv`: Marker gene scores per topic
  * `ranked_genescores_{sample}.csv`: Ranked gene lists

* **Annotation**:

  * Pathway enrichment results for each gene set
  * Gene set overlap tables for user-defined references

* **Visualization**:

  * Spatial maps of topics/programs (per sample)
  * QC and statistics for HVG selection
  * Niche network plots showing interactions between topics

---

## Example

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

## Citing spOT-NMF

If you use this package, please cite:

> Abdelkareem, A.O., et al. Optimal Transport-Based Matrix Factorization for Accurate Deconvolution in Spatial Transcriptomics. 2025. \[Preprint/BioRxiv/DOI\:XXXXXX]
> [GitHub: MorrissyLab/spOT-NMF](https://github.com/MorrissyLab/spOT-NMF)

---

## Contributing

Pull requests, issues, and feature suggestions are welcome!
See [CONTRIBUTING.md](./CONTRIBUTING.md) or open an issue to get started.

---

## License

GPL-3.0 License. See [LICENSE](./LICENSE) for details.

---

For questions or support, please open a GitHub issue

---

**spot-NMF: Optimal Transport-Based Matrix Factorization for Accurate Deconvolution of Spatial Transcriptomics**

---
