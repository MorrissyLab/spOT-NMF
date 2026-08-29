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
conda install -c conda-forge scanpy python-igraph leidenalg scikit-learn statsmodels pandas pygam scipy=1.12 adjusttext
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
spotnmf deconvolve --sample_name SAMPLE1 --adata_path ./data/sample1.h5ad --data_mode h5ad --results_dir ./results --k 5
spotnmf plot       --sample_name SAMPLE1 --adata_path ./data/sample1.h5ad --data_mode h5ad --results_dir ./results
spotnmf annotate   --sample_name SAMPLE1 --results_dir ./results --genome GRCh38
spotnmf network    --sample_name SAMPLE1 --results_dir ./results
```

> The `network` command reuses the per-spot usages written by `deconvolve`.

#### Niche inference

Niches are inferred by testing every program pair for co-occurrence enrichment
against a **spatial null**, rather than by applying fixed cutoffs:

```bash
spotnmf network --sample_name SAMPLE1 --results_dir ./results \
  --presence_method otsu --null torus --n_perm 1000 \
  --fdr 0.05 --min_log2_oe 1.0
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--presence_method` | `otsu` | How a program is called present in a spot. `otsu` thresholds each program on its own usage distribution; `quantile` reproduces the legacy fixed cutoff. |
| `--presence_quantile` | `0.90` | Quantile used by `--presence_method quantile`. Previously hard-coded and unreachable. |
| `--null` | `torus` | Null model. `torus` rigidly translates each program's map, preserving its spatial autocorrelation. `label` ignores spatial structure and is anticonservative. |
| `--n_perm` | `1000` | Permutations. |
| `--fdr` | `0.05` | Benjamini–Hochberg FDR across the P(P−1) program pairs. |
| `--min_log2_oe` | `1.0` | Minimum log2 observed/expected co-occurrence — 1.0 is two-fold. |
| `--min_prevalence_frac` | `0.01` | Minimum fraction of spots per program. Replaces the absolute `--n_bins`. |

**Why these replaced `--n_bins 1000 --edge_threshold 0.199`.** The old rule zeroed
every program below its own 90th percentile, which forced *every* program to be
present in exactly 10% of spots. A conditional co-occurrence probability of
0.199 was therefore ~2× the chance expectation — a reasonable target, but one
whose meaning changed silently if the quantile moved. On the CRC test matrix the
edges that rule retained span log2(O/E) 1.00–2.67, i.e. it *was* a two-fold
enrichment filter. `--min_log2_oe 1.0` states that target directly and holds it
fixed regardless of how presence is called. Separately, `--n_bins 1000` was an
absolute spot count: combined with the 10% presence cap it required ≥10,000
spots, so it silently produced empty graphs on standard Visium.

To reproduce previously published networks exactly, pass `--legacy_network`
together with the old `--usage_threshold`, `--n_bins` and `--edge_threshold`.

#### Robustness sweep

```bash
spotnmf network_sweep --sample_name SAMPLE1 --results_dir ./results
```

Writes `{sample}_niche_parameter_sweep.csv`: edge count, niche count, edge-set
Jaccard and partition ARI against the default, across presence rules, quantiles,
FDR levels and enrichment cutoffs. Read `ari`/`jaccard` **together with**
`n_edges` — both agreement measures saturate trivially on near-empty graphs, so
the result to report is the widest region where agreement is high *and* the
network is non-trivial. The `degenerate` column flags rows too sparse to interpret.

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
python cli.py plot --sample_name SAMPLE1 --adata_path ./data/sample1.h5ad --results_dir ./results
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
  build a program–program co-occurrence network from the usage matrix, test it against a
  spatial null, detect niches with Infomap, and render the network + connection heatmap
  (step-by-step and one-call).

---

## ⚙️ CLI Overview

| Command      | Description                                                  |
| ------------ | ------------------------------------------------------------ |
| `spotnmf`    | Full pipeline: deconvolution → annotation → spatial plotting |
| `deconvolve` | Run OT-NMF and save results                                  |
| `plot`       | Visualize spatial topic/program usage                        |
| `annotate`   | Enrich and annotate gene programs                            |
| `network`    | Infer and visualize niche networks from topic co-occurrence   |
| `network_sweep` | Parameter-robustness sweep for niche inference             |

Run `spotnmf <command> --help` for per-command options.

---

## 📁 Outputs

* `topics_per_spot_{sample}.csv` — topic/program usage per spot
* `genescores_per_topic_{sample}.csv` — gene scores per topic
* `ranked_genescores_{sample}.csv` — ranked marker genes per topic
* `{sample}_presence_calls.csv` — per-program presence threshold, realised prevalence, and whether the call was clamped
* `{sample}_cooccurrence_stats.csv` — every program pair with log2(O/E), log odds ratio, Jaccard, permutation z, p and BH q
* `{sample}_niche_parameter_sweep.csv` — robustness sweep (from `network_sweep`)
* Pathway enrichment and gene-set overlap tables
* Spatial plots & QC visualizations
* Network plots of topic–topic interactions

---

## 🔬 Reproducibility (Manuscript Notebooks)

The **main** branch provides the reusable software package.
The original Jupyter notebooks used to reproduce manuscript figures are maintained in the **`manuscript`** branch:

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
