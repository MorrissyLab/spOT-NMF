# spOT-NMF

**Optimal Transport-Based Matrix Factorization for Accurate Deconvolution of Spatial Transcriptomics**

`spOT-NMF` is a Python package for unsupervised deconvolution and discovery of
gene programs in spatial transcriptomics. It integrates **Optimal Transport (OT)**
into a non-negative matrix factorization (NMF) framework, enabling robust topic
modeling, high-resolution spatial deconvolution, and rich biological annotation.

This package supports the analyses in
**spOT-NMF: Optimal Transport-Based Matrix Factorization for Accurate Deconvolution
of Spatial Transcriptomics** — bioRxiv (2025),
DOI: [10.1101/2025.08.02.668292](https://doi.org/10.1101/2025.08.02.668292).

## Key features

- **OT-NMF deconvolution** — reference-free topic modeling with OT-regularized NMF.
- **HVG selection** — flexible, batch-aware highly variable gene selection.
- **Biological annotation** — automated enrichment and gene-set overlap of inferred programs.
- **Spatial visualization** — publication-quality spatial plots for topic/program usage.
- **Scalable & modular** — built for large datasets and multi-sample workflows.
- **CLI & Python API** — run from the command line or import in notebooks.

## Quick start

```bash
# PyTorch is installed separately so you can pick the right build (CPU or CUDA).
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install spot-nmf
```

```python
import spotnmf

adata = spotnmf.io.read(...)          # load spatial data
adata = spotnmf.hvg.select(adata)     # highly variable genes
model = spotnmf.models.spotnmf(...)   # OT-NMF deconvolution
```

See {doc}`installation` for the full setup (including `uv` and GPU builds), and the
{doc}`tutorials` for an end-to-end walkthrough on the bundled example dataset.

```{toctree}
:hidden:
:maxdepth: 2

installation
tutorials
api
```
