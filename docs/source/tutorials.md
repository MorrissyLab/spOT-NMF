# Tutorials

```{toctree}
:maxdepth: 1

tutorials/full_pipeline
tutorials/proximity_networks
tutorials/mosaicmpi-integration
```

- **Full pipeline** — load data, select highly variable genes, run the OT-NMF
  deconvolution, map programs spatially, recover marker genes, and validate the
  results against ground-truth cell types, end-to-end on the bundled example dataset.
- **Proximity (niche) networks** — turn the per-spot usage matrix into a
  program–program co-occurrence network, detect niches with Infomap, and render the
  network and in-/cross-niche connection heatmap, both step-by-step and via the
  one-call `plot_network_analysis` / `spotnmf network` API.
- **mosaicMPI integration** — use spOT-NMF as the factorization backend inside
  [mosaicMPI](https://github.com/MorrissyLab/mosaicMPI) by setting
  `algorithm="spotnmf"`, gaining consensus factorization, stability/error rank
  selection, usage heatmaps, and multi-dataset integration on a scanpy Visium slide.
