# API reference

Import the package as:

```python
import spotnmf
```

The API is organised by analysis stage, mirroring a typical spOT-NMF workflow:
read data → select genes → factorize → score & annotate programs → analyse
niche networks → evaluate → visualise. Click any function for its full signature
and documentation.

## Reading & writing — `spotnmf.io`

```{eval-rst}
.. currentmodule:: spotnmf.io

.. autosummary::
   :toctree: generated
   :nosignatures:

   read_adata
   read_spatial_data
   read_visium_hd
   read_visium_hd_sample
   load_spatial_image
   read_gmt_file
   load_experiment_result
   get_ground_truth
   save_ranked_genes
   save_list
   check_dir
```

## Highly variable genes — `spotnmf.hvg`

```{eval-rst}
.. currentmodule:: spotnmf.hvg

.. autosummary::
   :toctree: generated
   :nosignatures:

   compute_overdispersed_genes
   compute_overdispersed_genes_batches
   save_hvg_list
   load_hvg_list
```

## Deconvolution model — `spotnmf.models`

```{eval-rst}
.. currentmodule:: spotnmf.models

.. autosummary::
   :toctree: generated
   :nosignatures:

   spotnmf
   run_spotnmf
   seed_all
```

## Marker genes & scoring — `spotnmf.gscore`

```{eval-rst}
.. currentmodule:: spotnmf.gscore

.. autosummary::
   :toctree: generated
   :nosignatures:

   calculate_marker_genes_topics_df
   compute_glm_coefficients
```

## Annotation — `spotnmf.annotate`

```{eval-rst}
.. currentmodule:: spotnmf.annotate

.. autosummary::
   :toctree: generated
   :nosignatures:

   compute_genesets_annotation
   list_genesets
   annotate_with_benchmark
   sum_cell_types
   mean_cell_types
   annot_corr_heatmap
   benchmark_corr_silverstandard
```

## Pathway enrichment — `spotnmf.enrichment`

```{eval-rst}
.. currentmodule:: spotnmf.enrichment

.. autosummary::
   :toctree: generated
   :nosignatures:

   run_topics_pathway_enrichment
   program_gprofiler
   order_genesets
   plot_geneset_pval_heatmap
   plot_geneset_pval_clustermap
```

## Niche networks — `spotnmf.niche_networks`

```{eval-rst}
.. currentmodule:: spotnmf.niche_networks

.. autosummary::
   :toctree: generated
   :nosignatures:

   compute_pairwise_stats
   generate_node_attributes
   build_network_graph
   detect_communities_infomap
   get_node_positions
   plot_network_graph
   plot_network_analysis
   calculate_outgoing_and_incoming_connections
   plot_connection_heatmap
```

## Evaluation & benchmarking — `spotnmf.eval`

```{eval-rst}
.. currentmodule:: spotnmf.eval

.. autosummary::
   :toctree: generated
   :nosignatures:

   annotate_programs_by_ground_truth
   get_annotation_from_corr
   get_ranking_score
```

## Plotting — `spotnmf.pl`

```{eval-rst}
.. currentmodule:: spotnmf.pl

.. autosummary::
   :toctree: generated
   :nosignatures:

   plot_spatial_topic
   plot_spatial_all_topics
   plot_spatial_all_topics_aggr
   plot_df_heatmap
   plot_benchmark_methods_topics
```

## Utilities — `spotnmf.utils`

```{eval-rst}
.. currentmodule:: spotnmf.utils

.. autosummary::
   :toctree: generated
   :nosignatures:

   reference_dataset
   compute_ground_cost
   normalize_tensor
   entropy
   entropy_dual_loss
   ot_dual_loss
   early_stop
   clean_mixed_gene_names
```

## Command-line interface — `spotnmf.cli`

The package installs a `spotnmf` console command. These functions back its
subcommands and can also be called directly.

```{eval-rst}
.. currentmodule:: spotnmf.cli

.. autosummary::
   :toctree: generated
   :nosignatures:

   run_experiment
   plot_programs
   annotate_programs
   plot_networks
```
