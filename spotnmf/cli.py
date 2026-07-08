"""
spotnmf Experiment Runner

This script provides an interface for running spatial transcriptomics experiments 
using the spotnmf package. It supports topic modeling-based deconvolution, 
gene set enrichment, and annotation.

Functionality:
    - run_experiment: Main function for processing spatial data and saving outputs.
    - plot: Visualize inferred topics on spatial tissue maps.
    - network: Plot networks of gene programs.
    - annotate_programs: Perform enrichment and annotation of gene programs.
    - main: Command-line interface to execute various modes (deconvolution, plotting, annotation).
"""


import argparse
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import time
from datetime import datetime
import pandas as pd
pd.options.display.float_format = '{:f}'.format
from typing import Union

from spotnmf.models import run_spotnmf, run_spotnmf_consensus
from spotnmf.gscore import calculate_marker_genes_topics_df
from spotnmf import io, pl, annotate, enrichment, hvg, niche_networks


def plot_programs(results_dir, sample_name, adata_spatial, is_visium=True, genome=None, is_xenograft=False, is_aggr = True):
    """Plot topic usages spatially for a given sample.

    Reads the saved ``topics_per_spot_{sample_name}.csv`` usages and dispatches to
    the appropriate plotting routine depending on the data type. For xenograft
    Visium data, usages are scaled by the per-spot admixture ratio for the selected
    genome before plotting.

    Args:
        results_dir (str): Root directory containing the sample's result files.
        sample_name (str): Sample identifier; used to locate inputs and title plots.
        adata_spatial (anndata.AnnData): AnnData object with spatial coordinates and,
            for xenografts, an "admix" column in ``obs``.
        is_visium (bool): If True, use the Visium aggregate plotting paths; otherwise
            use the generic spatial-topics plot. Defaults to True.
        genome (str): Reference genome used to select the admixture ratio for
            xenograft scaling (e.g. "GRCh38", "mm10"). Defaults to None.
        is_xenograft (bool): If True, scale usages by the admixture ratio.
            Defaults to False.
        is_aggr (bool): If True, use the aggregate multi-sample plot; otherwise use
            the manuscript-style per-sample plot. Defaults to True.
    """
    print("Plotting spatial topics")
    results_path = os.path.join(results_dir, sample_name)
    rf_usages = pd.read_csv(os.path.join(results_path, f"topics_per_spot_{sample_name}.csv"), index_col=0)

    if is_visium:
        if is_xenograft:
            ratio_map = {
                "GRCh38": 1 - adata_spatial.obs["admix"],
                "mm10": adata_spatial.obs["admix"]
            }
            admix_ratios = ratio_map.get(genome)
            if admix_ratios is not None:
                rf_usages = rf_usages.mul(admix_ratios, axis=0)

        if is_aggr:
            pl.plot_spatial_all_topics_aggr(
                adata_spatial,
                rf_usages=rf_usages,
                results_dir_path=results_path,
                title_name=sample_name,
                same_legend=False,
                plot_topic=True,
                filter_th = 0.9
            )
        else:
            pl.plot_spatial_all_topics_aggr_manuscript(
                adata_spatial,
                rf_usages=rf_usages,
                results_dir_path=results_path,
                title_name=sample_name,
                same_legend=False,
                COLS=5, ROWS=10,
                filter_th = 0.9
            )
            
    else:
        pl.plot_spatial_all_topics(
            adata_spatial,
            rf_usages=rf_usages,
            results_dir_path=results_path,
            title_name=sample_name,
            is_show=False,
        )


def annotate_programs(results_dir, sample_name, genome):
    """Run gene set enrichment and annotate topics for a given sample.

    Loads the saved ``genescores_per_topic_{sample_name}.csv`` gene scores, runs
    pathway enrichment for a set of standard gene sets (GO:CC, GO:BP, KEGG, REAC),
    and computes gene-set annotations for each geneset available for the genome.

    Args:
        results_dir (str): Root directory containing the sample's result files and
            where enrichment/annotation outputs are written.
        sample_name (str): Sample identifier; used to locate inputs and title
            outputs.
        genome (str): Reference genome (e.g. "GRCh38", "mm10") used for enrichment
            and to select available gene sets.
    """
    print("Annotating programs with pathway enrichment")
    results_path = os.path.join(results_dir, sample_name)
    gene_scores_df = pd.read_csv(os.path.join(results_path, f"genescores_per_topic_{sample_name}.csv"), index_col=0)

    for gene_set in ["GO:CC", "GO:BP", "KEGG", "REAC"]:
        print(f"Running enrichment for {gene_set}")
        enrichment.run_topics_pathway_enrichment(
            gene_scores_df,
            gene_set=gene_set,
            results_dir_path=results_path,
            top_n_features=1000,
            genome=genome,
            experiment_title=sample_name,
        )

    print("Matching gene sets for annotation")
    for gene_set in annotate.list_genesets(genome=genome):
        annotate.compute_genesets_annotation(
            gene_scores_df,
            gene_set,
            results_dir_path=results_path,
            max_top_genes=100,
            ranking_method="rboext",
            experiment_title=f"{sample_name}_{gene_set}",
        )


def plot_networks(results_dir: str, sample_name: str, usage_threshold: Union[float, int], n_bins: int, edge_threshold: float, annot_file: Union[str, None]):
    """Plot niche networks for a given sample.

    Thin wrapper around ``niche_networks.plot_network_analysis`` that builds and
    plots co-occurrence/niche networks from the sample's saved results.

    Args:
        results_dir (str): Root directory containing the sample's result files.
        sample_name (str): Sample identifier.
        usage_threshold (float | int): Minimum topic usage for a spot to count
            toward the network.
        n_bins (int): Number of spatial bins used when building the network.
        edge_threshold (float): Minimum edge weight to retain in the network.
        annot_file (str | None): Optional path to an annotation file mapping topics
            to labels; None to skip annotation.
    """
    print("Plotting niche networks.")

    niche_networks.plot_network_analysis(
        results_dir=results_dir,
        sample_name=sample_name,
        usage_threshold=usage_threshold,
        n_bins=n_bins,
        edge_threshold=edge_threshold,
        annot_file=annot_file
    )


def run_experiment(
    adata_spatial,
    k: int,
    sample_name: str,
    results_dir: str,
    genome=None,
    filter_genes=True,
    hvg_file=None,
    annotate=False,
    plot=False,
    network=False,
    is_visium=True,
    is_aggr = False,
    is_xenograft=False,
    usage_threshold: Union[float, int] = 0,
    n_bins: int = 1000,
    edge_threshold: float = 0.199,
    annot_file: Union[str, None] = None,
    model_params={},
    consensus: bool = False,
    n_iter: int = 10,
    **kwargs,
):
    """Run a complete spotnmf experiment: gene selection, model training, gene ranking, and optional plotting/annotation/networks.

    Optionally selects highly variable (overdispersed) genes, fits the spotnmf topic
    model, saves the factorization matrices and per-topic gene scores, and then
    optionally annotates topics, plots spatial maps, and plots niche networks. Also
    writes a timing/loss record.

    Args:
        adata_spatial (anndata.AnnData): AnnData object with the spatial expression
            data to factorize.
        k (int): Number of topics/components to learn.
        sample_name (str): Sample identifier; used for the output directory and file
            names.
        results_dir (str): Root directory under which the sample's results are saved.
        genome (str): Reference genome (e.g. "GRCh38", "mm10"); if None and needed
            for annotation, it is read from ``adata_spatial.uns``. Defaults to None.
        filter_genes (bool): If True, perform overdispersed-gene selection before
            fitting. Defaults to True.
        hvg_file (str): Optional path to a precomputed highly variable gene list; if
            given, gene selection is loaded rather than computed. Defaults to None.
        annotate (bool): If True, run enrichment/annotation after fitting. Defaults
            to False.
        plot (bool): If True, plot spatial topic maps after fitting. Defaults to
            False.
        network (bool): If True, plot niche networks after fitting. Defaults to
            False.
        is_visium (bool): Whether the data is Visium; passed through to plotting.
            Defaults to True.
        is_aggr (bool): Whether the data is aggregated across libraries; controls
            batch-mode gene selection and aggregate plotting. Defaults to False.
        is_xenograft (bool): Whether the dataset is a xenograft model; passed through
            to plotting. Defaults to False.
        usage_threshold (float | int): Minimum topic usage used for network plotting.
            Defaults to 0.
        n_bins (int): Number of spatial bins used for network plotting. Defaults to
            1000.
        edge_threshold (float): Minimum edge weight for network plotting. Defaults to
            0.199.
        annot_file (str | None): Optional annotation file for network plotting.
            Defaults to None.
        model_params (dict): Extra keyword arguments forwarded to ``run_spotnmf``
            (or ``run_spotnmf_consensus`` when ``consensus=True``). Defaults to an
            empty dict.
        consensus (bool): If True, run consensus spOT-NMF -- ``n_iter`` replicate
            factorizations clustered into consensus programs with an
            optimal-transport usage refit (:func:`spotnmf.models.run_spotnmf_consensus`)
            -- instead of a single factorization. More robust and, on the STARmap
            deconvolution benchmark, more accurate, at ``n_iter``x the runtime.
            Defaults to False.
        n_iter (int): Number of replicate factorizations when ``consensus=True``.
            Defaults to 10.
        **kwargs: Additional keyword arguments forwarded to the gene-selection
            routines.

    Returns:
        dict: The factorization results (e.g. "topics_per_spot", "genes_per_topic")
        plus an "adata" entry holding the mutated ``adata_spatial`` with learned
        usages/spectra for plotting.

    Raises:
        ValueError: If fewer highly variable genes are selected than ``k``.
    """
    start_time = time.time()
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] Running experiment for '{sample_name}' with k = {k}")

    results_path = io.check_dir(os.path.join(results_dir, sample_name))

    # Filter genes
    if filter_genes:
        if hvg_file:
            print(f"Loading overdispersed genes from '{hvg_file}'")
            overdispersed_genes = hvg.load_hvg_list(hvg_file)

        else:
            if is_aggr:
                print(f"Computing overdispersed genes using batch mode")
                overdispersed_genes = hvg.compute_overdispersed_genes_batches(
                    adata_spatial,
                    batch_keys=["sample_id"],
                    union_agg=False,
                    **kwargs
                )
            else:
                print("Computing overdispersed genes (global mode)")
                adata_spatial, overdispersed_genes = hvg.compute_overdispersed_genes(
                    adata_spatial,
                    save_dir=results_path,
                    is_show=False,
                    n_top_genes=None,
                    **kwargs
                )
                # Save gene stats for inspection
                adata_spatial.var.to_csv(os.path.join(results_path, f"gene_stats_{sample_name}.csv"))
                
  

            # Save final list of selected genes
            hvg.save_hvg_list(overdispersed_genes, os.path.join(results_path, f"top_genes_{sample_name}.csv"))

        adata_spatial.var['highly_variable'] = adata_spatial.var.index.isin(overdispersed_genes)

        n_hvg = int(adata_spatial.var['highly_variable'].sum())
        if n_hvg < int(k):
            raise ValueError(
                f"Only {n_hvg} highly variable gene(s) were selected, which is fewer than "
                f"k={k} requested factors. The factorization needs more genes than factors. "
                f"Relax the gene-selection thresholds (e.g. increase 'alpha'), provide an "
                f"--hvg_file, or use a dataset with more expression variability."
            )

    # Run topic model (single factorization, or consensus over n_iter replicates).
    if consensus:
        print(f"Running consensus spOT-NMF over {n_iter} replicates")
        results, losses = run_spotnmf_consensus(
            adata_spatial, components=k, n_iter=n_iter, **model_params
        )
    else:
        results, losses = run_spotnmf(adata_spatial, components=k, **model_params)

    # Save matrices
    for key, df in results.items():
        df.to_csv(os.path.join(results_path, f"{key}_{sample_name}.csv"))

    # Calculate and save gene scores
    print("Calculating and saving gene scores")
    gene_scores_df = calculate_marker_genes_topics_df(adata_spatial, results["topics_per_spot"])
    gene_scores_df.to_csv(os.path.join(results_path, f"genescores_per_topic_{sample_name}.csv"))

    io.save_ranked_genes(gene_scores_df, os.path.join(results_path, f"ranked_genescores_{sample_name}.csv"))
    if "genes_per_topic" in results:
        io.save_ranked_genes(results["genes_per_topic"], os.path.join(results_path, f"ranked_genes_{sample_name}.csv"))

    # Annotate topics
    if annotate:
        if not genome:
            genome = adata_spatial.uns.params["genome"]
        annotate_programs(results_dir, sample_name, genome)

    # Plot spatial maps
    if plot:
        plot_programs(results_dir, sample_name, adata_spatial, is_visium=is_visium, genome=genome, is_xenograft=is_xenograft, is_aggr = is_aggr)

    # Plot networks
    if network:
        plot_networks(results_dir, sample_name, usage_threshold=usage_threshold, n_bins=n_bins, edge_threshold=edge_threshold, annot_file=annot_file)

    # Save timing
    duration = time.time() - start_time
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Completed in {duration:.2f} seconds.")
    with open(os.path.join(results_path, f"time_{sample_name}.txt"), "w") as f:
        f.write(f"{duration}\t{losses[-1]}\t")

    # Return the factorization matrices alongside the mutated AnnData, which now
    # holds the learned usages/spectra (obsm["W_OT"], uns["H_OT"]) for plotting.
    return {**results, "adata": adata_spatial}


def main():
    """Command-line entry point for running spotnmf experiments (spotnmf, deconvolve, plot, annotate, network)."""
    parser = argparse.ArgumentParser(description="Run spatial transcriptomics experiments with spotnmf.")
    parser.add_argument("run_type", choices=["spotnmf", "deconvolve", "plot", "annotate", "network"], help="Type of operation to perform.")
    parser.add_argument("--sample_name", required=True, help="Sample identifier.")
    parser.add_argument("--results_dir", required=True, help="Directory for saving results.")

    parser.add_argument("--adata_path", help="Path to AnnData file.")
    parser.add_argument("--k", type=int, help="Number of topics/components.")
    parser.add_argument("--genome", default="mm10", help="Reference genome (e.g., GRCh38, mm10).")
    parser.add_argument("--data_mode", default="visium", help="Data mode (e.g., visium, visium_hd).")
    parser.add_argument("--bin_size", default=16, help="Bin Size (e.g., 2, 8, 16, 24).")
    parser.add_argument("--is_xeno", action="store_true", help="Whether the dataset is a xenograft model.")
    parser.add_argument("--is_aggr", action="store_true", help="Whether data is aggregated across libraries.")
    parser.add_argument("--select_sample", default=None, help="Subset a specific sample.")
    parser.add_argument("--hvg_file", default=None, help="Precomputed highly variable genes file.")
    parser.add_argument("--usage_threshold", type=float, default=0, help="Usage threshold.")
    parser.add_argument("--n_bins", type=int, default=1000, help="Number of bins.")
    parser.add_argument("--edge_threshold", type=float, default=0.199, help="Edge threshold.")
    parser.add_argument("--annot_file", default=None, help="Annotation file.")

    parser.add_argument("--lr", type=float, default=0.01, help="Learning rate.")
    parser.add_argument("--h", type=float, default=0.01, help="H Regularizer parameter.")
    parser.add_argument("--w", type=float, default=0.01, help="W Regularizer parameter.")
    parser.add_argument("--eps", type=float, default=0.02, help="Entropy.")
    parser.add_argument("--normalize_rows", action=argparse.BooleanOptionalAction,
                        default=True, help="Normalize rows of input matrix (use --no-normalize_rows to disable).")

    parser.add_argument("--consensus", action="store_true",
                        help="Run consensus spOT-NMF (cluster --n_iter replicate "
                             "factorizations, then OT usage refit) instead of a single run.")
    parser.add_argument("--n_iter", type=int, default=10,
                        help="Number of replicate factorizations for --consensus (default: 10).")

    args = parser.parse_args()

    # Validate required arguments
    if args.run_type in {"spotnmf", "deconvolve"} and not (args.adata_path and args.k):
        parser.error("--adata_path and --k are required for 'spotnmf' or 'deconvolve'.")
    if args.run_type == "annotate" and not args.genome:
        parser.error("--genome is required for 'annotate'.")
    if args.run_type == "plot" and not args.adata_path:
        parser.error("--adata_path is required for 'plot'.")
    if args.run_type == "network" and (args.usage_threshold is None or args.n_bins is None or args.edge_threshold is None):
        parser.error("--usage_threshold, --n_bins, and --edge_threshold are required for 'network'.")

    is_visium = args.data_mode in {"visium", "visium_hd"}

    # Only commands that operate on the data need to load it; 'annotate' and
    # 'network' work purely from previously-saved result files.
    adata_spatial = None
    if args.run_type in {"spotnmf", "deconvolve", "plot"}:
        adata_spatial = io.read_adata(
            data_path=args.adata_path,
            data_mode=args.data_mode,
            genome=args.genome,
            is_aggr=args.is_aggr,
            is_xenograft=args.is_xeno,
            select_sample=args.select_sample,
            bin_size=args.bin_size
        )

    model_params = {
        "lr": args.lr,
        "h": args.h,
        "w": args.w,
        "eps": args.eps,
        "normalize_rows": args.normalize_rows,
    }

    if args.run_type == "spotnmf":
        run_experiment(
            adata_spatial, args.k, args.sample_name, args.results_dir,
            genome=args.genome, hvg_file=args.hvg_file,
            annotate=True, plot=True, network=True,
            is_visium=is_visium, is_xenograft=args.is_xeno, is_aggr=args.is_aggr,
            model_params=model_params,
            consensus=args.consensus, n_iter=args.n_iter,
        )
    elif args.run_type == "deconvolve":
        run_experiment(
            adata_spatial, args.k, args.sample_name, args.results_dir,
            genome=args.genome, hvg_file=args.hvg_file,
            annotate=False, plot=False, network=False,
            is_visium=is_visium, is_xenograft=args.is_xeno, is_aggr=args.is_aggr,
            model_params=model_params,
            consensus=args.consensus, n_iter=args.n_iter,
        )
    elif args.run_type == "plot":
        plot_programs(args.results_dir, args.sample_name, adata_spatial, is_visium=is_visium, genome=args.genome, is_xenograft=args.is_xeno, is_aggr=args.is_aggr)
    elif args.run_type == "annotate":
        annotate_programs(args.results_dir, args.sample_name, genome=args.genome)
    elif args.run_type == "network":
        plot_networks(args.results_dir, args.sample_name, args.usage_threshold, args.n_bins, args.edge_threshold, annot_file=args.annot_file)


if __name__ == "__main__":
    main()
