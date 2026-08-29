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
import numpy as np
from typing import Union

from spotnmf import run_spotnmf
from spotnmf.gscore import calculate_marker_genes_topics_df
from spotnmf import io, pl, annotate, enrichment, hvg, niche_networks


def plot_programs(results_dir, sample_name, adata_spatial, is_visium=True, genome=None, is_xenograft=False, is_aggr = True):
    """
    Plot topic usages spatially for a given sample.
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
    """
    Run gene set enrichment and annotate topics for a given sample.
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


def plot_networks(results_dir: str, sample_name: str, annot_file: Union[str, None] = None,
                  presence_method: str = "otsu", presence_quantile: float = 0.90,
                  fdr: float = 0.05, min_log2_oe: float = 1.0,
                  min_prevalence_frac: float = 0.01, null: str = "torus",
                  n_perm: int = 1000, seed: int = 0, legacy: bool = False,
                  usage_threshold: Union[float, int] = 0, n_bins: int = 0,
                  edge_threshold: Union[float, None] = None,
                  adata_spatial=None):
    """Plot niche networks for a given sample.

    Thin wrapper around ``niche_networks.plot_network_analysis``. Spatial
    coordinates are taken from ``adata_spatial`` when supplied, and otherwise
    parsed from Visium HD barcodes in the usage index; without either, the
    spatial null degrades to label permutation.

    Args:
        results_dir (str): Root directory containing the sample's result files.
        sample_name (str): Sample identifier.
        annot_file (str | None): Optional path to an annotation file mapping topics
            to labels; None to skip annotation.
        presence_method (str): Presence rule -- "otsu", "mixture" or "quantile".
        presence_quantile (float): Quantile for presence_method="quantile".
        fdr (float): Benjamini-Hochberg false discovery rate.
        min_log2_oe (float): Minimum log2 observed/expected enrichment.
        min_prevalence_frac (float): Minimum fraction of spots per program.
        null (str): Null model -- "torus", "block" or "label".
        n_perm (int): Number of permutations.
        seed (int): Random seed.
        legacy (bool): Reproduce the original fixed-threshold behaviour.
        usage_threshold (float | int): Legacy usage threshold (legacy mode only).
        n_bins (int): Legacy absolute spot count per edge (legacy mode only).
        edge_threshold (float | None): Legacy co-occurrence cutoff (legacy mode only).
        adata_spatial (anndata.AnnData | None): Optional AnnData carrying
            ``obsm["spatial"]`` coordinates for non-Visium-HD platforms.
    """
    print("Plotting niche networks.")

    coords = None
    if adata_spatial is not None and "spatial" in getattr(adata_spatial, "obsm", {}):
        coords = adata_spatial.obsm["spatial"]

    niche_networks.plot_network_analysis(
        results_dir=results_dir,
        sample_name=sample_name,
        annot_file=annot_file,
        presence_method=presence_method,
        presence_quantile=presence_quantile,
        fdr=fdr,
        min_log2_oe=min_log2_oe,
        min_prevalence_frac=min_prevalence_frac,
        null=null,
        n_perm=n_perm,
        coords=coords,
        seed=seed,
        legacy=legacy,
        usage_threshold=usage_threshold,
        n_bins=n_bins,
        edge_threshold=edge_threshold,
    )


def sweep_networks(results_dir: str, sample_name: str, annot_file: Union[str, None] = None,
                   n_perm: int = 200, null: str = "torus", seed: int = 0,
                   adata_spatial=None):
    """Write a parameter-robustness sweep for a sample's niche network.

    Runs :func:`~spotnmf.niche_stats.parameter_sweep` across presence rules,
    presence quantiles, FDR levels and enrichment cutoffs, and saves the table
    to ``<results_dir>/<sample_name>/<sample_name>_niche_parameter_sweep.csv``.
    This is what turns "why these thresholds?" into a reportable answer: the
    table shows how the edge set and the niche partition respond across the
    whole grid rather than at one chosen point.

    Args:
        results_dir (str): Root directory containing the sample's result files.
        sample_name (str): Sample identifier.
        annot_file (str | None): Optional annotation file mapping topics to labels.
        n_perm (int): Permutations per grid cell. Defaults to 200.
        null (str): Null model. Defaults to "torus".
        seed (int): Random seed. Defaults to 0.
        adata_spatial (anndata.AnnData | None): Optional AnnData carrying
            ``obsm["spatial"]`` coordinates for non-Visium-HD platforms.

    Returns:
        pandas.DataFrame: The sweep table.
    """
    from spotnmf import niche_stats as nstats

    print("Running niche parameter sweep.")
    results_path = os.path.join(results_dir, sample_name)
    usage = pd.read_csv(
        os.path.join(results_path, f"topics_per_spot_{sample_name}.csv"), index_col=0
    )
    if annot_file is not None:
        annot = pd.read_csv(annot_file)
        if "Annotation" not in annot.columns or "Program" not in annot.columns:
            raise ValueError("Annotation file must contain 'Program' and 'Annotation' columns.")
        annot_dict = dict(zip(annot["Program"], annot["Annotation"]))
        usage.columns = [
            col.replace("ot_", "ot") + f"_{annot_dict[col]}" if col in annot_dict
            else col.replace("ot_", "ot")
            for col in usage.columns
        ]

    coords = None
    if adata_spatial is not None and "spatial" in getattr(adata_spatial, "obsm", {}):
        coords = adata_spatial.obsm["spatial"]

    sweep = nstats.parameter_sweep(
        usage, coords=coords, n_perm=n_perm, null=null, seed=seed, sample=sample_name
    )
    out_file = os.path.join(results_path, f"{sample_name}_niche_parameter_sweep.csv")
    sweep.to_csv(out_file, index=False)
    print(f"Wrote {len(sweep)} parameter combinations to {out_file}")

    stable = sweep[(~sweep["degenerate"]) & (sweep["ari"] >= 0.8)]
    print(
        f"{len(stable)} of {len(sweep)} settings are non-degenerate with ARI >= 0.8 "
        "against the default. Read ari/jaccard together with n_edges: both agreement "
        "measures saturate trivially on near-empty graphs."
    )
    return sweep


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
    is_aggr = True,
    is_xenograft=False,
    usage_threshold: Union[float, int] = 0,
    n_bins: int = 1000,
    edge_threshold: float = 0.199,
    annot_file: Union[str, None] = None,
    presence_method: str = "otsu",
    presence_quantile: float = 0.90,
    fdr: float = 0.05,
    min_log2_oe: float = 1.0,
    min_prevalence_frac: float = 0.01,
    null: str = "torus",
    n_perm: int = 1000,
    seed: int = 0,
    legacy_network: bool = False,
    model_params={},
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
        presence_method (str): Niche presence rule -- "otsu", "mixture" or "quantile".
        presence_quantile (float): Quantile for presence_method="quantile". Defaults to 0.90.
        fdr (float): Benjamini-Hochberg false discovery rate for niche edges. Defaults to 0.05.
        min_log2_oe (float): Minimum log2 observed/expected co-occurrence. Defaults to 1.0.
        min_prevalence_frac (float): Minimum fraction of spots per program. Defaults to 0.01.
        null (str): Null model -- "torus", "block" or "label". Defaults to "torus".
        n_perm (int): Permutations for the niche null. Defaults to 1000.
        seed (int): Random seed for permutations and Infomap. Defaults to 0.
        legacy_network (bool): Reproduce the original fixed-threshold niche
            inference. Defaults to False.
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

    # Run topic model
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
        plot_networks(
            results_dir, sample_name, annot_file=annot_file,
            presence_method=presence_method, presence_quantile=presence_quantile,
            fdr=fdr, min_log2_oe=min_log2_oe, min_prevalence_frac=min_prevalence_frac,
            null=null, n_perm=n_perm, seed=seed, legacy=legacy_network,
            usage_threshold=usage_threshold, n_bins=n_bins,
            edge_threshold=edge_threshold, adata_spatial=adata_spatial,
        )

    # Save timing
    duration = time.time() - start_time
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Completed in {duration:.2f} seconds.")
    with open(os.path.join(results_path, f"time_{sample_name}.txt"), "w") as f:
        f.write(f"{duration}\t{losses[-1]}\t")


def main():
    """
    Command-line interface for running spotnmf experiments.
    """
    parser = argparse.ArgumentParser(description="Run spatial transcriptomics experiments with spotnmf.")
    parser.add_argument("run_type", choices=["spotnmf", "deconvolve", "plot", "annotate", "network", "network_sweep"], help="Type of operation to perform.")
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
    parser.add_argument("--annot_file", default=None, help="Annotation file.")

    # --- Niche inference -------------------------------------------------
    parser.add_argument("--presence_method", default="otsu", choices=["otsu", "mixture", "quantile"],
                        help="How to call a program present in a spot. 'otsu' (default) picks a "
                             "per-program threshold by maximizing between-class variance; "
                             "'quantile' reproduces the legacy fixed cutoff.")
    parser.add_argument("--presence_quantile", type=float, default=0.90,
                        help="Quantile used when --presence_method quantile. Default 0.90 -- the "
                             "value that used to be hard-coded and unreachable from the CLI.")
    parser.add_argument("--fdr", type=float, default=0.05,
                        help="Benjamini-Hochberg false discovery rate for edge selection.")
    parser.add_argument("--min_log2_oe", type=float, default=1.0,
                        help="Minimum log2 observed/expected co-occurrence. Default 1.0 = two-fold, "
                             "which is what the legacy --edge_threshold 0.199 encoded implicitly.")
    parser.add_argument("--min_prevalence_frac", type=float, default=0.01,
                        help="Both programs must be present in at least this fraction of spots. "
                             "Replaces the absolute --n_bins, which needed >=10,000 spots to pass.")
    parser.add_argument("--null", default="torus", choices=["torus", "block", "label"],
                        help="Null model for the permutation test. 'torus' (default) rigidly "
                             "translates each program's map, preserving its spatial "
                             "autocorrelation. 'label' ignores spatial structure and is "
                             "anticonservative -- use only as a sanity floor.")
    parser.add_argument("--n_perm", type=int, default=1000, help="Number of permutations.")
    parser.add_argument("--seed", type=int, default=0, help="Random seed for permutations and Infomap.")
    parser.add_argument("--legacy_network", action="store_true",
                        help="Reproduce the original fixed-threshold niche inference "
                             "(90th-percentile presence + --edge_threshold), for backwards "
                             "comparison only. Applies no null model or multiple-testing correction.")
    parser.add_argument("--usage_threshold", type=float, default=0, help="Legacy usage threshold.")
    parser.add_argument("--n_bins", type=int, default=1000, help="Legacy absolute bin count per edge.")
    parser.add_argument("--edge_threshold", type=float, default=0.199, help="Legacy edge threshold.")

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
    if args.run_type == "network" and args.legacy_network and (
        args.usage_threshold is None or args.n_bins is None or args.edge_threshold is None
    ):
        parser.error(
            "--usage_threshold, --n_bins and --edge_threshold are required with --legacy_network."
        )

    is_visium = args.data_mode in {"visium", "visium_hd"}

    adata_spatial = None
        
    if args.adata_path:
        if args.data_mode == "visium_hd":
            adata_spatial = io.read_visium_hd(
                adata_path=args.adata_path,
                bin_size=args.bin_size,
                genome=args.genome,
                is_aggr=args.is_aggr,
                is_xenograft=args.is_xeno
            )
        else:
            adata_spatial = io.read_spatial_data(
                adata_path=args.adata_path,
                genome=args.genome,
                is_xenograft=args.is_xeno,
                is_aggr=args.is_aggr,
                select_sample=args.select_sample
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
            model_params=model_params
        )
    elif args.run_type == "deconvolve":
        run_experiment(
            adata_spatial, args.k, args.sample_name, args.results_dir,
            genome=args.genome, hvg_file=args.hvg_file,
            annotate=False, plot=False, network=False,
            is_visium=is_visium, is_xenograft=args.is_xeno, is_aggr=args.is_aggr,
            model_params=model_params
        )
    elif args.run_type == "plot":
        plot_programs(args.results_dir, args.sample_name, adata_spatial, is_visium=is_visium, genome=args.genome, is_xenograft=args.is_xeno, is_aggr=args.is_aggr)
    elif args.run_type == "annotate":
        annotate_programs(args.results_dir, args.sample_name, genome=args.genome)
    elif args.run_type == "network_sweep":
        sweep_networks(
            args.results_dir, args.sample_name, annot_file=args.annot_file,
            n_perm=args.n_perm, null=args.null, seed=args.seed,
            adata_spatial=adata_spatial,
        )
    elif args.run_type == "network":
        plot_networks(
            args.results_dir, args.sample_name, annot_file=args.annot_file,
            presence_method=args.presence_method, presence_quantile=args.presence_quantile,
            fdr=args.fdr, min_log2_oe=args.min_log2_oe,
            min_prevalence_frac=args.min_prevalence_frac, null=args.null,
            n_perm=args.n_perm, seed=args.seed, legacy=args.legacy_network,
            usage_threshold=args.usage_threshold, n_bins=args.n_bins,
            edge_threshold=args.edge_threshold, adata_spatial=adata_spatial,
        )


if __name__ == "__main__":
    main()
