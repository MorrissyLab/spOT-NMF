## Code from MosaicMPI 
from types import SimpleNamespace
from typing import Union, Optional, Literal
from collections.abc import Collection
import os
import pandas as pd
import numpy as np
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import traceback

from spotnmf.io import check_dir
from spotnmf.utils import clean_mixed_gene_names

def run_topics_pathway_enrichment(rf_usages, gene_set, results_dir_path, top_n_features=1000, genome='mm10', experiment_title="experiment"):
    """Run pathway enrichment analysis for the topics in ``rf_usages``.

    Queries g:Profiler for each topic's top genes, saves enrichment results and
    a readable annotation summary as CSVs, and writes heatmap and clustermap
    plots of the -log10(p-value) scores to ``results_dir_path``.

    Args:
        rf_usages (pandas.DataFrame): Usage matrix with topics as columns and
            genes as rows.
        gene_set (str): Name of the gene set source to use for enrichment
            analysis.
        results_dir_path (str): Directory path in which to save results.
        top_n_features (int): Number of top-scoring genes per topic to consider
            for enrichment. Defaults to 1000.
        genome (str): Genome identifier; ``'GRCh38'`` maps to human, any other
            value (e.g. ``'mm10'``) maps to mouse. Defaults to ``'mm10'``.
        experiment_title (str): Title for the experiment, used in output file
            naming. Defaults to ``"experiment"``.
    """
    gene_set_safe = gene_set.replace(":", "_")
    rf_usages.index = clean_mixed_gene_names(rf_usages.index, genome)

    try:
        result = program_gprofiler(
            program_df=rf_usages,
            species="hsapiens" if genome == 'GRCh38' else 'mmusculus',
            n_hsg=top_n_features,
            gene_sets=[gene_set],
            min_termsize=10,
            max_termsize=2000
        )
    except Exception as e:
        print(f"Error in program_gprofiler: {e}")
        if len(rf_usages.index) < 500:
            print("Not enough genes for pathway enrichment. Minimum required: 500.")
        return

    output_path = check_dir(os.path.join(results_dir_path, "pathway_results_v7"))
    file_base = f"{experiment_title}_{gene_set_safe}_n{top_n_features}"

    # Save enrichment results
    result.gprofiler_output.to_csv(os.path.join(output_path, f"results_{file_base}.csv"), sep=",")
    result.summary.to_csv(os.path.join(output_path, f"results_summary_{file_base}.csv"), sep=",")

    # Prepare annotation summary with scores
    annotation_summary = {}
    for topic in rf_usages.columns:
        sorted_pathways = result.summary[("-log10pval", topic)].dropna().sort_values(ascending=False)
        annotation_summary[topic] = [f"{x[2]}:{score:.2f}:{x[1]}" for x, score in zip(sorted_pathways.index, sorted_pathways.values)]
    
    # Convert annotation summary to a DataFrame and save
    annotation_df = pd.DataFrame({k: pd.Series(v) for k, v in annotation_summary.items()})
    annotation_df.to_csv(os.path.join(output_path, f"readable_results_summary_{file_base}.csv"), sep=",")

    # Generate heatmap plot of -log10p values
    df = result.summary["-log10pval"].dropna(how="all").fillna(0)
    df = order_genesets(df)
    
    try:
        # Plot heatmap and legend
        fig, figlegend = plot_geneset_pval_heatmap(df=df, plot_title=file_base)
        figlegend.savefig(os.path.join(output_path, f"heatmap_legend_{file_base}.pdf"))
        fig.savefig(os.path.join(output_path, f"heatmap_plot_{file_base}.pdf"))
        plt.close(fig)
        plt.close(figlegend)

        # Plot clustermap, detailed heatmap and legend
        fig, fig_cluster, figlegend = plot_geneset_pval_clustermap(df=df, plot_title=file_base)
        fig.savefig(os.path.join(output_path, f"heatmap_detailed_{file_base}.pdf"))
        fig_cluster.savefig(os.path.join(output_path, f"clustermap_{file_base}.pdf"))
        plt.close(fig)
        plt.close(fig_cluster.figure)  # sns.clustermap returns a ClusterGrid, not a Figure
        plt.close(figlegend)
        plt.close("all")

    except Exception as e:
        print("The plot is huge and it can't be rendered.")
        traceback.print_exc()  # optional: prints detailed error info
        plt.close("all")


def program_gprofiler(program_df: pd.DataFrame,
                      species: Literal["hsapiens", "mmusculus"],
                      n_hsg: int = 1000,
                      gene_sets: Collection[str] = [],
                      no_iea: bool = False,
                      min_termsize: int = 10,
                      max_termsize: int = 2000,
                      batch_size: int = 20,
                      show_progress_bar: bool = True
                      ) -> SimpleNamespace:
    """Run g:Profiler functional enrichment on the top genes of each program.

    For every program (column) in ``program_df``, the top ``n_hsg`` genes by
    score are queried against g:Profiler in batches. Results are combined into a
    single output table and a pivoted summary of enrichment statistics, filtered
    by term size.

    Args:
        program_df (pandas.DataFrame): Program matrix with genes as rows and
            programs as columns. Column names may be single- or multi-level.
        species (Literal["hsapiens", "mmusculus"]): g:Profiler organism to query.
        n_hsg (int): Number of top-scoring genes per program to query. Defaults
            to 1000.
        gene_sets (Collection[str]): Enrichment data sources to query. Defaults
            to an empty list.
        no_iea (bool): Whether to exclude electronically inferred annotations.
            Recorded on the result namespace. Defaults to False.
        min_termsize (int): Minimum term size to retain in the summary. Defaults
            to 10.
        max_termsize (int): Maximum term size to retain in the summary. Defaults
            to 2000.
        batch_size (int): Number of programs to include per g:Profiler query.
            Defaults to 20.
        show_progress_bar (bool): Whether to display a progress bar while
            querying. Defaults to True.

    Returns:
        types.SimpleNamespace: Namespace with the query inputs (``background``,
        ``query``, ``gene_sets``, ``no_iea``), the concatenated
        ``gprofiler_output`` DataFrame (including a ``-log10pval`` column), and a
        pivoted ``summary`` DataFrame of enrichment statistics per program.
    """
    from gprofiler import GProfiler
    result = SimpleNamespace()
    result.background = program_df.dropna(how="all").index.to_list()  # all genes in program_df
    prog_names_str = program_df.columns.map(str)  # gProfiler multi-query only supports string names for queries
    if program_df.columns.nlevels == 1:
        prog_names_str_decoder = {progstr: [prog] for progstr, prog in zip(prog_names_str, program_df.columns)}
    else:
        prog_names_str_decoder = {progstr: list(prog) for progstr, prog in zip(prog_names_str, program_df.columns)}
    prog_level_names = program_df.columns.names

    # result.hsg = program_df.rank(ascending=False) <= n_hsg  # bool dataframe of high scoring genes across programs
    # result.query = {prog_str: genes[genes].index.to_list() for (prog, genes), prog_str in zip(result.hsg.items(), prog_names_str)}

    result.query = {}
    for program in program_df.columns:
        result.query[program] = program_df[program].sort_values(ascending=False).head(n_hsg).index.to_list()

    result.gene_sets = gene_sets
    result.no_iea = no_iea

    gp = GProfiler(return_dataframe=True)

    result.gprofiler_output = []
    batch_query = []
    for i, query in enumerate(tqdm(result.query.keys(), total=len(result.query), unit="program", desc="Querying g:Profiler", disable=not show_progress_bar), start=1):
        batch_query.append(query)
        if i % batch_size == 0 or i == len(result.query):
            batch_result = gp.profile(organism=species, query={q: result.query[q] for q in batch_query},
                                    sources=gene_sets, no_iea=False,
                                    domain_scope="annotated",
                                    measure_underrepresentation=False,
                                    no_evidences=False,
                                    user_threshold =0.05,
                                    significance_threshold_method="g_SCS")
            result.gprofiler_output.append(batch_result)
            batch_query = []
        
    result.gprofiler_output = pd.concat(result.gprofiler_output)
    result.gprofiler_output["-log10pval"] = np.minimum(-np.log10(result.gprofiler_output["p_value"]), 10)

    
    subset = ((result.gprofiler_output["term_size"] <= max_termsize) &
              (result.gprofiler_output["term_size"] >= min_termsize))
    result.summary = result.gprofiler_output[subset].pivot(index=["source", "native", "name", "description", "term_size"], columns="query")
    stats = ["-log10pval", "query_size", "intersection_size"]
    result.summary = result.summary[stats]
    result.summary.columns = pd.MultiIndex.from_tuples([([c[0]] + prog_names_str_decoder[c[1]]) for c in result.summary.columns], names=["stat"] + prog_level_names)

    # conform column order to input dataframe
    if program_df.columns.nlevels == 1:
        sorted_cols = pd.MultiIndex.from_tuples([(stat, prog) for stat in stats for prog in program_df.columns])
    else:
        sorted_cols = pd.MultiIndex.from_tuples([tuple([stat] + list(prog)) for stat in stats for prog in program_df.columns])
    result.summary = result.summary.reindex(columns=sorted_cols)
    return result

def order_genesets(df: pd.DataFrame):
    """Order gene sets by their most significant column and value.

    Assigns each gene set to the column in which it is most significant, then
    orders gene sets within each column by descending maximum value.

    Args:
        df (pandas.DataFrame): A gene set x program/sample matrix with
            -log10(p-value) values.

    Returns:
        pandas.DataFrame: The input matrix with rows reordered. Returned
        unchanged if it has no rows.
    """
    # sort gene sets by highest column and then highest value of that column
    if df.shape[0] > 0:
        stats = pd.DataFrame({"col": df.idxmax(axis=1), "max": df.max(axis=1)})
        ordered = []
        for col in df.columns:
            ordered.append(stats[stats["col"] == col].sort_values("max", ascending=False))
        ordered = pd.concat(ordered)
        ordered_df = df.loc[ordered.index]
    else:
        ordered_df = df
    return ordered_df


def plot_geneset_pval_heatmap(df: pd.DataFrame,
                              ax: Optional[Axes] = None,
                              axlegend: Optional[Axes] = None,
                              cmap: str = "Blues",
                              vmin: float = 0.,
                              vmax: float = 10.,
                              plot_title: str = "heatmap",
                              show_geneset_names: bool = False) -> Optional[Figure]:
    """Plot a heatmap of gene set -log10(p-value) scores.

    Draws a heatmap of ``df`` with a separate colorbar. When ``ax`` and
    ``axlegend`` are not supplied, new figures are created for the heatmap and
    the legend and returned.

    Args:
        df (pandas.DataFrame): Gene set x program matrix of -log10(p-value)
            values.
        ax (Optional[matplotlib.axes.Axes]): Axes on which to draw the heatmap.
            If None, a new figure and axes are created. Defaults to None.
        axlegend (Optional[matplotlib.axes.Axes]): Axes on which to draw the
            colorbar. If None, a new legend figure and axes are created.
            Defaults to None.
        cmap (str): Matplotlib colormap name. Defaults to ``"Blues"``.
        vmin (float): Minimum value of the color scale. Defaults to 0.0.
        vmax (float): Maximum value of the color scale. Defaults to 10.0.
        plot_title (str): Title for the heatmap. Defaults to ``"heatmap"``.
        show_geneset_names (bool): Whether to show gene set names as y-axis
            labels and adjust the figure size accordingly. Defaults to False.

    Returns:
        tuple[matplotlib.figure.Figure, matplotlib.figure.Figure]: The heatmap
        figure and the legend figure, returned only when both ``ax`` and
        ``axlegend`` are None; otherwise nothing is returned.
    """
    if ax is None:

        if show_geneset_names:
            figsize = [10 + df.shape[1]/4,  0.3 * df.shape[0]]
        else:
            figsize = [0.5 + df.shape[1]/4, 8]

        fig, ax_plot = plt.subplots(figsize=figsize, layout="constrained")
    else:
        ax_plot = ax
        
    if axlegend is None:
        figlegend, axlegend_plot = plt.subplots(figsize=[1, 3], layout="constrained")
    else:
        axlegend_plot = axlegend
    
    if df.shape[0] > 0:    
        sns.heatmap(df, cmap=cmap, vmin=vmin, vmax=vmax, ax=ax_plot, cbar_ax = axlegend_plot, yticklabels=show_geneset_names)
        ax_plot.tick_params(top=False,
                            bottom=True,
                            left=False,
                            right=False,
                            labelleft=False,
                            labelright=show_geneset_names,
                            labeltop= False,
                            labelbottom=True)

        ax_plot.set_xlabel("")
        ax_plot.set_ylabel("")
        ax_plot.set_title(plot_title)
        for _, spine in ax_plot.spines.items():
            spine.set_visible(True)
            spine.set_color('#aaaaaa')

    if ax is None and axlegend is None:
        return fig, figlegend
    

def plot_geneset_pval_clustermap(df: pd.DataFrame,
                            ax: Optional[Axes] = None,
                            axlegend: Optional[Axes] = None,
                            cmap: str = "Blues",
                            plot_title: str = "heatmap",
                            is_cluster: bool = True,
                            vmin: float = 0.,
                            vmax: float = 10.) -> Optional[Figure]:
    """Plot a heatmap and optional clustermap of gene set -log10(p-value) scores.

    Draws a heatmap of ``df`` with pathway labels (``source-name``) and, when
    ``is_cluster`` is True, an additional column-clustered clustermap.

    Args:
        df (pandas.DataFrame): Gene set x program matrix of -log10(p-value)
            values, indexed by pathway metadata including ``source`` and
            ``name``.
        ax (Optional[matplotlib.axes.Axes]): Axes on which to draw the heatmap.
            If None, a new figure and axes are created. Defaults to None.
        axlegend (Optional[matplotlib.axes.Axes]): Axes on which to draw the
            colorbar. If None, a new legend figure and axes are created.
            Defaults to None.
        cmap (str): Matplotlib colormap name. Defaults to ``"Blues"``.
        plot_title (str): Title for the plots. Defaults to ``"heatmap"``.
        is_cluster (bool): Whether to also produce a column-clustered
            clustermap. Defaults to True.
        vmin (float): Minimum value of the color scale. Defaults to 0.0.
        vmax (float): Maximum value of the color scale. Defaults to 10.0.

    Returns:
        tuple: The heatmap figure, the clustermap (a seaborn ``ClusterGrid``, or
        None when ``is_cluster`` is False), and the legend figure. The tuple is
        returned only when both ``ax`` and ``axlegend`` are None.
    """
    df_reset = df.reset_index()
    pathways_list = df_reset["source"] + "-" + df_reset["name"]
    figsize = [10 + df.shape[1]/4, 10 + df.shape[0]/4]

    if ax is None:
        fig, ax_plot = plt.subplots(figsize=figsize, layout="constrained")
    else:
        ax_plot = ax
        
    if axlegend is None:
        figlegend, axlegend_plot = plt.subplots(figsize=[1, 3], layout="constrained")
    else:
        axlegend_plot = axlegend

    if not df.empty:
        sns.heatmap(df, cmap=cmap, vmin=vmin, vmax=vmax, ax=ax_plot, cbar_ax=axlegend_plot, yticklabels=pathways_list)
        ax_plot.set(title=plot_title, xlabel="Program", ylabel="")
        ax_plot.tick_params(axis='y', labelrotation=0, labelright=True, labelleft=False)
        ax_plot.tick_params(axis='x', labelrotation=90)
        ax_plot.set_yticklabels([label.get_text()[:55] for label in ax_plot.get_yticklabels()])
        for _, spine in ax_plot.spines.items():
            spine.set_visible(True)
            spine.set_color('#aaaaaa')

    fig_cluster = None
    if is_cluster:
        fig_cluster = sns.clustermap(df, cmap=cmap, col_cluster=True, row_cluster=False, vmin=vmin, vmax=vmax, 
                                        figsize=figsize, yticklabels=pathways_list)
        fig_cluster.ax_heatmap.set(title=plot_title, xlabel="Program", ylabel="")
        fig_cluster.ax_heatmap.tick_params(axis='y', labelrotation=0, labelright=True, labelleft=False)
        fig_cluster.ax_heatmap.tick_params(axis='x', labelrotation=90)
        fig_cluster.ax_heatmap.set_yticklabels([label.get_text()[:55] for label in fig_cluster.ax_heatmap.get_yticklabels()])
        for _, spine in fig_cluster.ax_heatmap.spines.items():
            spine.set_visible(True)
            spine.set_color('#aaaaaa')
        fig_cluster.cax.set_visible(False)

    return fig, fig_cluster, figlegend if ax is None and axlegend is None else None
