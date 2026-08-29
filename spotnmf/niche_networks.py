import os
import random
import numpy as np
import pandas as pd
from typing import Union, Tuple
from joblib import Parallel, delayed
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import matplotlib.patheffects as path_effects
import matplotlib.patches as mpatches
from adjustText import adjust_text
import igraph as ig
import seaborn as sns

def compute_pairwise_stats(usage: pd.DataFrame, usage_threshold: Union[int, float], sample: str, save_path=None, file_prefix=None) -> pd.DataFrame:
    """
    Computes pairwise statistics for program interactions based on usage counts.

    For each pair of programs, calculates the number of samples where program_two is present or absent among those where program_one is present, using a specified usage threshold.
    Pairwise statistics are computed in parallel using joblib for efficiency.

    Parameters
    ----------
    usage : pandas.DataFrame
        DataFrame with programs as columns and samples as rows, containing usage counts for each program per sample.
    usage_threshold : int or float
        Threshold for usage counts; a program is considered present in a sample if its count exceeds this value.
    sample : str
        Condition name or sample identifier to annotate results.
    save_path : str, optional
        Directory to save the output file. If None, results are not saved.
    file_prefix : str, optional
        Prefix for the output file name.

    Returns
    -------
    pandas.DataFrame
        DataFrame with pairwise statistics for each program pair, including:
        - group1: Condition/sample identifier.
        - n: Number of samples considered.
        - program_one: First program in the pair.
        - program_two: Second program in the pair.
        - <sample>_P2pos: Number of samples where program_two is present (> threshold) among those where program_one is present.
        - <sample>_P2neg: Number of samples where program_two is absent (<= threshold) among those where program_one is present.
        - Additional columns for enrichment values and thresholds if computed.
    """

    def _pairwise(p1, p2):
        # Find samples where program_one (p1) is present (usage > threshold)
        program_one_present = cutoff[p1] > usage_threshold
        rows_with_p1 = cutoff.index[program_one_present]
        
        # Get usage values for program_two (p2) in those samples
        program_two_values = cutoff.loc[rows_with_p1, p2]
        
        # Count samples where program_two is present or absent among those
        program_two_present = program_two_values > usage_threshold
        pos = program_two_present.sum()
        neg = (~program_two_present).sum()
        
        # Return statistics for this pair
        return {
            'group1': sample, 
            'n': pos + neg,
            'program_one': p1, 
            'program_two': p2,
            f'{sample}_P2pos': pos, 
            f'{sample}_P2neg': neg
        }

    if save_path is not None:
        stats_fn = os.path.join(save_path, f"uthresh{usage_threshold}_{file_prefix}_Cell_Cell_Interaction_Enrichment.tsv")
        
        # Check if stats file already exists
        if os.path.exists(stats_fn):
            return pd.read_csv(stats_fn, sep='\t')
    
        # If not exists, compute pairwise statistics
        else:
            programs = usage.columns
            cutoff = usage

            # Parallel computation of pairwise statistics
            records = Parallel(n_jobs=-1)(delayed(_pairwise)(p1, p2) for p1 in programs for p2 in programs)
            stats_df = pd.DataFrame.from_records(records)
            stats_df.to_csv(stats_fn, sep='\t', index=False)
            return stats_df
    else:
        # If save_path is None, compute pairwise statistics without saving
        programs = usage.columns
        cutoff = usage

        records = Parallel(n_jobs=-1)(delayed(_pairwise)(p1, p2) for p1 in programs for p2 in programs)
        stats_df = pd.DataFrame.from_records(records)
        return stats_df


def generate_node_attributes(usage: pd.DataFrame, usage_threshold: Union[int, float], sample: str, prevalence_threshold: int = 100):
    """
    Generates node attributes for each program for use in network visualization.

    Calculates, for each program, the number and proportion of samples where usage exceeds the threshold, and whether the program is prevalent (>=prevalence_threshold samples above threshold).

    Parameters
    ----------
    usage : pandas.DataFrame
        DataFrame with programs as columns and samples as rows, containing usage counts for each program per sample.
    usage_threshold : int or float
        Threshold for usage counts; a program is considered present in a sample if its count exceeds this value.
    sample : str
        Condition name or sample identifier to annotate nodes.
    prevalence_threshold : int, optional
        Minimum number of samples above threshold for a program to be considered prevalent (default: 100).

    Returns
    -------
    pandas.DataFrame
        DataFrame with node attributes for each program, including:
        - sample_id: Identifier for the condition/sample.
        - total_spots: Total number of samples (rows) in the usage matrix.
        - program: Program name (column in usage).
        - num_samples_cps_gt_<prevalence_threshold>: Indicator (1/0) if the number of samples with usage > threshold is at least prevalence_threshold (column name is dynamic based on the prevalence_threshold value).
        - cpp: Number of samples where usage > threshold for the program.
        - proportion: Fraction of samples where usage > threshold for the program.
    """

    total_spots = len(usage)
    programs = usage.columns

    # Calculate the number of samples with counts per sample (cps) greater than the threshold, the counts per program (cpp), and the proportion of samples with cps greater than the threshold
    node_summary = []
    for prog in programs:
        cps = (usage[prog] > usage_threshold).sum()
        node_summary.append({
            'sample_id': sample,
            'total_spots': total_spots,
            'program': prog,
            'num_samples_cps_gt_prevalence_threshold': int(cps >= prevalence_threshold),
            'cpp': cps,
            'proportion': cps / total_spots
        })

    node_summary = pd.DataFrame(node_summary)
    node_summary = node_summary.drop_duplicates(subset=["sample_id", "program"])
    return node_summary


def build_network_graph(stats_df: pd.DataFrame, node_attrs: pd.DataFrame, sample: str,
                        n_bins: int = 0, weight_column: Union[str, None] = None) -> Tuple[nx.DiGraph, nx.DiGraph]:
    """Build a directed NetworkX graph from pairwise statistics and node attributes.

    Adds a node for each program (carrying its prevalence attribute) and an edge
    for each non-self program pair that passes ``n_bins``. Each edge stores the
    signed enrichment as ``weight_col``, its magnitude as ``weight`` (community
    detection and layout require non-negative weights), and its direction of
    effect as ``sign``.

    The signed value is preferred from ``log2_oe`` when
    :func:`~spotnmf.niche_stats.cooccurrence_stats` produced the table, and
    falls back to the legacy ``<sample>_val`` conditional probability otherwise.
    This distinction matters: ``<sample>_val`` lies in (0, 1) so taking its
    magnitude was a no-op, but ``log2_oe`` is signed, and collapsing it with
    ``abs`` would draw a strongly *mutually exclusive* pair as though it were a
    strongly co-occurring one. ``sign`` keeps the two apart so
    :func:`detect_communities_infomap` can ignore depletion.

    Args:
        stats_df (pandas.DataFrame): Pairwise statistics for program pairs.
        node_attrs (pandas.DataFrame): Node attributes for each program,
            including a ``num_samples_cps_gt_*`` prevalence column.
        sample (str): Condition name used to locate the legacy
            ``<sample>_val`` edge weight column.
        n_bins (int): Minimum ``n`` for an edge to be included. Defaults to 0,
            since edges are normally already filtered by
            :func:`~spotnmf.niche_stats.select_edges`.
        weight_column (str | None): Explicit signed-weight column. Auto-detected
            when None. Defaults to None.

    Returns
    -------
    graph : networkx.DiGraph
        Full directed graph with all nodes and edges.
    graph_filtered : networkx.DiGraph
        Subgraph containing only nodes with at least one edge.
    """

    graph = nx.DiGraph()

    # Add nodes with attributes
    prevalence_col = [col for col in node_attrs.columns if col.startswith('num_samples_cps_gt_')][0]
    for row in node_attrs.itertuples(index=False):
        graph.add_node(row.program, 
                       **{prevalence_col: getattr(row, prevalence_col)}, 
                       name=row.program)

    if weight_column is None:
        weight_column = "log2_oe" if "log2_oe" in stats_df.columns else f"{sample}_val"

    # Only keep edges that meet the threshold
    for row in stats_df.itertuples(index=False):
        if int(row.n) >= n_bins and row.program_one != row.program_two:
            w = getattr(row, weight_column)
            graph.add_edge(row.program_one, 
                           row.program_two, 
                           weight=abs(w), 
                           weight_col=w,
                           sign=(1 if w >= 0 else -1))
    
    # Compute the filtered graph by removing nodes with no edges
    graph_filtered = graph.subgraph([n for n in graph.nodes if graph.degree[n] > 0]).copy()

    return graph, graph_filtered


def detect_communities_infomap(graph: nx.DiGraph, seed: int = 0) -> nx.DiGraph:
    """Detect communities (clusters) in a graph using the Infomap algorithm via iGraph.

    This function modifies the input graph in-place by setting a node 'cluster' attribute
    to indicate community membership.

    Infomap is stochastic, so it is run against a seeded random source: without
    a seed, two runs on identical input could return different niche labels,
    which makes the reported niches irreproducible. Edges marked as depleted
    (``sign == -1``) are excluded, because mutual exclusion is evidence against
    two programs sharing a niche, not for it.

    Args:
        graph (networkx.DiGraph): The input directed graph. Its edges are
            expected to carry ``weight`` and ``weight_col`` attributes.
        seed (int): Seed for the Infomap random source. Defaults to 0.

    Returns
    -------
    networkx.DiGraph
        The input graph with updated cluster attributes.
    """
    # Convert NetworkX graph to DataFrame and extract weighted edge tuples
    edges_from_nx = nx.to_pandas_edgelist(graph)
    required_cols = {'source', 'target', 'weight', 'weight_col'}
    if not edges_from_nx.empty and 'sign' in edges_from_nx.columns:
        edges_from_nx = edges_from_nx[edges_from_nx['sign'] > 0]
    if edges_from_nx.empty or not required_cols.issubset(edges_from_nx.columns):
        # No weighted edges to cluster on; mark every node as unassigned (-1).
        for node in graph.nodes:
            graph.nodes[node]["cluster"] = -1
        return graph
    edge_tuples = list(edges_from_nx[['source', 'target', 'weight', 'weight_col']].itertuples(index=False, name=None))

    # Create igraph graph with node names as vertex IDs
    ig_graph = ig.Graph.TupleList(edge_tuples, directed=True, vertex_name_attr="name", edge_attrs=["weight", "weight_col"])

    # Run Infomap against a seeded RNG so the partition is reproducible.
    rng_state = random.Random(seed)
    try:
        ig.set_random_number_generator(rng_state)
        communities = ig_graph.community_infomap(trials=10, edge_weights="weight")
    finally:
        ig.set_random_number_generator(random)

    # Map node name to community
    membership = communities.membership
    names = ig_graph.vs["name"]
    mapping = dict(zip(names, membership))

    # Set cluster memberships as node attributes in-place
    nx.set_node_attributes(graph, mapping, "cluster")

    # Assign default cluster (-1) to nodes missing from mapping (e.g., isolated nodes)
    for node in graph.nodes:
        if "cluster" not in graph.nodes[node]:
            graph.nodes[node]["cluster"] = -1
    
    return graph


def get_node_positions(graph: Union[nx.Graph, nx.DiGraph], layout_algorithm: str = 'graphopt') -> dict:
    """
    Computes node positions for visualization using a specified iGraph layout algorithm.

    Parameters
    ----------
    graph : networkx.Graph or networkx.DiGraph
        The input graph.
    layout_algorithm : str, optional
        iGraph layout algorithm to use (e.g., 'circle', 'grid', 'star', 'graphopt').

    Returns
    -------
    pos : dict
        dictionary mapping node names to (x, y) coordinates.
    """
    # Convert NetworkX graph to DataFrame and extract weighted edge tuples
    edges_from_nx = nx.to_pandas_edgelist(graph)
    edge_tuples = list(edges_from_nx[['source', 'target', 'weight']].itertuples(index=False, name=None))

    # Create igraph graph from edge list
    ig_graph = ig.Graph.TupleList(edge_tuples, directed=True, vertex_name_attr="name", edge_attrs=["weight"])

    # Ensure all nodes (including isolated) are present
    existing_nodes = set(ig_graph.vs["name"])
    all_nodes = set(graph.nodes)
    isolated_nodes = all_nodes - existing_nodes

    for node in isolated_nodes:
        ig_graph.add_vertex(name=node)

    # Now generate layout
    layout = ig_graph.layout(layout_algorithm)
    name_mapping = ig_graph.vs["name"]

    # Convert to networkx-style position dictionary
    pos = {name_mapping[i]: tuple(coord) for i, coord in enumerate(layout)}

    return pos


def base_plot(graph: Union[nx.Graph, nx.DiGraph], pos: dict, title: Union[str, None] = None, suptitle: Union[str, None] = None):
    """
    Plots a network graph with nodes colored by cluster and edges colored by weight.

    Adds legends for clusters and edge weights, and adjusts node labels for clarity.

    Parameters
    ----------
    graph : networkx.Graph or networkx.DiGraph
        The input graph.
    pos : dict
        dictionary mapping node names to (x, y) coordinates.
    title : str, optional
        Title for the plot.
    suptitle : str, optional
        Super title for the plot.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    ax : matplotlib.axes.Axes
        The axes object.
    """
    
    fig, ax = plt.subplots(figsize=(20, 15))

    # Check if the graph is empty
    if graph.number_of_nodes() == 0 or graph.number_of_edges() == 0:
        print("Empty graph: no nodes to plot.")
        return fig, ax

    # Edge attributes
    edge_weights = [graph[u][v].get('weight') for u, v in graph.edges()]
    edge_colors = [graph[u][v].get('weight_col') for u, v in graph.edges()]

    # Scale width and alpha from the weights' own range rather than using the
    # raw values. The legacy weights were conditional probabilities in (0, 1),
    # where feeding them straight to `alpha` happened to work; enrichment
    # weights such as log2(O/E) are unbounded, and would silently clip.
    w = np.asarray(edge_weights, dtype=float)
    w_span = np.ptp(w)
    w_unit = (w - w.min()) / w_span if w_span > 0 else np.full_like(w, 0.5)
    edge_widths = 1.0 + 4.0 * w_unit
    edge_alphas = 0.25 + 0.65 * w_unit

    # Draw edges
    edge_cmap = LinearSegmentedColormap.from_list('edge_cmap', ['black', 'magenta'])
    edge_norm = mcolors.Normalize(vmin=min(edge_colors), vmax=max(edge_colors))

    nx.draw_networkx_edges(graph, pos, ax=ax,
                           edge_color=edge_colors, edge_cmap=edge_cmap, edge_vmin=edge_norm.vmin, edge_vmax=edge_norm.vmax,
                           width=edge_widths, alpha=edge_alphas,
                           arrows=True, arrowsize=20, arrowstyle='->', connectionstyle='arc3,rad=0.1')
    sm_edge = cm.ScalarMappable(norm=edge_norm, cmap=edge_cmap)
    sm_edge.set_array([])
    cbar = plt.colorbar(sm_edge, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("[Edge Color]", fontsize=12)
    
    # Draw nodes
    clusters = nx.get_node_attributes(graph, "cluster")
    cluster_ids = sorted(set(clusters.values()))
    custom_colors = ['#28112f', '#2e2d5f', '#275d6d', '#218e66', '#063b69', '#d53f24']
    cluster_color_map = {cid: custom_colors[i % len(custom_colors)] for i, cid in enumerate(cluster_ids)}
    node_colors = [cluster_color_map[clusters[node]] for node in graph.nodes]
    node_sizes = 100

    nx.draw_networkx_nodes(graph, pos, ax=ax, node_color=node_colors, node_size=node_sizes, alpha=1)

    # Legend
    cluster_legend = ax.legend(handles=[plt.Line2D([0], [0], marker='o', color='w', 
                                                    markerfacecolor=cluster_color_map[cid], 
                                                    markersize=10, label=f"Cluster {cid}")
                                        for cid in cluster_ids],
                                title="Node Clusters", 
                                bbox_to_anchor=(1, 1), 
                                loc='upper left', 
                                frameon=True,
                                fancybox=True,
                                shadow=True)
    ax.add_artist(cluster_legend)

    # Node labels with adjust_text
    node_labels = [graph.nodes[node].get('name', str(node)) for node in graph.nodes()]
    texts = []
    for (x, y), label in zip(pos.values(), node_labels):
        t = ax.text(x, y, label, 
                    fontsize=10, ha='center', va='center', color='black',
                    path_effects=[path_effects.Stroke(linewidth=2.5, foreground='white'), path_effects.Normal()])
        texts.append(t)
    
    adjust_text(texts, ax=ax)

    ax.set_axis_off()
    if title:
        ax.set_title(title, fontsize=20, fontweight='bold', pad=20)
    if suptitle:
        fig.suptitle(suptitle, fontsize=14, y=0.02)

    # Edge Width and Color Legend
    # Label the legend with the underlying weights, not the rescaled widths.
    probs = [0.25, 0.5, 0.75, 1.0]
    labels = [f"{q:.2f}" for q in np.quantile(w, probs)]
    lw_quantiles = np.quantile(edge_widths, probs)

    # Create dummy lines with corresponding linewidths
    edge_legend_lines = [Line2D([0], [0], color='black', lw=lw, label=label, alpha=0.8)
                         for lw, label in zip(lw_quantiles, labels)]

    # Add the legend to your plot
    edge_width_legend = ax.legend(handles=edge_legend_lines, title="Edge Weight", loc='lower left')
    ax.add_artist(edge_width_legend)
    
    return fig, ax


def plot_network_graph(graph: Union[nx.Graph, nx.DiGraph], pos: dict, sample: str, 
                          n_bins: int, edge_threshold: float, usage_threshold: Union[int, float],
                          save: bool = False, save_path: Union[str, None] = None, prefix: Union[str, None] = None):
    """
    Creates and saves network visualizations, including the full network and cluster-specific subplots.

    Parameters
    ----------
    graph : networkx.Graph or networkx.DiGraph
        The input graph.
    pos : dict
        dictionary mapping node names to (x, y) coordinates.
    sample : str
        Sample name string.
    n_bins : int
        Minimum number of bins for an edge to be included.
    edge_threshold : float
        Threshold for edge weights.
    usage_threshold : int or float
        Threshold for usage counts.
    save : bool, optional
        If True, saves the plots to PDF. If False, returns the figures.
    save_path : str, optional
        Directory to save the plots.
    prefix : str, optional
        Prefix for the saved file names.

    Returns
    -------
    fig1, ax1 : matplotlib.figure.Figure, matplotlib.axes.Axes
        Figure and axes for the whole network plot.
    fig2, ax2 : matplotlib.figure.Figure, matplotlib.axes.Axes
        Figure and axes for the split cluster plots.
    """
    
    # Create save directory if it doesn't exist
    if save and save_path is None:
        raise ValueError("save_path must be provided if save is True.")
    
    # Get node attributes for coloring
    clusters = [graph.nodes[node].get('cluster', 0) for node in graph.nodes()]
    
    # 1. Plot the whole network (cluster coloring)
    fig1, ax1 = base_plot(graph, pos, 
                          title = "Network Colored by Cluster",
                          suptitle = (f"Interactions with |{sample}| > {edge_threshold}"
                                      f" & Bins >= {n_bins} & \n{sample} Prevalent Programs Are Positive"))
        
    # 2. Additional cluster split plots (if multiple clusters exist)
    unique_clusters = set(clusters)
    fig2, ax2 = None, None
    if len(unique_clusters) > 1:
        # Create subplot for each cluster
        n_clusters = len(unique_clusters)
        fig2, ax2 = plt.subplots(1, n_clusters, figsize=(20*n_clusters, 15))
        if n_clusters == 1:
            ax2 = [ax2]
        
        for i, cluster in enumerate(sorted(unique_clusters)):
            # Filter graph for this cluster
            cluster_nodes = [node for node in graph.nodes() 
                            if graph.nodes[node].get('cluster', 0) == cluster]
            cluster_subgraph = graph.subgraph(cluster_nodes)
            cluster_pos = {node: pos[node] for node in cluster_nodes if node in pos}
            
            if cluster_subgraph.nodes():
                # Draw this cluster
                nx.draw_networkx_edges(cluster_subgraph, cluster_pos, ax=ax2[i], alpha=0.6)
                node_collection = nx.draw_networkx_nodes(cluster_subgraph, cluster_pos, ax=ax2[i], node_color="black")
                
                # Add labels
                node_names = [graph.nodes[node].get('name', str(node)) for node in cluster_nodes]
                for j, (node, (x, y)) in enumerate(cluster_pos.items()):
                    ax2[i].text(x, y, node_names[j], fontsize=10, ha='center', va='center',
                                bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
            
            ax2[i].set_title(f"Cluster {cluster}", fontsize=16)
            ax2[i].set_axis_off()
        
        fig2.suptitle("Network Split by Clusters", fontsize=20, fontweight='bold')
    
    if save and save_path:
        samples_split_dir = os.path.join(save_path, "samples_split")
        os.makedirs(samples_split_dir, exist_ok=True)

        # Create PDF file
        pdf_filename = os.path.join(samples_split_dir, 
                                f"usagethresh{usage_threshold}_{prefix}_nbins{n_bins}_{sample}_Program_Program_Interaction_Network_Edge_Cutoff.pdf")
        with PdfPages(pdf_filename) as pdf:
            pdf.savefig(fig1, bbox_inches='tight')
            if fig2 is not None:
                pdf.savefig(fig2, bbox_inches='tight')
    else:
        return fig1, ax1, fig2, ax2


def plot_network_analysis(results_dir: str, sample_name: str,
                          usage_threshold: Union[float, int] = 0,
                          n_bins: int = 0,
                          edge_threshold: Union[float, None] = None,
                          annot_file: Union[str, None] = None,
                          presence_method: str = "otsu",
                          presence_quantile: float = 0.90,
                          fdr: float = 0.05,
                          min_log2_oe: float = 1.0,
                          min_prevalence_frac: float = 0.01,
                          null: str = "torus",
                          n_perm: int = 1000,
                          coords: Union[np.ndarray, None] = None,
                          seed: int = 0,
                          legacy: bool = False,
                          save_tables: bool = True):
    """Run the full niche inference and visualization pipeline for a sample.

    Loads per-spot usages, calls program presence per program, tests every
    program pair for co-occurrence enrichment against a spatial null, selects
    edges by false discovery rate and effect size, detects niches with Infomap,
    and writes the network plots and the in-group/out-group connection heatmap.

    This replaces the previous fixed-threshold rule (zero every program below
    its own 90th percentile, then keep pairs whose conditional co-occurrence
    probability exceeds 0.199). Those two constants were not independent: the
    per-column quantile forced every program to be present in exactly 10% of
    spots, so ``0.199`` silently meant "roughly two-fold enrichment over
    chance". On the CRC test matrix the edges that rule retained span
    ``log2_oe`` 1.00 to 2.67 -- i.e. it *was* a two-fold enrichment filter,
    just not written down as one, and its meaning would have changed silently
    had the quantile moved. The default ``min_log2_oe=1.0`` states the same
    two-fold target explicitly and holds it fixed regardless of how presence is
    called.

    Args:
        results_dir (str): Directory containing per-sample subfolders used to
            load and save results.
        sample_name (str): Name of the sample to analyze; also names the input
            subfolder and file.
        usage_threshold (Union[float, int]): Legacy usage threshold, used only
            when ``legacy`` is True. Defaults to 0.
        n_bins (int): Legacy absolute minimum spot count per edge, used only
            when ``legacy`` is True. Defaults to 0.
        edge_threshold (Union[float, None]): Legacy conditional-probability
            cutoff, used only when ``legacy`` is True. Defaults to None.
        annot_file (Union[str, None]): Path to an annotation file with
            'Program' and 'Annotation' columns used to rename program columns.
            If None, no annotations are applied.
        presence_method (str): How to call a program present in a spot --
            ``"otsu"`` (default), ``"mixture"`` or ``"quantile"``. See
            :func:`~spotnmf.niche_stats.call_presence`.
        presence_quantile (float): Quantile for ``presence_method="quantile"``.
            Defaults to 0.90, the previously hard-coded value.
        fdr (float): Benjamini-Hochberg false discovery rate for edge
            selection. Defaults to 0.05.
        min_log2_oe (float): Minimum log2 observed/expected co-occurrence.
            Defaults to 1.0 (two-fold).
        min_prevalence_frac (float): Minimum fraction of spots in which both
            programs must be present. Defaults to 0.01.
        null (str): Null model -- ``"torus"`` (default), ``"block"`` or
            ``"label"``. See :func:`~spotnmf.niche_stats.add_permutation_pvalues`.
        n_perm (int): Number of permutations. Defaults to 1000.
        coords (Union[numpy.ndarray, None]): Spatial coordinates aligned to the
            usage index. Parsed from Visium HD barcodes when None.
        seed (int): Random seed for permutations and Infomap. Defaults to 0.
        legacy (bool): Reproduce the original fixed-threshold behaviour exactly.
            Defaults to False.
        save_tables (bool): Write the per-program presence calls and the full
            pairwise statistics alongside the figures. Defaults to True.

    Raises:
        ValueError: If ``annot_file`` is provided but lacks the required
            'Program' and 'Annotation' columns, or if ``legacy`` is True and
            ``edge_threshold`` is None.
    """
    from spotnmf import niche_stats as nstats

    # Load the relevant files
    results_path = os.path.join(results_dir, sample_name)
    rf_usages = pd.read_csv(os.path.join(results_path, f"topics_per_spot_{sample_name}.csv"), index_col=0)

    # Load the annotation file if it exists
    if annot_file is not None:
        # Load and clean annotation
        annot = pd.read_csv(annot_file)

        if "Annotation" not in annot.columns or "Program" not in annot.columns:
            raise ValueError("Annotation file must contain 'Program' and 'Annotation' columns, in that order.")

        # Rename columns based on annotations
        annot_dict = dict(zip(annot["Program"], annot["Annotation"]))

        rf_usages.columns = [col.replace("ot_", "ot") + f"_{annot_dict[col]}" 
                            if col in annot_dict else col.replace("ot_", "ot")
                            for col in rf_usages.columns]

    if legacy:
        if edge_threshold is None:
            raise ValueError("edge_threshold is required when legacy=True.")
        stats_df, presence_info = _legacy_edges(
            rf_usages, sample_name, usage_threshold, n_bins, edge_threshold, results_path
        )
        weight_column = f"{sample_name}_val"
        graph_n_bins = n_bins
        # Legacy behaviour: node attributes came from the 90th-percentile-zeroed
        # usages, so reproduce that here rather than using the raw matrix.
        node_source = rf_usages.copy()
        for _col in node_source.columns:
            _thr = node_source[_col].quantile(0.90)
            node_source.loc[node_source[_col] < _thr, _col] = 0
        node_prevalence_threshold = 100
    else:
        if coords is None:
            coords = nstats.coords_from_index(rf_usages.index)
            if coords is None:
                print(
                    "Could not parse spatial coordinates from the usage index; falling back to "
                    "label permutation. This ignores spatial autocorrelation and is "
                    "anticonservative -- pass coords=adata.obsm['spatial'] for a valid test."
                )

        presence, presence_info = nstats.call_presence(
            rf_usages, method=presence_method, quantile=presence_quantile, seed=seed
        )
        print("Realised prevalence per program: "
              f"min={presence_info['prevalence'].min():.4f} "
              f"median={presence_info['prevalence'].median():.4f} "
              f"max={presence_info['prevalence'].max():.4f}")

        all_stats = nstats.cooccurrence_stats(presence, sample=sample_name)
        all_stats = nstats.add_permutation_pvalues(
            all_stats, presence, coords=coords, null=null, n_perm=n_perm, seed=seed
        )
        stats_df = nstats.select_edges(
            all_stats, fdr=fdr, min_log2_oe=min_log2_oe,
            min_prevalence_frac=min_prevalence_frac,
        )
        # `build_network_graph` still reads `<sample>_val`; point it at the
        # enrichment statistic so the legacy column name keeps working.
        stats_df[f"{sample_name}_val"] = stats_df["log2_oe"]
        weight_column = "log2_oe"
        graph_n_bins = 0
        # Node prevalence must come from the presence calls. Consensus usages
        # have no exact zeros, so counting `raw_usage > 0` would report every
        # program as present in every spot.
        node_source = presence.astype(int)
        node_prevalence_threshold = max(int(min_prevalence_frac * len(rf_usages)), 1)

        print(f"Selected {len(stats_df)} of {len(all_stats)} ordered pairs "
              f"(null={null}, n_perm={n_perm}, FDR<{fdr}, log2(O/E)>={min_log2_oe}).")

        if save_tables:
            presence_info.to_csv(os.path.join(results_path, f"{sample_name}_presence_calls.csv"))
            all_stats.to_csv(
                os.path.join(results_path, f"{sample_name}_cooccurrence_stats.csv"), index=False
            )

    # Generate node attributes
    node_attrs = generate_node_attributes(
        node_source, usage_threshold, sample_name,
        prevalence_threshold=node_prevalence_threshold,
    )

    # Build the network graph
    graph, graph_filtered = build_network_graph(
        stats_df, node_attrs, sample_name, graph_n_bins, weight_column=weight_column
    )
    print("Network graph built.")

    # No edges survive the thresholds (common for small datasets) -- nothing to plot.
    if graph.number_of_edges() == 0:
        print(
            f"No edges survived selection for sample '{sample_name}'. "
            "Relax --fdr or --min_log2_oe, raise --n_perm, or check the realised "
            "prevalences printed above."
        )
        return

    label = f"FDR<{fdr}, log2(O/E)>={min_log2_oe}" if not legacy else f"|val| > {edge_threshold}"

    # Cluster the graph using Infomap and save the results
    graph = detect_communities_infomap(graph, seed=seed)
    pos = get_node_positions(graph, layout_algorithm='graphopt')

    plot_network_graph(graph=graph, pos=pos, sample=sample_name, 
                       n_bins=graph_n_bins, edge_threshold=label, 
                       usage_threshold=usage_threshold,
                       save=True, save_path=results_path, prefix=sample_name)

    # Cluster the filtered graph (no nodes lacking any edges) using Infomap and save the results
    graph_filtered = detect_communities_infomap(graph_filtered, seed=seed)
    pos_filtered = get_node_positions(graph_filtered, layout_algorithm='graphopt')

    plot_network_graph(graph=graph_filtered, pos=pos_filtered, sample=sample_name,
                       n_bins=graph_n_bins, edge_threshold=label, 
                       usage_threshold=usage_threshold,
                       save=True, save_path=results_path, prefix=str(sample_name + "_filtered"))
    print("Network analysis plots saved.")
    
    # Plot the inside-outside connections heatmap
    group_connections, columnannot, rowannot = calculate_outgoing_and_incoming_connections(graph)
    group_connections = np.log2(group_connections + 1)

    plot_connection_heatmap(group_connections,
                            rowannot=rowannot, columnannot=columnannot, 
                            figsize=(12, 10), cmap="RdBu_r",
                            cluster_rows=True, cluster_cols=True,
                            suptitle="In-group and Out-group Connections Heatmap", legend_title="Log2(n.edges + 1)",
                            save=True, save_path=results_path, prefix=sample_name)
    print("Network connections heatmap saved.")


def _legacy_edges(rf_usages: pd.DataFrame, sample_name: str, usage_threshold, n_bins: int,
                  edge_threshold: float, results_path: str):
    """Reproduce the original fixed-threshold edge selection exactly.

    Retained so previously published networks stay reproducible and can be
    included as one cell of :func:`~spotnmf.niche_stats.parameter_sweep`. Not
    recommended for new analyses: it applies no null model and no
    multiple-testing correction across the P(P-1) program pairs tested.

    Args:
        rf_usages (pandas.DataFrame): Spots-by-programs usage matrix.
        sample_name (str): Sample identifier.
        usage_threshold (Union[float, int]): Usage above which a program counts
            as present, applied after the 90th-percentile zeroing.
        n_bins (int): Minimum absolute spot count for an edge.
        edge_threshold (float): Conditional-probability cutoff.
        results_path (str): Directory used to cache pairwise statistics.

    Returns:
        Tuple[pandas.DataFrame, pandas.DataFrame]: The retained edges and a
        presence-call summary.
    """
    usages = rf_usages.copy()
    for col in usages.columns:
        thresh = usages[col].quantile(0.90)
        usages.loc[usages[col] < thresh, col] = 0

    stats_df = compute_pairwise_stats(
        usage=usages, usage_threshold=usage_threshold, sample=sample_name,
        save_path=results_path, file_prefix=sample_name
    )
    stats_df[f"{sample_name}_val"] = (
        (stats_df[f"{sample_name}_P2pos"] + 1)
        / (stats_df[f"{sample_name}_P2pos"] + stats_df[f"{sample_name}_P2neg"] + 2)
    )
    stats_df = stats_df[stats_df[f"{sample_name}_val"].abs() > edge_threshold]

    presence = usages > usage_threshold
    info = pd.DataFrame({"prevalence": presence.mean(axis=0)})
    info.index.name = "program"
    return stats_df, info


def calculate_outgoing_and_incoming_connections(graphobj):
    """
    Calculate in-group and out-group connections for each node in the graph.

    Parameters
    ----------
    graphobj : nx.Graph or nx.DiGraph
        The input graph object containing nodes and edges.
    
    Returns
    -------
    group_connections : pd.DataFrame
        DataFrame with in-group and out-group connections for each node.
    columnannot : pd.DataFrame
        DataFrame with column annotations for the connections.
    rowannot : pd.DataFrame
        DataFrame with row annotations for the nodes.
    """
    graph = graphobj

    # Node and cluster info
    node_data = [{'program': n, 'cluster': graph.nodes[n].get('cluster')} for n in graph.nodes()]
    nodeDF = pd.DataFrame(node_data)
    allclusters = nodeDF['cluster'].unique()

    # Edge info
    edgeDF = pd.DataFrame(list(graph.edges()), columns=['X1', 'X2'])

    in_group_connections = pd.DataFrame()
    out_group_connections = pd.DataFrame()

    for i, row in nodeDF.iterrows():
        ref_program = row['program']
        node_cluster = row['cluster']

        cluster_members = nodeDF[nodeDF['cluster'] == node_cluster]['program'].tolist()
        cluster_members = [m for m in cluster_members if m != ref_program]

        in_group_tempDF = pd.DataFrame(0, index=[ref_program], columns=[f"in_Cluster{cl}" for cl in allclusters])
        out_group_tempDF = pd.DataFrame(0, index=[ref_program], columns=[f"out_Cluster{cl}" for cl in allclusters])

        # In-group edges
        edgeDF_in_outgoing = edgeDF[(edgeDF['X1'] == ref_program) & (edgeDF['X2'].isin(cluster_members))]
        edgeDF_in_incoming = edgeDF[(edgeDF['X2'] == ref_program) & (edgeDF['X1'].isin(cluster_members))]
        edgeDF_in = pd.concat([edgeDF_in_outgoing, edgeDF_in_incoming]).drop_duplicates()

        in_group_tempDF[f"in_Cluster{node_cluster}"] = len(edgeDF_in)

        # Out-group edges
        other_clusters = [cl for cl in allclusters if cl != node_cluster]
        for other_cluster in other_clusters:
            other_cluster_members = nodeDF[nodeDF['cluster'] == other_cluster]['program'].tolist()

            edgeDF_out_outgoing = edgeDF[(edgeDF['X1'] == ref_program) & (edgeDF['X2'].isin(other_cluster_members))]
            edgeDF_out_incoming = edgeDF[(edgeDF['X2'] == ref_program) & (edgeDF['X1'].isin(other_cluster_members))]
            edgeDF_out = pd.concat([edgeDF_out_outgoing, edgeDF_out_incoming]).drop_duplicates()

            out_group_tempDF[f"out_Cluster{other_cluster}"] = len(edgeDF_out)

        in_group_connections = pd.concat([in_group_connections, in_group_tempDF])
        out_group_connections = pd.concat([out_group_connections, out_group_tempDF])

    group_connections = pd.concat([in_group_connections, out_group_connections], axis=1)

    columnannot = pd.DataFrame({'V1': group_connections.columns})
    columnannot[['connections', 'cluster']] = columnannot['V1'].str.split('_Cluster', expand=True)

    rowannot = nodeDF.copy()

    return group_connections, columnannot, rowannot


def plot_connection_heatmap(group_connections, rowannot: Union[pd.DataFrame, None] = None, columnannot: Union[pd.DataFrame, None] = None,
                            custom_colors=None, figsize=(14, 10), cmap="coolwarm",
                            cluster_rows: bool = True, cluster_cols: bool = True,
                            suptitle: str = "In-group and Out-group Connection Heatmap", legend_title: str = "n.edges",
                            save: bool = False, save_path: Union[str, None] = None, prefix: Union[str, None] = None) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot a clustered heatmap of group connections with optional row and column cluster annotations and legends.

    Parameters
    ----------
    group_connections : pd.DataFrame
        Square DataFrame of connection values to plot.
    rowannot : pd.DataFrame, optional
        DataFrame with at least ['program', 'cluster'] columns for row annotations.
    columnannot : pd.DataFrame, optional
        DataFrame with at least ['V1', 'cluster'] columns for column annotations.
    custom_colors : list of str, optional
        List of colors to map cluster IDs to colors.
    figsize : tuple, optional
        Figure size.
    cmap : str, optional
        Colormap for heatmap.
    cluster_rows : bool, optional
        Whether to cluster rows.
    cluster_cols : bool, optional
        Whether to cluster columns.
    suptitle : str, optional
        Title for the whole figure.
    legend_title : str, optional
        Title for the heatmap colorbar.
    save : bool, optional
        Whether to save the figure.
    save_path : str or Path, optional
        Folder to save figure.
    prefix : str, optional
        Subfolder or prefix for save_path.

    Returns
    -------
    fig : plt.Figure
        Figure object.
    ax_heatmap : plt.Axes
        Heatmap axes.
    """
    # Define color map for clusters (default: 0 to 5)
    if custom_colors is None:
        custom_colors = ['#28112f', '#2e2d5f', '#275d6d', '#218e66', '#063b69', '#d53f24']
    cluster_ids = list(range(len(custom_colors)))
    cluster_color_map = dict(zip(cluster_ids, custom_colors))

    # Row colors
    if rowannot is not None:
        row_clusters = rowannot.set_index('program').reindex(group_connections.index)['cluster'].astype(int)
        row_colors = row_clusters.map(cluster_color_map)
    else:
        row_colors = None
        row_clusters = None

    # Column colors
    if columnannot is not None:
        col_clusters = columnannot.set_index('V1').reindex(group_connections.columns)['cluster'].astype(int)
        col_colors = col_clusters.map(cluster_color_map)
    else:
        col_colors = None
        col_clusters = None

    # Create the clustermap
    g = sns.clustermap(
        group_connections,
        figsize=figsize,
        cmap=cmap,
        row_cluster=cluster_rows,
        col_cluster=cluster_cols,
        row_colors=row_colors,
        col_colors=col_colors,
        xticklabels=True,
        yticklabels=True,
        dendrogram_ratio=(0.1, 0.1),  # optional: reduce dendrogram size
        cbar_pos=(1, 0.3, 0.05, 0.2)  # optional: move colorbar
    )

    # Function to create a legend
    def create_cluster_legend(cluster_labels, color_map):
        unique_clusters = sorted(set(cluster_labels))
        handles = [
            mpatches.Patch(color=color_map[cl], label=f"Cluster {cl}")
            for cl in unique_clusters if cl in color_map
        ]
        return handles

    # Create legends
    row_legend_handles = create_cluster_legend(row_clusters, cluster_color_map) if row_clusters is not None else []

    # Combine both legends
    all_handles = row_legend_handles

    # Draw combined legend to the right of the figure
    if all_handles:
        g.fig.legend(handles=all_handles,
                     title="Clusters",
                     loc='center left',
                     bbox_to_anchor=(1, 0.7), frameon=False)

    g.fig.suptitle(suptitle, fontsize=14, y=1.05)

    if g.cax is not None:
        g.cax.set_title(legend_title, fontsize=12)

    if save and save_path is not None:
        if prefix is not None:
            filename = f"{prefix}_inside_outside_connection_heatmap.pdf"
        else:
            filename = "inside_outside_connection_heatmap.pdf"
        g.savefig(os.path.join(save_path, filename), bbox_inches='tight')
    else:
        return g.fig, g.ax_heatmap