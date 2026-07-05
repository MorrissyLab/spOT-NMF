import os
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
    """Compute pairwise statistics for program interactions based on usage counts.

    For each ordered pair of programs, counts the number of samples where
    program_two is present or absent among the samples where program_one is
    present, using a specified usage threshold. Pairwise statistics are computed
    in parallel using joblib for efficiency. If ``save_path`` is provided and a
    cached results file already exists, it is loaded and returned instead of
    being recomputed.

    Args:
        usage (pandas.DataFrame): DataFrame with programs as columns and samples
            as rows, containing usage counts for each program per sample.
        usage_threshold (Union[int, float]): Threshold for usage counts; a
            program is considered present in a sample if its count exceeds this
            value.
        sample (str): Condition name or sample identifier used to annotate
            results and name the output columns.
        save_path (str, optional): Directory to save (or load) the output file.
            If None, results are computed and returned without being saved.
            Defaults to None.
        file_prefix (str, optional): Prefix for the output file name. Defaults to
            None.

    Returns:
        pandas.DataFrame: DataFrame with pairwise statistics for each program
        pair, including the columns: ``group1`` (condition/sample identifier),
        ``n`` (number of samples considered), ``program_one`` (first program in
        the pair), ``program_two`` (second program in the pair),
        ``<sample>_P2pos`` (number of samples where program_two is present, i.e.
        > threshold, among those where program_one is present), and
        ``<sample>_P2neg`` (number of samples where program_two is absent, i.e.
        <= threshold, among those where program_one is present).
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
    """Generate node attributes for each program for use in network visualization.

    Calculates, for each program, the number and proportion of samples where
    usage exceeds the threshold, and whether the program is prevalent (at least
    ``prevalence_threshold`` samples above the usage threshold). Duplicate rows
    for the same sample and program are dropped.

    Args:
        usage (pandas.DataFrame): DataFrame with programs as columns and samples
            as rows, containing usage counts for each program per sample.
        usage_threshold (Union[int, float]): Threshold for usage counts; a
            program is considered present in a sample if its count exceeds this
            value.
        sample (str): Condition name or sample identifier used to annotate nodes.
        prevalence_threshold (int): Minimum number of samples above the usage
            threshold for a program to be considered prevalent. Defaults to 100.

    Returns:
        pandas.DataFrame: DataFrame with node attributes for each program,
        including the columns: ``sample_id`` (identifier for the
        condition/sample), ``total_spots`` (total number of samples/rows in the
        usage matrix), ``program`` (program name), ``num_samples_cps_gt_prevalence_threshold``
        (indicator, 1 or 0, of whether the number of samples with usage >
        threshold is at least ``prevalence_threshold``), ``cpp`` (number of
        samples where usage > threshold for the program), and ``proportion``
        (fraction of samples where usage > threshold for the program).
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


def build_network_graph(stats_df: pd.DataFrame, node_attrs: pd.DataFrame, sample: str, n_bins: int) -> Tuple[nx.DiGraph, nx.DiGraph]:
    """Build a directed NetworkX graph from pairwise statistics and node attributes.

    Adds a node for each program (carrying its prevalence attribute) and an edge
    between each program pair whose sample count meets the minimum bin threshold
    and that is not a self-loop. Edge weights are taken from the ``<sample>_val``
    enrichment column, with the absolute value stored as ``weight`` and the
    signed value as ``weight_col``.

    Args:
        stats_df (pandas.DataFrame): DataFrame containing pairwise statistics for
            program pairs, including a ``<sample>_val`` column of enrichment
            values.
        node_attrs (pandas.DataFrame): DataFrame containing node attributes for
            each program, including a ``num_samples_cps_gt_*`` prevalence column.
        sample (str): Condition name used to locate the ``<sample>_val`` edge
            weight column.
        n_bins (int): Minimum number of bins (samples) for an edge to be
            included.

    Returns:
        Tuple[networkx.DiGraph, networkx.DiGraph]: A tuple ``(graph,
        graph_filtered)`` where ``graph`` is the full directed graph with all
        nodes and edges, and ``graph_filtered`` is the subgraph containing only
        nodes with at least one edge.
    """

    graph = nx.DiGraph()

    # Add nodes with attributes
    prevalence_col = [col for col in node_attrs.columns if col.startswith('num_samples_cps_gt_')][0]
    for row in node_attrs.itertuples(index=False):
        graph.add_node(row.program, 
                       **{prevalence_col: getattr(row, prevalence_col)}, 
                       name=row.program)

    # Only keep edges that meet the threshold
    for row in stats_df.itertuples(index=False):
        if int(row.n) >= n_bins and row.program_one != row.program_two:
            w = getattr(row, f"{sample}_val")
            graph.add_edge(row.program_one, 
                           row.program_two, 
                           weight=abs(w), 
                           weight_col=w)
    
    # Compute the filtered graph by removing nodes with no edges
    graph_filtered = graph.subgraph([n for n in graph.nodes if graph.degree[n] > 0]).copy()

    return graph, graph_filtered


def detect_communities_infomap(graph: nx.DiGraph) -> nx.DiGraph:
    """Detect communities (clusters) in a graph using the Infomap algorithm via iGraph.

    Modifies the input graph in place by setting a ``cluster`` node attribute
    indicating community membership. Nodes that are isolated or otherwise absent
    from the clustering (including the case where the graph has no weighted
    edges) are assigned a cluster of ``-1``.

    Args:
        graph (networkx.DiGraph): The input directed graph. Its edges are
            expected to carry ``weight`` and ``weight_col`` attributes.

    Returns:
        networkx.DiGraph: The input graph with updated ``cluster`` node
        attributes.
    """
    # Convert NetworkX graph to DataFrame and extract weighted edge tuples
    edges_from_nx = nx.to_pandas_edgelist(graph)
    required_cols = {'source', 'target', 'weight', 'weight_col'}
    if edges_from_nx.empty or not required_cols.issubset(edges_from_nx.columns):
        # No weighted edges to cluster on; mark every node as unassigned (-1).
        for node in graph.nodes:
            graph.nodes[node]["cluster"] = -1
        return graph
    edge_tuples = list(edges_from_nx[['source', 'target', 'weight', 'weight_col']].itertuples(index=False, name=None))

    # Create igraph graph with node names as vertex IDs
    ig_graph = ig.Graph.TupleList(edge_tuples, directed=True, vertex_name_attr="name", edge_attrs=["weight", "weight_col"])

    # Run Infomap
    communities = ig_graph.community_infomap(trials=10, edge_weights="weight")

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
    """Compute node positions for visualization using an iGraph layout algorithm.

    Converts the graph to iGraph, ensures isolated nodes are included, and
    applies the requested layout algorithm to produce 2D coordinates.

    Args:
        graph (Union[networkx.Graph, networkx.DiGraph]): The input graph.
        layout_algorithm (str): iGraph layout algorithm to use (e.g., 'circle',
            'grid', 'star', 'graphopt'). Defaults to 'graphopt'.

    Returns:
        dict: Dictionary mapping node names to (x, y) coordinate tuples.
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
    """Draw a network graph with nodes colored by cluster and edges colored by weight.

    Adds legends for clusters and edge weights and adjusts node labels for
    clarity. If the graph has no nodes or edges, an empty figure is returned.

    Args:
        graph (Union[networkx.Graph, networkx.DiGraph]): The input graph.
        pos (dict): Dictionary mapping node names to (x, y) coordinates.
        title (Union[str, None]): Title for the plot axes. Defaults to None.
        suptitle (Union[str, None]): Super title for the figure. Defaults to None.

    Returns:
        Tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]: The figure and
        axes objects.
    """
    
    fig, ax = plt.subplots(figsize=(20, 15))

    # Check if the graph is empty
    if graph.number_of_nodes() == 0 or graph.number_of_edges() == 0:
        print("Empty graph: no nodes to plot.")
        return fig, ax

    # Edge attributes
    edge_weights = [graph[u][v].get('weight') for u, v in graph.edges()]
    edge_colors = [graph[u][v].get('weight_col') for u, v in graph.edges()]

    # Normalize edge weights for width and alpha
    edge_widths = np.array(edge_weights) + 2
    edge_alphas = np.array(edge_weights)

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
    # Compute quantiles from edge width data
    quantiles = np.quantile(edge_widths, [0.25, 0.5, 0.75, 1])
    labels = [f"{q:.2f}" for q in quantiles]

    # Create dummy lines with corresponding linewidths
    edge_legend_lines = [Line2D([0], [0], color='black', lw=w, label=label, alpha=0.8)
                         for w, label in zip(quantiles, labels)]

    # Add the legend to your plot
    edge_width_legend = ax.legend(handles=edge_legend_lines, title="Edge Weight", loc='lower left')
    ax.add_artist(edge_width_legend)
    
    return fig, ax


def plot_network_graph(graph: Union[nx.Graph, nx.DiGraph], pos: dict, sample: str, 
                          n_bins: int, edge_threshold: float, usage_threshold: Union[int, float],
                          save: bool = False, save_path: Union[str, None] = None, prefix: Union[str, None] = None):
    """Create network visualizations of the full network and cluster-specific subplots.

    Plots the whole network colored by cluster and, when more than one cluster
    is present, a second figure with one subplot per cluster. When ``save`` is
    True the figures are written to a PDF; otherwise they are returned. The
    ``edge_threshold``, ``n_bins``, and ``usage_threshold`` values are used to
    annotate the plot titles and build the output file name.

    Args:
        graph (Union[networkx.Graph, networkx.DiGraph]): The input graph.
        pos (dict): Dictionary mapping node names to (x, y) coordinates.
        sample (str): Sample name used in titles and file names.
        n_bins (int): Minimum number of bins for an edge to be included, used for
            annotation and file naming.
        edge_threshold (float): Threshold for edge weights, used for annotation.
        usage_threshold (Union[int, float]): Threshold for usage counts, used for
            file naming.
        save (bool): If True, save the plots to a PDF; if False, return the
            figures and axes. Defaults to False.
        save_path (Union[str, None]): Directory to save the plots. Required when
            ``save`` is True. Defaults to None.
        prefix (Union[str, None]): Prefix for the saved file name. Defaults to
            None.

    Returns:
        Tuple[matplotlib.figure.Figure, matplotlib.axes.Axes, matplotlib.figure.Figure, matplotlib.axes.Axes]:
        When ``save`` is False, returns ``(fig1, ax1, fig2, ax2)`` where the
        first pair is the whole-network plot and the second pair is the
        cluster-split plot (``fig2`` and ``ax2`` are None if fewer than two
        clusters exist). Returns None when ``save`` is True.

    Raises:
        ValueError: If ``save`` is True but ``save_path`` is None.
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


def plot_network_analysis(results_dir: str, sample_name: str, usage_threshold: Union[float, int], n_bins: int, edge_threshold: float, annot_file: Union[str, None]):
    """Run the full network analysis and visualization pipeline for a sample.

    Loads per-spot usage data from ``results_dir``, optionally renames program
    columns using an annotation file, thresholds usages at the 90th percentile,
    computes pairwise statistics, builds and clusters a network graph, and saves
    the network plots and the in-group/out-group connection heatmap. If no edges
    survive the thresholds, the plotting steps are skipped.

    Args:
        results_dir (str): Directory containing per-sample subfolders used to
            load and save results.
        sample_name (str): Name of the sample to analyze; also names the input
            subfolder and file.
        usage_threshold (Union[float, int]): Threshold for usage counts to
            consider a program present in a sample.
        n_bins (int): Minimum number of bins (samples) for an edge to be included
            in the network.
        edge_threshold (float): Threshold for edge weights used to filter the
            network.
        annot_file (Union[str, None]): Path to an annotation file with 'Program'
            and 'Annotation' columns used to rename program columns. If None, no
            annotations are applied.

    Raises:
        ValueError: If ``annot_file`` is provided but lacks the required
            'Program' and 'Annotation' columns.
    """
    
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

    # Filter usages based on the usage threshold
    for col in rf_usages.columns:
        thresh = rf_usages[col].quantile(0.90)
        rf_usages.loc[rf_usages[col] < thresh, col] = 0

    # Compute pairwise statistics
    stats_df = compute_pairwise_stats(usage=rf_usages, usage_threshold=usage_threshold, sample=sample_name, save_path=results_path)
    print("Stats DataFrame created and saved.")

    stats_df[f"{sample_name}_val"] = (stats_df[f"{sample_name}_P2pos"] + 1) / (stats_df[f"{sample_name}_P2pos"] + stats_df[f"{sample_name}_P2neg"] + 2)

    # Generate binary column for the edge threshold
    stats_df[f"gt{edge_threshold}"] = (stats_df[f"{sample_name}_val"].abs() > edge_threshold)

    # Subset the stats_df for the specified edge threshold
    stats_df = stats_df[stats_df[f"gt{edge_threshold}"] == 1]

    # Generate node attributes
    node_attrs = generate_node_attributes(rf_usages, usage_threshold, sample_name)

    # Build the network graph
    graph, graph_filtered = build_network_graph(stats_df, node_attrs, sample_name, n_bins)
    print("Network graph built.")

    # No edges survive the thresholds (common for small datasets) — nothing to plot.
    if graph.number_of_edges() == 0:
        print(
            f"No edges meet the thresholds (n_bins={n_bins}, edge_threshold={edge_threshold}) "
            f"for sample '{sample_name}'; skipping network plots. "
            f"Lower --n_bins / --edge_threshold or use a larger dataset."
        )
        return

    # Cluster the graph using Infomap and save the results
    graph = detect_communities_infomap(graph)
    pos = get_node_positions(graph, layout_algorithm='graphopt')

    plot_network_graph(graph=graph, pos=pos, sample=sample_name, 
                       n_bins=n_bins, edge_threshold=edge_threshold, 
                       usage_threshold=usage_threshold,
                       save=True, save_path=results_path, prefix=sample_name)

    # Cluster the filtered graph (no nodes lacking any edges) using Infomap and save the results
    graph_filtered = detect_communities_infomap(graph_filtered)
    pos_filtered = get_node_positions(graph_filtered, layout_algorithm='graphopt')

    plot_network_graph(graph=graph_filtered, pos=pos_filtered, sample=sample_name,
                       n_bins=n_bins, edge_threshold=edge_threshold, 
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


def calculate_outgoing_and_incoming_connections(graphobj):
    """Calculate in-group and out-group connection counts for each node in the graph.

    For every node, counts the edges (both incoming and outgoing) that connect
    it to other members of its own cluster (in-group) and to members of each
    other cluster (out-group), using the node ``cluster`` attribute.

    Args:
        graphobj (Union[networkx.Graph, networkx.DiGraph]): The input graph
            object whose nodes carry a ``cluster`` attribute.

    Returns:
        Tuple[pandas.DataFrame, pandas.DataFrame, pandas.DataFrame]: A tuple
        ``(group_connections, columnannot, rowannot)`` where ``group_connections``
        holds the per-node in-group and out-group connection counts (one row per
        program, columns named ``in_Cluster<id>`` / ``out_Cluster<id>``),
        ``columnannot`` maps each connection column to its connection type and
        cluster, and ``rowannot`` holds the per-node program and cluster
        annotations.
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
    """Plot a clustered heatmap of group connections with cluster annotations.

    Draws a seaborn clustermap of the connection values with optional row and
    column color bars keyed by cluster, a cluster legend, and a titled colorbar.
    When ``save`` is True and ``save_path`` is provided, the figure is written to
    a PDF; otherwise the figure and heatmap axes are returned.

    Args:
        group_connections (pandas.DataFrame): DataFrame of connection values to
            plot, with programs as the index and connection columns.
        rowannot (Union[pandas.DataFrame, None]): DataFrame with at least
            'program' and 'cluster' columns used for row annotations. If None, no
            row colors are drawn. Defaults to None.
        columnannot (Union[pandas.DataFrame, None]): DataFrame with at least 'V1'
            and 'cluster' columns used for column annotations. If None, no column
            colors are drawn. Defaults to None.
        custom_colors (list, optional): List of colors mapping cluster IDs to
            colors. If None, a built-in six-color palette is used. Defaults to
            None.
        figsize (tuple): Figure size. Defaults to (14, 10).
        cmap (str): Colormap for the heatmap. Defaults to "coolwarm".
        cluster_rows (bool): Whether to cluster (reorder) rows. Defaults to True.
        cluster_cols (bool): Whether to cluster (reorder) columns. Defaults to
            True.
        suptitle (str): Title for the whole figure. Defaults to "In-group and
            Out-group Connection Heatmap".
        legend_title (str): Title for the heatmap colorbar. Defaults to
            "n.edges".
        save (bool): Whether to save the figure to a PDF. Defaults to False.
        save_path (Union[str, None]): Folder in which to save the figure.
            Required for saving. Defaults to None.
        prefix (Union[str, None]): Prefix for the saved file name. Defaults to
            None.

    Returns:
        Tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]: The figure object
        and the heatmap axes. Returns None when the figure is saved instead of
        returned.
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