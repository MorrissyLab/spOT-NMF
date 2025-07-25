import os
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import networkx as nx
from infomap import Infomap
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from adjustText import adjust_text


def compute_pairwise_stats(usage, usage_threshold, cond1, save_path, file_prefix):
    stats_fn = os.path.join(save_path, f"uthresh{usage_threshold}_{file_prefix}_Cell_Cell_Interaction_Enrichment.tsv")
    
    # Check if stats file already exists
    if os.path.exists(stats_fn):
        return pd.read_csv(stats_fn, sep='\t')
    
    # If not exists, compute pairwise statistics
    else:
        programs = usage.columns
        cutoff = usage

        def pairwise(p1, p2):
            rows = cutoff.index[cutoff[p1] > usage_threshold]
            vals = cutoff.loc[rows, p2]
            pos = np.sum(vals > usage_threshold)
            neg = np.sum(vals <= usage_threshold)
            return {'group1': cond1, 
                    'n': pos + neg,
                    'program_one': p1, 
                    'program_two': p2,
                    f'{cond1}_P2pos': pos, 
                    f'{cond1}_P2neg': neg}

        # Parallel computation of pairwise statistics
        records = Parallel(n_jobs=-1)(delayed(pairwise)(p1, p2) for p1 in programs for p2 in programs)
        stats_df = pd.DataFrame.from_records(records)
        stats_df.to_csv(stats_fn, sep='\t', index=False)
        return stats_df

# TODO: This function needs to be tested with the Xenograft data and it won't work with multiple samples
def generate_node_attributes(usage, usage_threshold, cond1):
    total_spots = len(usage)
    programs = usage.columns

    # Calculate the number of samples with counts per sample (cps) greater than the threshold, the counts per program (cpp), and the proportion of samples with cps greater than the threshold
    node_summary = []
    for prog in programs:
        cps = (usage[prog] > usage_threshold).sum()
        node_summary.append({
            'sample_id': cond1,
            'total_spots': total_spots,
            'program': prog,
            'num_samples_cps_gt_100': int(cps >= 100),
            'cpp': cps,
            'proportion': cps / total_spots
        })

    node_summary = pd.DataFrame(node_summary)
    node_summary = (node_summary.drop_duplicates(subset=["sample_id", "program"]).loc[:, ["sample_id", "total_spots", "program", "num_samples_cps_gt_100", "cpp", "proportion"]])
    return node_summary


def build_network_graph(stats_df, node_attrs, cond1, n_bins):
    graph = nx.DiGraph()

    # Add nodes with attributes
    for row in node_attrs.itertuples(index=False):
        graph.add_node(row.program, 
                       num_samples_cps_gt_100=row.num_samples_cps_gt_100, 
                       name=row.program)

    # Only keep edges that meet the threshold
    for row in stats_df.itertuples(index=False):
        if row.n >= n_bins and row.program_one != row.program_two:
            w = getattr(row, f"{cond1}_val")
            graph.add_edge(row.program_one, 
                           row.program_two, 
                           weight=abs(w), 
                           weight_col=w)

    return graph


def detect_communities_infomap(graph):
    if not graph.edges():
        nx.set_node_attributes(graph, {node: 0 for node in graph.nodes()}, 'cluster')
        return graph
    
    ig_graph = ig.Graph.from_networkx(graph)
    ig_graph.vs["name"] = list(graph.nodes())
    communities = ig_graph.community_infomap()

    partition = communities.membership
    mapping = {node: cluster for node, cluster in zip(ig_graph.vs["name"], partition)}

    nx.set_node_attributes(graph, mapping, "cluster")
    return graph

    """Create a colormap similar to scico's hawaii palette"""
    colors = ['#3B0F70', '#8B0A50', '#D94824', '#F57C00', '#FFC107']
    return LinearSegmentedColormap.from_list('hawaii', colors[::-1], N=n)


def create_berlin_colormap():
    """Create a colormap similar to scico's berlin palette"""
    colors = ['#0571B0', '#92C5DE', '#F7F7F7', '#F4A582', '#CA0020']
    return LinearSegmentedColormap.from_list('berlin', colors[::-1])


def base_plot(graph, pos, sample):
    """
    Base plotting function for network visualization
    
    Parameters:
    - graph: NetworkX graph object
    - pos: dictionary of node positions
    - node_color_values: array/list of values for node coloring
    - sample: sample name string
    - cond1: condition 1 name
    """
    
    fig, ax = plt.subplots(1, 1, figsize=(20, 15))
    
    # Get edge attributes
    edge_weights = [graph[u][v].get('weight', 1) for u, v in graph.edges()]
    edge_colors = [graph[u][v].get('weight_col', 0) for u, v in graph.edges()]
    
    # Get node attributes
    node_sizes = [graph.nodes[node].get('num_samples_cps_gt_100', 100) for node in graph.nodes()]
    node_names = [graph.nodes[node].get('name', str(node)) for node in graph.nodes()]
    
    # Normalize edge weights for width and alpha
    if edge_weights:
        edge_widths = np.array(edge_weights)
        if max(edge_widths) > min(edge_widths):
            edge_widths = 0.3 + 2.7 * (edge_widths - min(edge_widths)) / (max(edge_widths) - min(edge_widths))
        else:
            edge_widths = np.ones_like(edge_widths) * 1.5
        
        edge_alphas = np.array(edge_weights)
        if max(edge_alphas) > min(edge_alphas):
            edge_alphas = 0.2 + 0.8 * (edge_alphas - min(edge_alphas)) / (max(edge_alphas) - min(edge_alphas))
        else:
            edge_alphas = np.ones_like(edge_alphas) * 0.6
    else:
        edge_widths = [1.0]
        edge_alphas = [0.6]
    
    # Draw edges with arrows
    edge_collection = None
    if graph.edges():
        edge_cmap = LinearSegmentedColormap.from_list('edge_cmap', ['black', 'magenta'])
        edge_norm = mcolors.Normalize(vmin=min(edge_colors), vmax=max(edge_colors))

        edge_collection = nx.draw_networkx_edges(graph, pos, ax=ax,
                                                 edge_color=edge_colors,
                                                 edge_cmap=edge_cmap,
                                                 edge_vmin=edge_norm.vmin,
                                                 edge_vmax=edge_norm.vmax,
                                                 width=edge_widths,
                                                 alpha=edge_alphas,
                                                 arrows=True,
                                                 arrowsize=20,
                                                 arrowstyle='->',
                                                 connectionstyle='arc3,rad=0.1')

        sm_edge = cm.ScalarMappable(norm=edge_norm, cmap=edge_cmap)
        sm_edge.set_array([])
        edge_cbar = plt.colorbar(sm_edge, ax=ax, shrink=0.6, pad=0.02)
        edge_cbar.set_label(f'{sample} [Edge Color]', fontsize=12)
    else:
        edge_collection = None
    
    # Draw nodes
    # node_cmap = create_berlin_colormap()  # You must define this or use an existing cmap
    # node_norm = mcolors.Normalize(vmin=min(node_color_values), vmax=max(node_color_values))

    node_collection = nx.draw_networkx_nodes(
        graph, pos, ax=ax,
        node_color="black",
        node_size=[max(s * 100, 50) for s in node_sizes],
        alpha=0.8
    )

    # sm_node = cm.ScalarMappable(norm=node_norm, cmap=node_cmap)
    # sm_node.set_array([])
    # node_cbar = plt.colorbar(sm_node, ax=ax, shrink=0.6, pad=0.08)
    # node_cbar.set_label('Node Color', fontsize=12)
    
    # Add node labels with repelling (approximate)
    texts = []
    for i, (node, (x, y)) in enumerate(pos.items()):
        if i < len(node_names):
            text = ax.text(x, y, node_names[i], 
                          fontsize=10, ha='center', va='center',
                          bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8, edgecolor='none'))
            texts.append(text)
    
    # Try to adjust text positions to avoid overlap
    try:
        adjust_text(texts, ax=ax)
    except:
        pass  # If adjust_text fails, keep original positions
    
    # Remove axes
    ax.set_axis_off()
    
    plt.tight_layout()
    return fig, ax


def plot_network_analysis(graph, graph_filtered, pos, pos_filtered, sample, n_bins, 
                          edge_threshold, cond1, save_path, usage_threshold, prefix):
    """
    Main plotting function that creates all network visualizations
    """
    
    # Create save directory if it doesn't exist
    samples_split_dir = os.path.join(save_path, "samples_split")
    os.makedirs(samples_split_dir, exist_ok=True)
    
    # Get node attributes for coloring
    fold_changes = [graph.nodes[node].get('fold_change', 1) for node in graph.nodes()]
    clusters = [graph.nodes[node].get('cluster', 0) for node in graph.nodes()]
    
    fold_changes_filtered = [graph_filtered.nodes[node].get('fold_change', 0) for node in graph_filtered.nodes()]
    clusters_filtered = [graph_filtered.nodes[node].get('cluster', 0) for node in graph_filtered.nodes()]
    
    # Create PDF file
    pdf_filename = os.path.join(samples_split_dir, 
                               f"usagethresh{usage_threshold}_{prefix}nbins{n_bins}_{sample}_Program_Program_Interaction_Network_Edge_Cutoff.pdf")
    
    with PdfPages(pdf_filename) as pdf:
        # 1. Plot the whole network (fold change coloring)
        fig1, ax1 = base_plot(graph, pos, fold_changes)
        
        # Set colormap for nodes (berlin palette, capped at -5, 5)
        if ax1.collections:
            node_collection = [c for c in ax1.collections if hasattr(c, 'set_cmap')][-1]
            if node_collection:
                node_collection.set_cmap(create_berlin_colormap())
                node_collection.set_clim(-5, 5)
        
        ax1.set_title(f"Program-Program Interaction Network: {sample}", 
                     fontsize=20, fontweight='bold', pad=20)
        
        subtitle = (f"Interactions with |{sample}| >{float(edge_threshold.replace('gt',''))/10} "
                   f"& Bins >= {n_bins} & P-Value <= 0.05\n{cond1} Prevalent Programs Are Positive")
        fig1.suptitle(subtitle, fontsize=14, y=0.02)
        
        pdf.savefig(fig1, bbox_inches='tight')
        plt.close(fig1)
        
        # 2. Plot without isolated nodes
        fig2, ax2 = base_plot(graph_filtered, pos_filtered, fold_changes_filtered)
        
        # Set colormap for nodes
        node_collection = ax2.collections[-1] if ax2.collections else None
        if node_collection:
            node_collection.set_cmap(create_berlin_colormap())
            node_collection.set_clim(-5, 5)
        
        ax2.set_title(f"Program-Program Interaction Network: {sample}", 
                     fontsize=20, fontweight='bold', pad=20)
        
        subtitle_filtered = (f"Interactions with |{sample}| >{float(edge_threshold.replace('gt',''))/10} "
                           f"& Bins >= {n_bins} & P-Value <= 0.05\n{cond1} Prevalent Programs Are Positive\n"
                           f"Isolated Nodes Removed")
        fig2.suptitle(subtitle_filtered, fontsize=14, y=0.02)
        
        pdf.savefig(fig2, bbox_inches='tight')
        plt.close(fig2)
        
        # 3. Plot the whole network (cluster coloring)
        fig3, ax3 = base_plot(graph, pos, clusters)
        
        # Set colormap for clusters
        node_collection = ax3.collections[-1] if ax3.collections else None
        if node_collection:
            node_collection.set_cmap(plt.cm.viridis)
        
        ax3.set_title("Network Colored by Cluster", fontsize=20, fontweight='bold', pad=20)
        
        cluster_subtitle = (f"Interactions with |{sample}| >{float(edge_threshold.replace('gt',''))/10} "
                          f"& Bins >= {n_bins} & P-Value <= 0.05\n{cond1} Prevalent Programs Are Positive")
        fig3.suptitle(cluster_subtitle, fontsize=14, y=0.02)
        
        pdf.savefig(fig3, bbox_inches='tight')
        plt.close(fig3)
        
        # 4. Plot without isolated nodes (cluster coloring)
        fig4, ax4 = base_plot(graph_filtered, pos_filtered, clusters_filtered)
        
        # Set colormap for clusters
        node_collection = ax4.collections[-1] if ax4.collections else None
        if node_collection:
            node_collection.set_cmap(plt.cm.viridis)
        
        ax4.set_title("Network Colored by Cluster", fontsize=20, fontweight='bold', pad=20)
        
        cluster_subtitle_filtered = (f"Interactions with |{sample}| >{float(edge_threshold.replace('gt',''))/10} "
                                   f"& Bins >= {n_bins} & P-Value <= 0.05\n{cond1} Prevalent Programs Are Positive\n"
                                   f"Isolated Nodes Removed")
        fig4.suptitle(cluster_subtitle_filtered, fontsize=14, y=0.02)
        
        pdf.savefig(fig4, bbox_inches='tight')
        plt.close(fig4)
        
        # 5. Additional cluster split plots (if multiple clusters exist)
        unique_clusters = set(clusters_filtered)
        if len(unique_clusters) > 1:
            # Create subplot for each cluster
            n_clusters = len(unique_clusters)
            fig5, axes = plt.subplots(1, n_clusters, figsize=(20*n_clusters, 15))
            if n_clusters == 1:
                axes = [axes]
            
            for i, cluster in enumerate(sorted(unique_clusters)):
                # Filter graph for this cluster
                cluster_nodes = [node for node in graph_filtered.nodes() 
                               if graph_filtered.nodes[node].get('cluster', 0) == cluster]
                cluster_subgraph = graph_filtered.subgraph(cluster_nodes)
                cluster_pos = {node: pos_filtered[node] for node in cluster_nodes if node in pos_filtered}
                cluster_fold_changes = [graph_filtered.nodes[node].get('fold_change', 0) for node in cluster_nodes]
                
                if cluster_subgraph.nodes():
                    # Draw this cluster
                    nx.draw_networkx_edges(cluster_subgraph, cluster_pos, ax=axes[i], alpha=0.6)
                    node_collection = nx.draw_networkx_nodes(
                        cluster_subgraph, cluster_pos, ax=axes[i],
                        node_color=cluster_fold_changes,
                        cmap=create_berlin_colormap(),
                        vmin=-5, vmax=5
                    )
                    
                    # Add labels
                    node_names = [graph_filtered.nodes[node].get('name', str(node)) for node in cluster_nodes]
                    for j, (node, (x, y)) in enumerate(cluster_pos.items()):
                        axes[i].text(x, y, node_names[j], fontsize=10, ha='center', va='center',
                                   bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
                
                axes[i].set_title(f"Cluster {cluster}", fontsize=16)
                axes[i].set_axis_off()
            
            fig5.suptitle("Network Split by Clusters", fontsize=20, fontweight='bold')
            pdf.savefig(fig5, bbox_inches='tight')
            plt.close(fig5)


def calculate_outgoing_and_incoming_connections(graphobj):
    """
    Calculate the number of outgoing and incoming connections for each node in a graph.
    """
    graph = graphobj

    # Get node and cluster info
    nodeDF = pd.DataFrame({
        "program": list(graph.nodes),
        "cluster": [graph.nodes[n].get("cluster", None) for n in G.nodes]
    })

    allclusters = sorted(nodeDF["cluster"].unique())

    edgeDF = pd.DataFrame(list(G.edges), columns=["X1", "X2"])

    in_group_connections = []
    out_group_connections = []

    for i, row in nodeDF.iterrows():
        ref_program = row["program"]
        node_cluster = row["cluster"]

        cluster_members = nodeDF[nodeDF["cluster"] == node_cluster]["program"].tolist()
        cluster_members = [m for m in cluster_members if m != ref_program]

        in_group_row = pd.Series(0, index=[f"in_Cluster{cl}" for cl in allclusters], name=ref_program)
        out_group_row = pd.Series(0, index=[f"out_Cluster{cl}" for cl in allclusters], name=ref_program)

        # Intra-cluster (in-group) connections
        outgoing = edgeDF[(edgeDF["X1"] == ref_program) & (edgeDF["X2"].isin(cluster_members))]
        incoming = edgeDF[(edgeDF["X2"] == ref_program) & (edgeDF["X1"].isin(cluster_members))]
        in_group_row[f"in_Cluster{node_cluster}"] = len(pd.concat([outgoing, incoming]).drop_duplicates())

        # Inter-cluster (out-group) connections
        other_clusters = [c for c in allclusters if c != node_cluster]
        for oc in other_clusters:
            other_members = nodeDF[nodeDF["cluster"] == oc]["program"].tolist()
            outgoing = edgeDF[(edgeDF["X1"] == ref_program) & (edgeDF["X2"].isin(other_members))]
            incoming = edgeDF[(edgeDF["X2"] == ref_program) & (edgeDF["X1"].isin(other_members))]
            out_group_row[f"out_Cluster{oc}"] = len(pd.concat([outgoing, incoming]).drop_duplicates())

        in_group_connections.append(in_group_row)
        out_group_connections.append(out_group_row)

    in_group_df = pd.DataFrame(in_group_connections)
    out_group_df = pd.DataFrame(out_group_connections)
    group_connections = pd.concat([in_group_df, out_group_df], axis=1)

    # Build column annotation
    colnames = group_connections.columns.tolist()
    columnannot = pd.DataFrame({"V1": colnames})
    columnannot[["connections", "cluster"]] = columnannot["V1"].str.extract(r"^(in|out)_Cluster(.+)$")

    # Row annotation (nodeDF)
    rowannot = nodeDF.set_index("program")

    # Matrix for heatmap
    conn_matrix = group_connections.fillna(0).astype(int)
    conn_matrix_log = np.log2(conn_matrix + 1)

    # Return everything (including matrix and placeholders for heatmaps)
    data_dict = {
        "connections_matrix": conn_matrix,
        "log_connections_matrix": conn_matrix_log,
        "column_annotations": columnannot,
        "row_annotations": rowannot
    }
    return data_dict

# def cc_interaction_networks(usage_matrix: pd.DataFrame, usage_threshold: float, 
#                             metadata: pd.DataFrame,
#                             cond1: str, cond2: str, n_bins: int, edge_threshold: str,
#                             connections_heatmap: bool, plot_pval_heatmap: bool,
#                             save_folder: str, file_prefix: str):
#     """
#     Create interaction networks based on usage matrix and metadata.
#     """
#     # Align usage to metadata
#     usage = usage_matrix.loc[metadata.index]
#     save_path = os.path.join(save_folder)
#     os.makedirs(save_path, exist_ok=True)

#     stats_fn = os.path.join(save_path, f"uthresh{usage_threshold}_{file_prefix}_Cell_Cell_Interaction_Enrichment.tsv")

#     if os.path.exists(stats_fn):
#         stats_df = pd.read_csv(stats_fn, sep='\t')
#     else:
#         programs = usage.columns
#         cutoff = usage

#         def pairwise(p1, p2):
#             rows = cutoff.index[cutoff[p1] > usage_threshold]
#             vals = cutoff.loc[rows, p2]
#             pos = np.sum(vals > usage_threshold)
#             neg = np.sum(vals <= usage_threshold)
#             # Using fisher_exact on contingency table [[pos, neg], [...] ]
#             # Here no "other group", just recording counts
#             return {'group1': cond1, 'n': pos + neg,
#                     'program_one': p1, 'program_two': p2,
#                     f'{cond1}_P2pos': pos, f'{cond1}_P2neg': neg}

#         records = Parallel(n_jobs=-1)(delayed(pairwise)(p1, p2) for p1 in programs for p2 in programs)
#         stats_df = pd.DataFrame.from_records(records)
#         stats_df.to_csv(stats_fn, sep='\t', index=False)

#     # Add pseudocount and optional threshold flags
#     def add_pc(pos, neg): return (pos + 1) / (pos + 1 + neg + 1)

#     stats_df[f"{cond1}_val"] = add_pc(stats_df[f"{cond1}_P2pos"], stats_df[f"{cond1}_P2neg"])

#     for i in range(1,6):
#         stats_df[f"gt0{i}"] = (stats_df[f"{cond1}_val"].abs() > (i/10)).astype(int)

#     # Node summary
#     cutoff = usage.copy()
#     cutoff['sample'] = cond1
#     total_spots = len(cutoff)
#     programs = [c for c in cutoff.columns if c.startswith('ot')]

#     node_summary = []
#     for prog in programs:
#         cps = (cutoff[prog] > usage_threshold).sum()
#         node_summary.append({'sample_id': cond1,
#                              'total_spots': total_spots,
#                              'program': prog,
#                              'num_samples_cps_gt_100': int(cps >= 100),
#                              'cpp': cps,
#                              'proportion': cps / total_spots})
#     node_summary = pd.DataFrame(node_summary)
#     node_attrs = node_summary.groupby('program').agg({'num_samples_cps_gt_100': 'first'}).reset_index()

#     # Build and plot network
#     G = nx.DiGraph()
    
#     # Add nodes
#     for row in node_attrs.itertuples(index=False):
#         G.add_node(row.program, 
#                    num_samples_cps_gt_100=row.num_samples_cps_gt_100,
#                    name=row.program)

#     # Add edges based on filters
#     for row in stats_df.itertuples(index=False):
#         threshold_met = getattr(row, edge_threshold, 0) == 1
#         if threshold_met and row.n >= n_bins and row.program_one != row.program_two:
#             weight = abs(getattr(row, f"{cond1}_val"))
#             G.add_edge(row.program_one, row.program_two, 
#                       weight=weight, weight_col=weight)

#     # Community detection with Infomap
#     if G.edges():
#         try:
#             # Create mapping from node names to integers
#             node2id = {node: i for i, node in enumerate(G.nodes())}
#             id2node = {i: node for node, i in node2id.items()}

#             # Initialize Infomap and add links
#             im = Infomap()
#             for u, v in G.edges():
#                 im.add_link(node2id[u], node2id[v])
            
#             im.run()

#             # Assign cluster labels
#             partition = {id2node[node.node_id]: node.module_id for node in im.nodes}
#             nx.set_node_attributes(G, partition, 'cluster')
#             print(f"Found {len(set(partition.values()))} clusters")
#         except Exception as e:
#             print(f"Infomap clustering failed: {e}. Using default clustering.")
#             nx.set_node_attributes(G, {node: 0 for node in G.nodes()}, 'cluster')
#     else:
#         print("No edges found, assigning all nodes to cluster 0")
#         nx.set_node_attributes(G, {node: 0 for node in G.nodes()}, 'cluster')

#     # Remove isolated nodes for filtered graph
#     Gf = G.subgraph([n for n, deg in G.degree() if deg > 0]).copy()

#     # Generate layouts
#     pos = nx.spring_layout(G, k=1, iterations=50) if G.nodes() else {}
#     pos_filtered = nx.spring_layout(Gf, k=1, iterations=50) if Gf.nodes() else {}

#     # Plot networks
#     plot_network_analysis(G, Gf, pos, pos_filtered, cond1, n_bins, 
#                           edge_threshold, cond1, save_path, usage_threshold, file_prefix)

#     # Connections heatmap
#     if connections_heatmap and G.edges():
#         inc = {}
#         out = {}
#         for u in G.nodes():
#             out[u] = G.out_degree(u, weight='weight')
#             inc[u] = G.in_degree(u, weight='weight')
        
#         conn_df = pd.DataFrame({
#             'program': list(G.nodes()), 
#             'incoming': [inc[u] for u in G.nodes()], 
#             'outgoing': [out[u] for u in G.nodes()]
#         })
        
#         conn_filename = os.path.join(save_path, "samples_split", f"connections_{cond1}.tsv")
#         conn_df.to_csv(conn_filename, sep='\t', index=False)
        
#         plt.figure(figsize=(12, 8))
#         sns.heatmap(conn_df.set_index('program'), cmap='viridis', annot=True, fmt='.2f')
#         plt.title("Incoming/Outgoing Connections")
#         plt.tight_layout()
#         plt.savefig(os.path.join(save_path, "samples_split", f"connections_{cond1}.png"), dpi=300)
#         plt.close()

#     # P-value heatmap (if p-values are available)
#     if plot_pval_heatmap and 'p' in stats_df.columns:
#         # Pivot data for heatmap
#         heat_df = stats_df.pivot(index='program_two', columns='program_one', values='p')

#         # Raw P-value Heatmap
#         plt.figure(figsize=(20, 16))
#         sns.heatmap(heat_df, cmap="Blues", square=True, cbar_kws={'label': 'p-value'},
#                     linewidths=0.5, linecolor='lightgrey', xticklabels=True, yticklabels=True)
#         plt.title(f'Fisher Exact P-Value: {cond1} vs {cond2}\n[Raw P-values]')
#         plt.xlabel('Program One')
#         plt.ylabel('Program Two')
#         plt.tight_layout()
#         plt.savefig(os.path.join(save_path, f"{file_prefix}_raw_pvalue_heatmap.png"), dpi=300)
#         plt.close()

#         # Log P-value heatmaps
#         logp_df = -np.log10(heat_df.replace(0, 1e-10))  # Replace 0 with small value
#         logp_df = logp_df.clip(upper=4)  # Cap at -log10(0.0001)

#         def plot_logp_heatmap(df, cluster_rows=False, cluster_cols=False, name_suffix=""):
#             plt.figure(figsize=(20, 16))
#             cg = sns.clustermap(df,
#                                 cmap='Reds',
#                                 row_cluster=cluster_rows,
#                                 col_cluster=cluster_cols,
#                                 linewidths=0.5,
#                                 linecolor='gray',
#                                 xticklabels=True,
#                                 yticklabels=True,
#                                 cbar_kws={'label': r'-log10(p-value) [Capped @ 4]'},
#                                 figsize=(20, 16))
#             plt.suptitle(f'-log10(P) Heatmap: {cond1} vs {cond2}', y=1.02)
#             file_name = f"{file_prefix}_logp_heatmap_{name_suffix}.png"
#             cg.savefig(os.path.join(save_path, file_name), dpi=300)
#             plt.close()

#         # Generate all 4 logP heatmaps
#         plot_logp_heatmap(logp_df, False, False, "no_clustering")
#         plot_logp_heatmap(logp_df, False, True, "col_clustering")
#         plot_logp_heatmap(logp_df, True, False, "row_clustering")
#         plot_logp_heatmap(logp_df, True, True, "full_clustering")

#         # Calculate connection analysis
#         if G.nodes():
#             # 1. Pivot data (like pivot_wider in R)
#             heat_df = stats_df.pivot(index='program_two', columns='program_one', values='p')

#             # 2. Matrix form for seaborn/heatmap use
#             heat_array = heat_df.to_numpy()
#             program_two_labels = heat_df.index.tolist()
#             program_one_labels = heat_df.columns.tolist()

#             # 3. Raw P-value Heatmap
#             plt.figure(figsize=(20, 16))
#             sns.heatmap(heat_df, cmap="Blues", square=True, cbar_kws={'label': 'p-value'},
#                         linewidths=0.5, linecolor='lightgrey', xticklabels=True, yticklabels=True)
#             plt.title(f'Fisher Exact P-Value: {cond1} vs {cond2}\n[Raw P-values]')
#             plt.xlabel('Program One')
#             plt.ylabel('Program Two')
#             plt.tight_layout()
#             plt.savefig(os.path.join(save_path, f"{prefix}_raw_pvalue_heatmap.png"), dpi=300)
#             plt.close()

#             # 4. Capped -log10(P) matrix
#             logp_df = -np.log10(heat_df)
#             logp_df.replace([np.inf, -np.inf], 0, inplace=True)
#             logp_df = logp_df.clip(upper=4)  # cap at -log10(0.0001)

#             # 5. Colormap for capped log10
#             reds = sns.color_palette("Reds", as_cmap=True)

#             # Variants: no clustering, column, row, both
#             def plot_logp_heatmap(df, cluster_rows=False, cluster_cols=False, name_suffix=""):
#                 sns.set(font_scale=0.8)
#                 cg = sns.clustermap(df,
#                                     cmap=reds,
#                                     row_cluster=cluster_rows,
#                                     col_cluster=cluster_cols,
#                                     linewidths=0.5,
#                                     linecolor='gray',
#                                     xticklabels=True,
#                                     yticklabels=True,
#                                     cbar_kws={'label': r'-log10(p-value) [Capped @ 4]'},
#                                     figsize=(20, 16))
#                 plt.suptitle(f'-log10(P) Heatmap: {cond1} vs {cond2}', y=1.02)
#                 file_name = f"{prefix}_logp_heatmap_{name_suffix}.png"
#                 cg.savefig(os.path.join(save_path, file_name), dpi=300)
#                 plt.close()

#             # Generate all 4 logP heatmaps
#             plot_logp_heatmap(logp_df, False, False, "no_clustering")
#             plot_logp_heatmap(logp_df, False, True, "col_clustering")
#             plot_logp_heatmap(logp_df, True, False, "row_clustering")
#             plot_logp_heatmap(logp_df, True, True, "full_clustering")