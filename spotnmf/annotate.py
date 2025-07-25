import os
import pandas as pd
from .io import load_experiment_result, check_dir, get_ground_truth, save_ranked_genes, read_gmt_file
from sklearn.metrics.pairwise import cosine_similarity
import seaborn as sns
import numpy as np
import rbo
from .pl import plot_df_heatmap
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import random
from tqdm import tqdm
from matplotlib.patches import Patch

def list_genesets(genome = None):
    parent_directory = os.path.dirname(os.path.dirname(__file__))
    geneset_dir = os.path.join(parent_directory, 'data', 'genesets')
    available_genesets =  [f.replace(".gmt", "") for f in os.listdir(geneset_dir) if f.endswith(".gmt")]

    if(genome == "mm10"):
        available_genesets = [f for f in available_genesets if f.startswith("mouse")]
    elif(genome == "GRCh38"):
        available_genesets = [f for f in available_genesets if f.startswith("human")]

    return available_genesets
    

def compute_genesets_annotation(rf_usages, gene_set, results_dir_path, max_top_genes=100, ranking_method="rboext", experiment_title="experiment"):
    """
    Compute gene set annotations by calculating similarity scores between 
    reference gene sets and ranked gene usage lists.
    
    Parameters:
    -----------
    rf_usages : pd.DataFrame
        DataFrame with columns as different programs and rows as genes, with values indicating gene usage scores.
    gene_set : str
        Name of the gene set file (without extension) to be used for annotation.
    max_top_genes : int, optional, default=100
        Maximum number of top genes to consider from each program.
    ranking_method : str, optional, default="rboext"
        Method to calculate similarity scores. Options: "rboext" (extrapolated Rank-Biased Overlap),
        "rbo" (rank-biased overlap), "mgs" (marker gene scores).
    experiment_title : str, optional, default="experiment"
        Title for the experiment, used in output filenames.
    
    Returns:
    --------
    None
        Saves gene score results and plots heatmaps in specified output directory.
    """

    # Define paths and directories

    parent_directory = os.path.dirname(os.path.dirname(__file__))
    geneset_dir = os.path.join(parent_directory, 'data', 'genesets')
    geneset_file = os.path.join(geneset_dir, gene_set + ".gmt")

    # Check if geneset file exists, return available gene sets if not found
    if not os.path.isfile(geneset_file):
        available_genesets = list_genesets()
        print(f"Gene set '{gene_set}' not found. Available gene sets: {available_genesets}")
        return

    # Read gene set data
    df_geneset = read_gmt_file(geneset_file)
    
    # Prepare output directory
    output_path = check_dir(os.path.join(results_dir_path, "genesets_results"))
    

    df_gene_scores = pd.DataFrame(index=df_geneset.columns.values, columns=rf_usages.columns.values)
    for q_program in df_geneset.columns: 
        for p_program in rf_usages.columns:
            query_list = df_geneset[q_program].dropna().unique().tolist()
            n_top_genes = min(len(query_list)*2, max_top_genes)

            program_list = rf_usages[p_program].sort_values(ascending=False).head(n_top_genes).index.to_list()
            
            results_score = get_ranking_score(query_list, program_list, rank_type = ranking_method)
            df_gene_scores.at[q_program, p_program] = results_score


    df_gene_scores.to_csv(os.path.join(output_path, f"genesets_scores_{experiment_title}.csv"))
    save_ranked_genes(df_gene_scores, results_file = os.path.join(output_path, f"ranked_genesets_{experiment_title}.csv"), n_top_genes=5)
    plot_df_heatmap(df_gene_scores,title_name =experiment_title, x_name = "Topic", y_name = "Geneset", results_dir_path= output_path,is_cluster=True)   


def get_ranking_score(query_list, program_list, rank_type = 'rbo'):
    intersected_genes = [(x,i+1) for i, x in enumerate(program_list) if x in query_list]
    if(rank_type == 'rbo'):
        score = rbo.RankingSimilarity(query_list, program_list).rbo()
    elif(rank_type == "rboext"):
        score = rbo.RankingSimilarity(query_list, program_list).rbo_ext()
    elif(rank_type == 'mgs'):
        score = sum(1 / value_tuple[1] for value_tuple in intersected_genes if value_tuple[1] != 0)
    
    return score


def sum_cell_types(adata_sc, cluster_key="celltype"):
    cell_annotations = adata_sc.obs[[cluster_key]]
    cell_annotations.index.name = "spots"

    counts_df = adata_sc.to_df()
    counts_df.columns = adata_sc.var.index
    counts_df.index = adata_sc.obs.index
    counts_df.index.name = "spots"

    cell_type_gene_df  = pd.merge(counts_df, adata_sc.obs[cluster_key], on='spots').groupby(cluster_key).sum().T
    
    return cell_type_gene_df

def mean_cell_types(adata_sc, cluster_key="celltype"):
    cell_annotations = adata_sc.obs[[cluster_key]]
    cell_annotations.index.name = "spots"

    counts_df = adata_sc.to_df()
    counts_df.columns = adata_sc.var.index
    counts_df.index = adata_sc.obs.index
    counts_df.index.name = "spots"

    # Merge counts with annotations
    merged_df = pd.merge(counts_df, adata_sc.obs[cluster_key], on='spots')

    # Group by cluster and compute mean
    cell_type_gene_df = merged_df.groupby(cluster_key).mean().T

    return cell_type_gene_df
    

def get_annotation_from_corr(df_corr):
    # Initialize list for annotation
    annotation_df = []
    # Flatten the DataFrame and sort by correlation values in descending order
    sorted_correlations = df_corr.stack().sort_values(ascending=False)

    # Initialize a set to track used cell types
    used_celltypes = set()
    used_program = set()
    # Iterate through the sorted correlations
    for (program, celltype), correlation in sorted_correlations.items():
        if celltype not in used_celltypes and program not in used_program:
            annotation_df.append([program, celltype, correlation])
            used_celltypes.add(celltype)
            used_program.add(program)  # Mark this celltype as used

    # Convert annotation to DataFrame
    annotation_df = pd.DataFrame(annotation_df, columns=["program", "celltype", "score"])
    return annotation_df


def annotate_programs_by_ground_truth(genes_topics_df, ground_truth_cell_type_gene_df, correlation_type='pearson', top_n_features=500):
    # Create empty DataFrame for correlation or scores
    df_corr = pd.DataFrame(index=genes_topics_df.columns, columns=ground_truth_cell_type_gene_df.columns)
    
    if correlation_type in ['rbo', 'mgs']:
        # Calculate gene score for each combination of predicted and ground truth programs
        for predicted_program in genes_topics_df.columns:
            for ground_truth_program in ground_truth_cell_type_gene_df.columns:
                query_list = genes_topics_df[predicted_program].sort_values(ascending=False).head(top_n_features).index.to_list()
                program_list = ground_truth_cell_type_gene_df[ground_truth_program].sort_values(ascending=False).head(top_n_features).index.to_list()
                score = get_ranking_score(query_list=query_list, program_list=program_list, rank_type=correlation_type)

                df_corr.at[predicted_program, ground_truth_program] = score

    elif correlation_type in ['pearson', 'spearman']:
        # Align indices and calculate Pearson correlation
        common_genes = ground_truth_cell_type_gene_df.index
        predicted_aligned = genes_topics_df.reindex(common_genes).fillna(0)
        
        for cell_type in ground_truth_cell_type_gene_df.columns:
            corr_cof = predicted_aligned.corrwith(ground_truth_cell_type_gene_df[cell_type], method=correlation_type).values
            df_corr.loc[:, cell_type] = corr_cof

    ## For Spots 
    elif correlation_type == "cosine":
        common_spots = ground_truth_cell_type_gene_df.index

        ## Note: fill spots that are missing by zeros , Note: Alternative to make the common spots just subset those genes. 
        predicted_aligned = genes_topics_df.reindex(common_spots).fillna(0) 

        for cell_type in ground_truth_cell_type_gene_df.columns:
            for topic in predicted_aligned.columns:

                column_value = cosine_similarity(ground_truth_cell_type_gene_df[cell_type].values.reshape(1, -1), predicted_aligned[topic].values.reshape(1, -1))
                df_corr.loc[topic, cell_type] = column_value[0,0]


    return df_corr

def annotate_with_benchmark(sample_results_dir, adata_spatial, exp_name, correlation_type = "pearson", top_n_features = 500):
    ### Used infered genes per topic and correlate it with main df
    sample_name = adata_spatial.uns["dataset_name"]

    results_dir_path = os.path.join(sample_results_dir, f'{exp_name}_{sample_name}')

    if(correlation_type == "cosine"):
        ## ground truth spots and spots_df
        spots_df = load_experiment_result(results_dir_path, sample_name=sample_name, exp_name = exp_name, mode ="spots")
        ground_truth = get_ground_truth(adata_spatial, mode="spots")
        df_corr = annotate_programs_by_ground_truth(spots_df,  ground_truth, correlation_type, top_n_features)
    else:
        ## ground truth genes and predicted genes
        gene_scores_df = load_experiment_result(results_dir_path, sample_name=sample_name, exp_name = exp_name, mode ="genescores")
        ground_truth = get_ground_truth(adata_spatial, mode="genes")
        df_corr = annotate_programs_by_ground_truth(gene_scores_df,  ground_truth, correlation_type, top_n_features)
    
    
    correlation_path = check_dir(os.path.join(results_dir_path,"correlations"))
    df_corr.to_csv(os.path.join(correlation_path, f"correlation_{correlation_type}_{exp_name}_{sample_name}.csv"))

    annotation_df = get_annotation_from_corr(df_corr)
    annotation_df.to_csv(os.path.join(results_dir_path, f"annotation_{exp_name}_{sample_name}.csv"))


def annot_corr_heatmap(
    results_dir, adata_spatial, experiments_list, correlation_type, is_triangle=False, is_show=True
):
    def plot_heatmap(ax, df, title, x_label, y_label, mask=None):
        """Helper function to plot a single heatmap."""
        sns.heatmap(
            df.astype(float),
            mask=mask,
            annot=False,
            cmap="Blues",
            ax=ax,
            square=True,
            cbar=False,
            vmin=min_corr,
            vmax=max_corr,
        )
        ax.set_title(title, fontsize=12)
        ax.set_xlabel(x_label, fontsize=10)
        ax.set_ylabel(y_label, fontsize=10)
        ax.set_xticks(np.arange(len(df.columns)) + 0.5)
        ax.set_xticklabels(df.columns, rotation=90, fontsize=8)
        ax.set_yticks(np.arange(len(df.index)) + 0.5)
        ax.set_yticklabels(df.index, rotation=0, fontsize=8)

    sample_name = adata_spatial.uns["dataset_name"]
    fig_width = 5 * (len(experiments_list) + 1)
    fig, axes = plt.subplots(
        nrows=1,
        ncols=len(experiments_list) + 1,
        figsize=(fig_width, 6),
        gridspec_kw={"wspace": 0.6},
    )

    plt_dir_corr = check_dir(os.path.join(results_dir, sample_name, "analysis", "correlations"))

    min_corr, max_corr = 0, 1  # Fixed legend limits
    ground_truth = get_ground_truth(adata_spatial, mode="spots")
    df_corr = annotate_programs_by_ground_truth(ground_truth, ground_truth, correlation_type)
    sorted_columns = sorted(df_corr.columns)
    df_corr = df_corr[sorted_columns].reindex(sorted_columns)
    df_corr.to_csv(os.path.join(plt_dir_corr, f"{sample_name}_{correlation_type}.csv"))

    # Mask for upper triangle if required
    mask = np.triu(np.ones_like(df_corr, dtype=bool), k=1) if is_triangle else None

    # Plot the ground truth heatmap
    plot_heatmap(
        axes[0],
        df_corr,
        title="Ground Truth",
        x_label="Cell Type",
        y_label="Cell Type",
        mask=mask,
    )

    # Process and plot each experiment's heatmap
    for i, exp_name in enumerate(experiments_list, start=1):
        annotate_path = os.path.join(results_dir, sample_name, f"{exp_name}_{sample_name}")
        corr_file = os.path.join(
            annotate_path,
            "correlations",
            f"correlation_{correlation_type}_{exp_name}_{sample_name}.csv",
        )
        df_corr = pd.read_csv(corr_file, index_col=0)
        df_annotation = get_annotation_from_corr(df_corr).sort_values(by="celltype")
        df_reordered = df_corr[df_annotation["celltype"].values]
        df_reordered = df_reordered.reindex(df_annotation["program"].values).T

        mask = np.triu(np.ones_like(df_reordered, dtype=bool), k=1) if is_triangle else None
        plot_heatmap(
            axes[i],
            df_reordered,
            title=exp_name,
            x_label="Programs (Predicted)",
            y_label="Cell Type (Ground Truth)",
            mask=mask,
        )

    # Add a single shared colorbar
    cbar = fig.colorbar(
        plt.cm.ScalarMappable(cmap="Blues", norm=plt.Normalize(vmin=min_corr, vmax=max_corr)),
        ax=axes,
        orientation="vertical",
        fraction=0.02,
        pad=0.02,
    )
    cbar.set_label("Correlation", fontsize=10)
    cbar.ax.tick_params(labelsize=8)

    # Adjust layout and save the figure
    plt.tight_layout(rect=[0, 0, 0.9, 0.95])
    plt.suptitle(f"Cosine Similarity (Correlations) {sample_name}", fontsize=16, y=1.02)
    plt.savefig(
        os.path.join(plt_dir_corr, f"{sample_name}_annotation_correlation_heatmap_{correlation_type}.pdf"),
        dpi=300,
    )
    if is_show:
        plt.show()
    plt.close()

def benchmark_corr_silverstandard(results_dir_path, silver_standard_df,
                                  rf_usages, correlation_type = "cosine",
                                    title = "Correlation Heatmap", annotation_df=None, is_show=False):
    
    plt_dir_path = check_dir(os.path.join(results_dir_path, "SS_Annotation"))
    common_index = silver_standard_df.index.intersection(rf_usages.index)
    # Reindex both DataFrames to keep only common indices
    silver_standard_df = silver_standard_df.loc[common_index]
    rf_usages = rf_usages.loc[common_index]


    df_corr = pd.DataFrame(index=silver_standard_df.columns, columns=rf_usages.columns)
    for cell_type in rf_usages.columns:
        for topic in silver_standard_df.columns:
            if correlation_type == "cosine":
                column_value = cosine_similarity(
                    rf_usages[cell_type].values.reshape(1, -1),
                    silver_standard_df[topic].values.reshape(1, -1)
                )[0, 0]
            elif correlation_type == "pearson":
                column_value, _ = pearsonr(rf_usages[cell_type], silver_standard_df[topic])
            elif correlation_type == "spearman":
                column_value, _ = spearmanr(rf_usages[cell_type], silver_standard_df[topic])
            else:
                raise ValueError("Unsupported correlation type: choose from 'cosine', 'pearson', 'spearman'")
            
            df_corr.loc[topic, cell_type] = column_value

    df_annotate = get_annotation_from_corr(df_corr)
    sorted_programs = df_annotate.sort_values(by="celltype")["program"].values
    sorted_columns = sorted(df_corr.columns)
    df_corr = df_corr[sorted_columns].reindex(sorted_programs)

    df_corr.to_csv(os.path.join(plt_dir_path, f"ss_correlation_{title}_{correlation_type}.csv"))
    df_annotate.to_csv(os.path.join(plt_dir_path, f"ss_annotation_{title}_{correlation_type}.csv"))

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(15, 15))

    # Plot heatmap
    sns.heatmap(
        df_corr.astype(float),
        annot=False,
        cmap="Blues",
        ax=ax,
        square=True,
        cbar=False,
        vmin=0,
        vmax=1,
    )
    x_label = "Predicted Programs"
    y_label = "Silver Standard Annotation"
    ax.set_title(title, fontsize=12)
    ax.set_xlabel(x_label, fontsize=10)
    ax.set_ylabel(y_label, fontsize=10)

    # Ticks
    ax.set_xticks(np.arange(len(df_corr.columns)) + 0.5)
    ax.set_xticklabels(df_corr.columns, rotation=90, fontsize=8)
    ax.set_yticks(np.arange(len(df_corr.index)) + 0.5)
    ax.set_yticklabels(df_corr.index, rotation=0, fontsize=8)
    # Show plot
    plt.savefig(
        os.path.join(plt_dir_path, f"ss_correlation_heatmap_{title}_{correlation_type}.pdf"),
        dpi=300,
    )
    if is_show:
        plt.show()
    plt.close()


    unique_cats = annotation_df["annotation"].unique()
    palette = sns.color_palette("Set2", n_colors=len(unique_cats))
    category_colors = dict(zip(unique_cats, palette))
    row_colors = annotation_df["annotation"].map(category_colors)

    # Create Clustermap
    g = sns.clustermap(
        df_corr.astype(float),
        row_colors=row_colors,
        cmap="Blues",
        figsize=(15, 15),
        vmin=0,
        vmax=1,
        col_cluster=False,
        cbar_pos=(0.90, 0.80, 0.03, 0.18),  # (x, y, width, height)
        # if you wanted to disable the colorbar:
        # cbar_pos=None,
        # or to tweak it:
        # cbar_kws={"shrink": 0.5, "label": "correlation"}
    )

    # 4. Create small legend handles
    handles = [
        Patch(facecolor=category_colors[cat], edgecolor="none", label=cat)
        for cat in unique_cats
    ]

    # 5. Add a compact legend inside the plotting area
    g.ax_heatmap.legend(
        handles=handles,
        title="annotation",
        loc="upper right",
        bbox_to_anchor=(1.02, 1.15),   # y = 1.15 puts it 15% above the top edge
        borderaxespad=0.0,
        frameon=False,
        fontsize=8,
        title_fontsize=9,
        handlelength=1,
        handleheight=1
    )

    # Apply reordered labels
    g.ax_heatmap.set_xticks(np.arange(len(df_corr.columns)) + 0.5)
    g.ax_heatmap.set_yticks(np.arange(len(df_corr.index)) + 0.5)
    g.ax_heatmap.set_yticklabels(df_corr.index[g.dendrogram_row.reordered_ind], rotation=0, fontsize=8)

    # Save and show clustermap
    plt.savefig(os.path.join(plt_dir_path, f"ss_clustermap_{title}_{correlation_type}.pdf"), dpi=300)
    if is_show:
        plt.show()
    plt.close()


def benchmark_corr_silverstandard2(results_dir_path, counts_df, rf_scores, cell_indices_dict, title = "Correlation Heatmap", correlation_type="pearson", is_show=False):
    # Get unique cell types and topics
    plt_dir_path = check_dir(os.path.join(results_dir_path, "SS_Annotation2"))

    cell_types = list(cell_indices_dict.keys())
    topics = rf_scores.columns
    
    # Initialize an empty DataFrame for storing correlations
    df_corr = pd.DataFrame(index=topics, columns=cell_types, dtype=float)

    # Loop through each topic and cell type
    for cell_type in tqdm(cell_types, desc="Cell types"):
        cell_indices_all = cell_indices_dict[cell_type]
        plasma_counts_clean = counts_df.loc[cell_indices_all]
        for topic in topics:
            aligned_rf_clean = rf_scores[[topic]]
    
            # Align indices
            common_cells = aligned_rf_clean.index.intersection(plasma_counts_clean.columns)
    
            plasma_counts_aligned = plasma_counts_clean[common_cells]
            aligned_rf_aligned = aligned_rf_clean.loc[common_cells, topic]
    
            # Compute correlation
            correlations = plasma_counts_aligned.corrwith(aligned_rf_aligned, axis=1)
            df_corr.loc[topic, cell_type] = correlations.mean()

    # Create the DataFrame and pivot it for visualization

    df_annotate = get_annotation_from_corr(df_corr)
    sorted_programs = df_annotate.sort_values(by="celltype")["program"].values
    sorted_columns = sorted(df_corr.columns)
    df_corr = df_corr[sorted_columns].reindex(sorted_programs)
    df_annotate.to_csv(os.path.join(plt_dir_path, f"ss_annotation_{title}_{correlation_type}.csv"))
    df_corr.to_csv(os.path.join(plt_dir_path, f"ss_correlation_{title}_{correlation_type}.csv"))

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(15, 15))

    # Plot heatmap
    sns.heatmap(
        df_corr.astype(float),
        annot=False,
        cmap="Blues",
        ax=ax,
        square=True,
        cbar=False,
        vmin=0,
        vmax=1,
    )
    x_label = "Predicted Programs"
    y_label = "Silver Standard Annotation"
    ax.set_title(title, fontsize=12)
    ax.set_xlabel(x_label, fontsize=10)
    ax.set_ylabel(y_label, fontsize=10)

    # Ticks
    ax.set_xticks(np.arange(len(df_corr.columns)) + 0.5)
    ax.set_xticklabels(df_corr.columns, rotation=90, fontsize=8)
    ax.set_yticks(np.arange(len(df_corr.index)) + 0.5)
    ax.set_yticklabels(df_corr.index, rotation=0, fontsize=8)
    # Show plot
    plt.savefig(
        os.path.join(plt_dir_path, f"ss_correlation_heatmap_{title}_{correlation_type}.pdf"),
        dpi=300,
    )
    if is_show:
        plt.show()
    plt.close()

    # Create Clustermap
    g = sns.clustermap(
        df_corr.astype(float),
        cmap="Blues",
        figsize=(15, 15),
        cbar=True,
        vmin=0,
        vmax=1
    )

    # Get reordered indices
    row_order = g.dendrogram_row.reordered_ind
    col_order = g.dendrogram_col.reordered_ind

    # Apply reordered labels
    g.ax_heatmap.set_xticks(np.arange(len(df_corr.columns)) + 0.5)
    g.ax_heatmap.set_xticklabels(df_corr.columns[col_order], rotation=90, fontsize=8)
    g.ax_heatmap.set_yticks(np.arange(len(df_corr.index)) + 0.5)
    g.ax_heatmap.set_yticklabels(df_corr.index[row_order], rotation=0, fontsize=8)

    # Save and show clustermap
    plt.savefig(os.path.join(plt_dir_path, f"ss_clustermap_{title}_{correlation_type}.pdf"), dpi=300)
    if is_show:
        plt.show()
    plt.close()


