import pandas as pd
import rbo
from sklearn.metrics.pairwise import cosine_similarity

def get_ranking_score(query_list, program_list, rank_type='rbo'):
    """
    Compute a ranking-based similarity score between two ranked gene lists.

    Parameters
    ----------
    query_list : list of str
        The reference (ground truth) ranked list of genes.
    program_list : list of str
        The predicted ranked list of genes.
    rank_type : str, default='rbo'
        Type of ranking score to compute. Options:
        - 'rbo'    : Rank Biased Overlap (standard)
        - 'rboext' : Extended RBO
        - 'mgs'    : Mean Gene Score (inverse-rank sum of matches)

    Returns
    -------
    score : float
        The computed ranking similarity score.
    """
    # Find intersection and ranks for MGS
    intersected_genes = [(x, i + 1) for i, x in enumerate(program_list) if x in query_list]

    if rank_type == 'rbo':
        score = rbo.RankingSimilarity(query_list, program_list).rbo()
    elif rank_type == 'rboext':
        score = rbo.RankingSimilarity(query_list, program_list).rbo_ext()
    elif rank_type == 'mgs':
        score = sum(1 / rank for _, rank in intersected_genes if rank != 0)
    else:
        raise ValueError(f"Unsupported rank_type: {rank_type}")
    return score

def get_annotation_from_corr(df_corr):
    """
    Given a correlation/similarity matrix, assign each predicted program to the best-matching cell type
    (one-to-one assignment, greedy).

    Parameters
    ----------
    df_corr : pd.DataFrame
        Correlation or similarity matrix (programs x cell types).

    Returns
    -------
    annotation_df : pd.DataFrame
        DataFrame with columns ['program', 'celltype'] representing the assignments.
    """
    annotation = []
    # Flatten and sort all correlations/similarities descendingly
    sorted_correlations = df_corr.stack().sort_values(ascending=False)
    used_celltypes = set()
    used_programs = set()
    for (program, celltype), _ in sorted_correlations.items():
        if celltype not in used_celltypes and program not in used_programs:
            annotation.append([program, celltype])
            used_celltypes.add(celltype)
            used_programs.add(program)
    return pd.DataFrame(annotation, columns=["program", "celltype"])

def annotate_programs_by_ground_truth(
    genes_topics_df,
    ground_truth_cell_type_gene_df,
    correlation_type='pearson',
    top_n_features=500
):
    """
    Annotate predicted gene programs with reference cell types using different similarity/correlation metrics.

    Parameters
    ----------
    genes_topics_df : pd.DataFrame
        DataFrame of gene scores per program (genes x programs).
    ground_truth_cell_type_gene_df : pd.DataFrame
        DataFrame of gene scores per ground truth cell type (genes x cell types).
    correlation_type : str, default='pearson'
        Metric for similarity/correlation: 'pearson', 'spearman', 'rbo', 'mgs', 'cosine'.
    top_n_features : int, default=500
        Number of top genes/features to use for ranking-based metrics.

    Returns
    -------
    df_corr : pd.DataFrame
        Matrix of pairwise scores between predicted programs and ground truth cell types.
    """
    df_corr = pd.DataFrame(index=genes_topics_df.columns, columns=ground_truth_cell_type_gene_df.columns)

    if correlation_type in ['rbo', 'mgs', 'rboext']:
        # Ranking-based metrics: operate on top-N gene lists
        for pred_prog in genes_topics_df.columns:
            for gt_prog in ground_truth_cell_type_gene_df.columns:
                query_list = genes_topics_df[pred_prog].sort_values(ascending=False).head(top_n_features).index.to_list()
                program_list = ground_truth_cell_type_gene_df[gt_prog].sort_values(ascending=False).head(top_n_features).index.to_list()
                score = get_ranking_score(query_list, program_list, rank_type=correlation_type)
                df_corr.at[pred_prog, gt_prog] = score

    elif correlation_type in ['pearson', 'spearman']:
        # Vector correlations (across all genes)
        common_genes = ground_truth_cell_type_gene_df.index
        predicted_aligned = genes_topics_df.reindex(common_genes).fillna(0)
        for gt_celltype in ground_truth_cell_type_gene_df.columns:
            corr_values = predicted_aligned.corrwith(
                ground_truth_cell_type_gene_df[gt_celltype], method=correlation_type
            )
            df_corr[gt_celltype] = corr_values.values

    elif correlation_type == "cosine":
        # Cosine similarity across spot/feature vectors
        common_spots = ground_truth_cell_type_gene_df.index
        predicted_aligned = genes_topics_df.reindex(common_spots).fillna(0)
        for gt_celltype in ground_truth_cell_type_gene_df.columns:
            for topic in predicted_aligned.columns:
                sim = cosine_similarity(
                    ground_truth_cell_type_gene_df[gt_celltype].values.reshape(1, -1),
                    predicted_aligned[topic].values.reshape(1, -1)
                )
                df_corr.at[topic, gt_celltype] = sim[0, 0]
    else:
        raise ValueError(f"Unsupported correlation_type: {correlation_type}")

    return df_corr
