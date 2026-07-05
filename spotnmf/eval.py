import pandas as pd
import rbo
from sklearn.metrics.pairwise import cosine_similarity

def get_ranking_score(query_list, program_list, rank_type='rbo'):
    """Compute a ranking-based similarity score between two ranked gene lists.

    Args:
        query_list (list): Reference (ground truth) ranked list of gene names.
        program_list (list): Predicted ranked list of gene names.
        rank_type (str): Type of ranking score to compute. Defaults to
            ``'rbo'``. Options:

            - ``'rbo'``: Rank-Biased Overlap (standard).
            - ``'rboext'``: Extended (extrapolated) Rank-Biased Overlap.
            - ``'mgs'``: Mean Gene Score (sum of inverse ranks of genes in
              ``program_list`` that also appear in ``query_list``).

    Returns:
        float: The computed ranking similarity score.

    Raises:
        ValueError: If ``rank_type`` is not one of the supported options.
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
    """Greedily assign each program to its best-matching cell type.

    Flattens the correlation/similarity matrix, sorts all pairs in descending
    order, and assigns programs to cell types one-to-one, ensuring each program
    and each cell type is used at most once.

    Args:
        df_corr (pandas.DataFrame): Correlation or similarity matrix indexed by
            program with cell types as columns.

    Returns:
        pandas.DataFrame: Assignments with columns ``['program', 'celltype']``.
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
    """Score predicted programs against reference cell types.

    Computes a pairwise similarity/correlation matrix between predicted gene
    programs and reference (ground truth) cell types using the chosen metric.
    Ranking-based metrics operate on the top-N feature lists, ``'pearson'`` and
    ``'spearman'`` correlate full aligned gene vectors, and ``'cosine'``
    computes cosine similarity between aligned vectors.

    Args:
        genes_topics_df (pandas.DataFrame): Gene (or spot) scores per predicted
            program, with features as the index and programs as columns.
        ground_truth_cell_type_gene_df (pandas.DataFrame): Reference scores per
            cell type, with features as the index and cell types as columns.
        correlation_type (str): Similarity/correlation metric. Defaults to
            ``'pearson'``. One of ``'pearson'``, ``'spearman'``, ``'rbo'``,
            ``'rboext'``, ``'mgs'`` or ``'cosine'``.
        top_n_features (int): Number of top features used for ranking-based
            metrics (``'rbo'``, ``'rboext'``, ``'mgs'``). Defaults to 500.

    Returns:
        pandas.DataFrame: Pairwise scores indexed by program with cell types as
        columns.

    Raises:
        ValueError: If ``correlation_type`` is not supported.
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
