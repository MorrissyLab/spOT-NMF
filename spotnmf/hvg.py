import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import scipy.sparse as sp
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from pygam import s, LinearGAM 
from scipy.stats import chi2, f
from statsmodels.formula.api import ols
from statsmodels.sandbox.stats.multicomp import multipletests

# Patch for pygam compatibility with scipy sparse matrices
def to_array(self):
    return self.toarray()
sp.spmatrix.A = property(to_array)


def compute_overdispersed_genes_batches(adata, batch_keys=["sample_id"], union_agg=False, **kwargs):
    """
    Computes a union set of overdispersed genes across all samples in adata, 
    using a specified batch_key for grouping.

    Parameters:
    - adata: AnnData object
    - batch_keys: list[str], column(s) in adata.obs to group samples by
    - union_agg: bool, whether to also compute on aggregated data
    - **kwargs: additional keyword arguments passed to compute_overdispersed_genes

    Returns:
    - A list of unique overdispersed genes across samples
    """
    odg_list = []
    
    if union_agg:
        # Global computation
        print(f"Calculating ODG of aggregated samples")
        _, overdispersed_genes = compute_overdispersed_genes(adata.copy(), **kwargs)
        odg_list.extend(overdispersed_genes)

    # Per-sample computation
    for batch_key in batch_keys:
        for sample in adata.obs[batch_key].unique():
            print(f"Calculating ODG of {sample}")
            _, overdispersed_genes = compute_overdispersed_genes(
                adata[adata.obs[batch_key] == sample].copy(), **kwargs
            )
            odg_list.extend(overdispersed_genes)

    return list(set(odg_list))
    
def compute_overdispersed_genes(
    adata_spatial,
    save_dir = None,
    gam_k=5,
    alpha=0.05,
    max_adjusted_variance=1e3,
    min_adjusted_variance=1e-3,
    removeAbove=1.0,
    removeBelow=0.05,
    min_counts=1,
    min_genes=1,
    n_top_genes=None,
    use_unadjusted_pvals=False,
    is_spatial=False,
    is_show=False,
    verbose=True,
):
    """
    Filter data and Identifies overdispersed genes from spatial transcriptomic data.

    Parameters:
    -----------
    adata_spatial : AnnData
        Annotated data matrix containing spatial data.
    gam_k : int, optional (default=5)
        Smoothing parameter for generalized additive model (GAM).
    alpha : float, optional (default=0.05)
        Significance level for identifying overdispersed genes.
    max_adjusted_variance : float, optional (default=1e3)
        Maximum allowable variance for adjusted values.
    min_adjusted_variance : float, optional (default=1e-3)
        Minimum allowable variance for adjusted values.
    removeAbove : float, optional (default=1.0)
        Threshold above which genes are excluded based on max cell count.
    removeBelow : float, optional (default=0.05)
        Threshold below which genes are excluded based on min cell count.
    min_counts : int, optional (default=1)
        Minimum counts threshold for filtering spots.
    min_genes : int, optional (default=1)
        Minimum number of genes for filtering spots.
    n_top_genes : int or None, optional (default=None)
        Maximum number of top overdispersed genes to return.
    use_unadjusted_pvals : bool, optional (default=False)
        If True, use unadjusted p-values instead of adjusted p-values.
    is_show : bool, optional (default=True)
        If True, plot mean-variance relationships and other statistics.
    verbose : bool, optional (default=True)
        If True, display progress and debugging information.

    Returns:
    --------
    adata_spatial : AnnData
        The updated AnnData object with overdispersion information in `var` attribute.
    odg : list of str
        List of names of overdispersed genes.
    """
    print(f'Selecting Genes with alpha {alpha} and use_unadjusted_pvals {use_unadjusted_pvals} is_spatial {is_spatial}')
    # Initial data summary

    print(f"Initial data: {adata_spatial.X.shape[0]} spots, {adata_spatial.X.shape[1]} genes")

    # Filter spots based on minimum counts and genes per spot
    sc.pp.filter_cells(adata_spatial, min_counts=min_counts, inplace=True)
    sc.pp.filter_cells(adata_spatial, min_genes=min_genes, inplace=True)
    print(f"Spots after filtering: {adata_spatial.X.shape[0]}")

    # Filter genes based on occurrence thresholds across spots
    sc.pp.filter_genes(adata_spatial, min_cells=int(removeBelow * len(adata_spatial)), inplace=True)
    sc.pp.filter_genes(adata_spatial, max_cells=int(removeAbove * len(adata_spatial)), inplace=True)
    print(f"Genes after filtering: {adata_spatial.X.shape[1]}")

    # Prepare data matrix and gene labels
    if(is_spatial):
        import SpatialDE as sd
        df, df_test_detailed = sd.test(adata_spatial)
        df.set_index("gene", inplace=True)
        df.index.names = ['genes']
        df.sort_values("padj", inplace=True, ascending=True)
        df['is_odg'] = df['pval'].lt(alpha) if use_unadjusted_pvals else df['padj'].lt(alpha)

        # Extract overdispersed genes
        odg = list(df[df['is_odg'] == True].index)
        if verbose:
            print(f"Identified {len(odg)} overdispersed genes.")

    else:
        mat = adata_spatial.X
        filtered_genes = adata_spatial.var.index.values
        n_spots = len(adata_spatial)
        n_genes = len(filtered_genes)

        # Compute log mean and variance for each gene
        dfm = np.log(np.mean(mat, axis=0))
        dfv = np.log(np.var(mat, axis=0))
        df = pd.DataFrame({'m': dfm,'v': dfv, 'res': -np.inf}, index=filtered_genes)
        vi = np.where(np.isfinite(dfv))[0]  # indices of valid variance values
        
        # Adjust GAM smoothing parameter if insufficient data points
        if len(vi) < gam_k * 1.5:
            gam_k = 1

        # Fit model to variance
        if gam_k < 2:
            if verbose:
                print("Using linear model due to insufficient data for GAM ...")
            model = ols('v ~ m', data=df.iloc[vi]).fit()
            df.loc[df.index[vi], 'res'] = model.resid.values
        else:
            if verbose:
                print(f"Using GAM model with k={gam_k} ...")
            model = LinearGAM(s(0, n_splines=gam_k)).fit(df.iloc[vi]['m'], df.iloc[vi]['v'])
            df.loc[df.index[vi], 'res'] = model.deviance_residuals(df.iloc[vi]['m'].values, df.iloc[vi]['v'].values)

        # Calculate p-values and adjust for multiple testing
        df['p_value'] = df.apply(lambda row: 1 - f.cdf(np.exp(row['res']), n_spots, n_spots), axis=1) # Fix variable spots
        df['p_value_adj'] = multipletests(df['p_value'].values, method='fdr_bh')[1]
        df['chi_qv'] = chi2.ppf(1 - df['p_value'].values, n_genes - 1) / n_genes

        # Scale gene expression variances
        df['gene_scale_factor'] = np.sqrt(np.maximum(min_adjusted_variance, np.minimum(max_adjusted_variance, df['chi_qv'])) / np.exp(df['v']))
        df['gene_scale_factor'] = np.where(np.isfinite(df['gene_scale_factor']), df['gene_scale_factor'], 0)

        # Identify overdispersed genes based on significance level
        df['is_odg'] = df['p_value'].lt(alpha) if use_unadjusted_pvals else df['p_value_adj'].lt(alpha)


        # Extract overdispersed genes
        odg = list(df[df['is_odg'] == True].index)
        if verbose:
            print(f"Identified {len(odg)} overdispersed genes.")

        # Plot distributions of genes per spot and spots per gene
        if(save_dir):
            os.makedirs(save_dir, exist_ok=True)

            fig, axs = plt.subplots(1, 2, figsize=(10, 5))
            axs[0].hist(np.log10(mat.sum(axis=0) + 1), bins=20)
            axs[0].set_title("Genes Per Spot Distribution")
            
            axs[1].hist(np.log10(mat.sum(axis=1) + 1), bins=20)
            axs[1].set_title("Spots Per Gene Distribution")

            plt.savefig(os.path.join(save_dir, f"odg_genes_spots_distribution_plot.pdf"), dpi=300)
            if is_show:
                plt.show()
            plt.close()

            # Optional plotting of relationships
            fig, axes = plt.subplots(1, 2, figsize=(12, 6))
            # Mean vs variance plot
            sns.kdeplot(x=df['m'], y=df['v'], cmap="Blues", ax=axes[0], bw_adjust=0.5, fill=True)
            axes[0].scatter(x=df['m'], y=df['v'], color="black", s=0.5)
            axes[0].set(xlabel='Mean', ylabel='Variance', title='Mean vs Variance (Genes)')
            
            # Linear regression model plot
            grid = np.linspace(df['m'].iloc[vi].min(), df['m'].iloc[vi].max(), 1000)
            model = LinearRegression().fit(df[['m']].values, df['v'].values)
            axes[0].plot(grid, model.predict(grid.reshape(-1, 1)), color="blue")
            
            if len(odg) > 0:
                axes[0].scatter(df[df['is_odg'] == True]['m'], df[df['is_odg'] == True]['v'], color="red", s=5, marker='.')

            # Mean vs chi_qv plot
            sns.kdeplot(x=df['m'].iloc[vi], y=df['chi_qv'].iloc[vi], cmap="Blues", ax=axes[1], bw_adjust=0.5, fill=True)
            axes[1].scatter(x=df['m'].iloc[vi], y=df['chi_qv'].iloc[vi], color="black", s=0.5)
            axes[1].set(xlabel='Mean', ylabel='Chi_qv', title='Mean vs Chi_qv')
            axes[1].axhline(y=1, color='gray', linestyle='--')
            axes[1].scatter(df[df['is_odg'] == True]['m'], df[df['is_odg'] == True]['chi_qv'], color="red", s=5, marker='.')
            
            plt.savefig(os.path.join(save_dir, f"odg_mean_variance_plot.pdf"), dpi=300)
            if is_show:
                plt.show()
            plt.close('all')


        
        # Update AnnData object with overdispersion results
        adata_spatial.norm_mat = mat * df['gene_scale_factor'].T.values


    adata_spatial.var = adata_spatial.var.join(df)
    adata_spatial.var["is_odg"] = adata_spatial.var["is_odg"].fillna(False)

    # Select Top genes
    if n_top_genes:
        odg = odg[:n_top_genes]

    return adata_spatial, odg

def save_hvg_list(l, file):
    """
    Save a list of highly variable genes (HVGs) to a CSV file.

    Parameters
    ----------
    l : list
        List of gene names or identifiers to be saved.
    file : str
        Path to the CSV file where the list will be saved.

    Returns
    -------
    None
    """
    pd.DataFrame(l).to_csv(file, index=False, header=False)


def load_hvg_list(file_path):
    """
    Load a list of highly variable genes (HVGs) from a CSV file.

    Parameters
    ----------
    file_path : str
        Path to the CSV file containing the list of HVGs.

    Returns
    -------
    list
        List of gene names or identifiers.
    """
    return pd.read_csv(file_path, header=None)[0].tolist()