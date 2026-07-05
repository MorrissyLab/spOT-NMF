import numpy as np
from scipy.stats import zscore
import statsmodels.api as sm
import pandas as pd
pd.options.display.float_format = '{:f}'.format


def compute_glm_coefficients(Z, U, family=sm.families.Gaussian()):
    """Compute regression coefficients for each gene using a Generalized Linear Model (GLM).

    A separate GLM is fit for each gene independently.

    Args:
        Z (numpy.ndarray): Z-scored gene expression matrix of shape
            (num_samples, num_genes).
        U (numpy.ndarray): Un-normalized consensus usage matrix of shape
            (num_samples, num_gep).
        family: statsmodels GLM family. Defaults to Gaussian (equivalent to OLS).

    Returns:
        numpy.ndarray: GLM coefficients for each gene, of shape
        (num_genes, num_gep + 1).
    """

    B = []  # Initialize list to store coefficients for each gene
    # Fit a GLM for each gene independently
    for j in range(Z.shape[1]):
        # Extract the target variable for the j-th gene (column vector)
        target = Z[:, j]
        # Define and fit the GLM model
        model = sm.GLM(target, U, family=family)
        result = model.fit()

        # Append the estimated coefficients to B
        B.append(result.params)

    # Convert B to a numpy array for consistency, with shape (num_genes, num_gep + 1)
    B = np.array(B)
    return B


def calculate_marker_genes_topics_df(adata_spatial, rf_usages, model_type="ols"):
    """Calculate regression coefficients of marker genes for each topic.

    Performs Ordinary Least Squares (OLS) regression or a Generalized Linear Model
    to compute regression coefficients of marker genes for each topic, based on the
    spatial data and reference usages.

    Args:
        adata_spatial (anndata.AnnData): An AnnData object containing spatial
            transcriptomic data with the gene expression matrix (``X``) and spot
            metadata (``obs``).
        rf_usages (pandas.DataFrame): Reference usages for each topic, indexed by
            spot identifiers, aligning with ``adata_spatial.obs.index``.
        model_type (str): Type of regression model, either ``"ols"`` or ``"glm"``.
            Defaults to ``"ols"``.

    Returns:
        pandas.DataFrame: Regression coefficients for each gene-topic pair. Rows
        represent genes and columns represent topics.
    """
    
    # Normalize gene expression matrix (z-score normalization: (T - mean) / std) for each gene
    adata_spatial = adata_spatial[rf_usages.index, :]
    Z = zscore(adata_spatial.X)
    U = rf_usages.reindex(adata_spatial.obs.index).values

    if(model_type == "ols"):
        # Calculate OLS regression coefficients: B = (U^T * U)^-1 * U^T * Z
        B = (np.linalg.pinv(U) @ Z).T
    elif(model_type == "glm"):
        B = compute_glm_coefficients(Z, U, family=sm.families.Gaussian())
        print('h')
    
    # Create DataFrame of regression coefficients, with genes as rows and topics as columns
    usage_coef_df = pd.DataFrame(B, index=adata_spatial.var.index, columns=rf_usages.columns)
    
    return usage_coef_df
