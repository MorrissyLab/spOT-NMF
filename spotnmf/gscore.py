import numpy as np
from scipy.stats import zscore
import statsmodels.api as sm
import pandas as pd
pd.options.display.float_format = '{:f}'.format


def compute_glm_coefficients(Z, U, family=sm.families.Gaussian()):
    """
    Computes the regression coefficients for each gene using a Generalized Linear Model (GLM).
    
    Parameters:
    - Z: numpy array of shape (num_samples, num_genes), the z-scored gene expression matrix.
    - U: numpy array of shape (num_samples, num_gep), the un-normalized consensus usage matrix.
    - family: statsmodels GLM family, default is Gaussian (OLS will be similliar).
    
    Returns:
    - B: numpy array of shape (num_genes, num_gep + 1), GLM coefficients for each gene.
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
    """
    Performs Ordinary Least Squares (OLS) regression or General Linar Model to calculate regression coefficients 
    of marker genes for each topic, based on the spatial data and reference usages.
    
    Parameters:
    - adata_spatial: AnnData
        An AnnData object containing spatial transcriptomic data with gene expression matrix (X) 
        and spot metadata (obs).
    - rf_usages: DataFrame
        A DataFrame with reference usages for each topic, indexed by spot identifiers, 
        aligning with `adata_spatial.obs.index`.
    - model_type: string
        type of regression model
    
    Returns:
    - usage_coef_df: DataFrame
        A DataFrame containing the regression coefficients for each gene-topic pair. 
        Rows represent genes, and columns represent topics, with coefficients computed 
        using OLS regression.
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
