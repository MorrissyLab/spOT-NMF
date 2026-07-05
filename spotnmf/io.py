import os
import gc
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.image import imread
import anndata as ad

def check_dir(dir_path):
    """Create the directory if it doesn't exist, and return its path.

    Args:
        dir_path (str): Path of the directory to create.

    Returns:
        str: The same directory path.
    """
    os.makedirs(dir_path, exist_ok=True)
    return dir_path

def read_adata(
    data_path,
    data_mode,
    genome = "mm10",
    is_aggr=False,
    is_xenograft=False,
    select_sample=None,
    bin_size=None
):
    """Read an AnnData object based on the data mode and parameters.

    Dispatches to the Visium HD or Visium spatial readers based on
    ``data_mode``, falling back to reading a plain ``.h5ad`` file.

    Args:
        data_path (str): Path to the AnnData or spatial data.
        data_mode (str): Mode of the data, e.g. 'visium_hd', 'visium', or 'h5ad'.
        genome (str): Genome reference name. Defaults to "mm10".
        is_aggr (bool): Whether the data is aggregated. Defaults to False.
        is_xenograft (bool): Whether the data is from a xenograft model. Defaults to False.
        select_sample (str): Specific sample name to select. Defaults to None.
        bin_size (int): Bin size for spatial data. Defaults to None.

    Returns:
        anndata.AnnData: Loaded AnnData object.
    """
    if data_mode == "visium_hd":
        adata_spatial = read_visium_hd(
            adata_path=data_path,
            bin_size=bin_size,
            genome=genome,
            is_aggr=is_aggr,
            is_xenograft=is_xenograft
        )
        return adata_spatial

    elif data_mode == "visium":
        adata_spatial = read_spatial_data(
            adata_path=data_path,
            genome=genome,
            is_xenograft=is_xenograft,
            is_aggr=is_aggr,
            select_sample=select_sample
        )
        return adata_spatial 
       
    else:
        adata = sc.read_h5ad(data_path)
        return adata
    

def load_spatial_image(image_path, mode='original'):
    """Load a spatial image and return it in the desired color mode.

    Args:
        image_path (str): Path to the image file.
        mode (str): Mode for processing the image, one of 'grayscale',
            'original', or 'interpolated'. Defaults to 'original'.

    Returns:
        numpy.ndarray: Image array of shape (H, W, 3) in uint8 format.

    Raises:
        ValueError: If ``mode`` is not 'grayscale', 'original', or 'interpolated'.
    """
    img = imread(image_path)
    if img.dtype == np.uint8:
        img = img / 255.0
    if img.shape[-1] == 4:
        img = img[:, :, :3]

    if mode == 'original':
        return (img * 255).astype(np.uint8)

    # Convert to grayscale (luminosity method)
    gray = np.dot(img[..., :3], [0.2989, 0.5870, 0.1140])
    if mode == 'grayscale':
        gray_rgb = np.stack([gray] * 3, axis=-1)
        return (gray_rgb * 255).astype(np.uint8)
    elif mode == 'interpolated':
        cornsilk_rgb = np.array([255, 248, 220]) / 255.0
        white_rgb = np.array([1.0, 1.0, 1.0])
        tinted = (1 - gray[..., None]) * cornsilk_rgb + gray[..., None] * white_rgb
        return (tinted.clip(0, 1) * 255).astype(np.uint8)
    else:
        raise ValueError("mode must be 'grayscale', 'original', or 'interpolated'")

def read_visium_hd_sample(
    sample_adata_path, bin_size=16, genome=None,
    library_id="visium_hd_sample", is_xenograft=False, image_color="grayscale"
):
    """Read a single Visium HD sample and return an AnnData object.

    Args:
        sample_adata_path (str): Path to the sample directory.
        bin_size (int): Spatial bin size in microns. Defaults to 16.
        genome (str): Reference genome, used only for xenograft support. Defaults to None.
        library_id (str): Library/sample identifier. Defaults to "visium_hd_sample".
        is_xenograft (bool): Whether the data is from a xenograft model. Defaults to False.
        image_color (str): How to load tissue images, one of 'grayscale',
            'original', or 'interpolated'. Defaults to "grayscale".

    Returns:
        anndata.AnnData: Annotated data object for the sample.
    """
    bin_sampled_path = os.path.join(
        sample_adata_path, f"binned_outputs/square_{int(bin_size):03d}um"
    )

    adata = sc.read_10x_h5(
        os.path.join(bin_sampled_path, "filtered_feature_bc_matrix.h5"),
        genome=genome if is_xenograft else None
    )

    # Join spatial coordinates
    df_spatial = pd.read_parquet(os.path.join(bin_sampled_path, "spatial/tissue_positions.parquet"))
    df_spatial.set_index("barcode", inplace=True)
    adata.obs = adata.obs.join(df_spatial, how="left")
    adata.obsm["spatial"] = adata.obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy().astype(float)
    adata.obs.drop(columns=["pxl_row_in_fullres", "pxl_col_in_fullres"], inplace=True)

    # Load images and scalefactors
    adata.uns["spatial"] = {
        library_id: {
            "images": {
                res: load_spatial_image(os.path.join(sample_adata_path, "spatial", f"tissue_{res}_image.png"), mode=image_color)
                for res in ["hires", "lowres"]
            },
            "scalefactors": json.loads(
                Path(os.path.join(bin_sampled_path, "spatial/scalefactors_json.json")).read_bytes()
            )
        }
    }

    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.X = adata.X.toarray()
    adata.obs["X"] = adata.obsm["spatial"][:, 1]
    adata.obs["Y"] = adata.obsm["spatial"][:, 0]
    adata.obs["sample_id"] = library_id
    return adata

def read_visium_hd(
    adata_path, bin_size=16, genome=None, is_aggr=False, is_xenograft=False, image_color="grayscale"
):
    """Read Visium HD data from a directory, with optional aggregation.

    Args:
        adata_path (str): Path to the root directory of HD data (single sample
            or multi-sample).
        bin_size (int): Bin size for HD samples. Defaults to 16.
        genome (str): Genome for xenograft mapping. Defaults to None.
        is_aggr (bool): If True, aggregate all samples in the directory. Defaults to False.
        is_xenograft (bool): Xenograft data flag. Defaults to False.
        image_color (str): How to load tissue images, one of 'grayscale',
            'original', or 'interpolated'. Defaults to "grayscale".

    Returns:
        anndata.AnnData: Combined annotated data object.
    """
    if is_aggr:
        adatas_concat = None
        for sample_id in os.listdir(adata_path):
            sample_adata_path = os.path.join(adata_path, sample_id)
            adata = read_visium_hd_sample(
                sample_adata_path, bin_size=bin_size, genome=genome,
                library_id=sample_id, is_xenograft=is_xenograft, image_color=image_color
            )
            if adatas_concat is None:
                adatas_concat = adata
            else:
                adatas_concat = ad.concat([adatas_concat, adata], merge="same")
            del adata
            gc.collect()
    else:
        sample_id = os.path.basename(os.path.normpath(adata_path))
        adatas_concat = read_visium_hd_sample(
            adata_path, bin_size=bin_size, genome=genome,
            library_id=sample_id, is_xenograft=is_xenograft, image_color=image_color
        )

    adatas_concat.var_names_make_unique()
    adatas_concat.obs_names_make_unique()
    return adatas_concat

def read_spatial_data(
    adata_path, genome=None, is_xenograft=False,
    is_aggr=True, is_spatial=True, select_sample=None
):
    """Read spatial transcriptomics data and return an AnnData object.

    Args:
        adata_path (str): Directory with spatial transcriptomics data.
        genome (str): Genome reference for xenograft. Defaults to None.
        is_xenograft (bool): Xenograft flag; adds cell type classification data. Defaults to False.
        is_aggr (bool): Whether the data is aggregated across samples. Defaults to True.
        is_spatial (bool): Whether to load spatial data. Defaults to True.
        select_sample (str): Only keep samples whose IDs start with this string. Defaults to None.

    Returns:
        anndata.AnnData: Processed AnnData object.
    """
    sample_data_path = os.path.join(adata_path, "filtered_feature_bc_matrix.h5")
    if os.path.exists(sample_data_path):
        adata = sc.read_10x_h5(sample_data_path, genome=genome if is_xenograft else None)
    else:
        adata = sc.read_10x_mtx(os.path.join(adata_path, "filtered_feature_bc_matrix"))

    # Aggregation support
    if is_aggr:
        adata.obs['sample_id'] = adata.obs.index.to_series().apply(lambda x: int(x.split('-')[-1]) - 1)
        aggregation_map = pd.read_csv(os.path.join(adata_path, 'aggregation.csv'))["library_id"].to_dict()
        adata.obs['sample_id'] = adata.obs['sample_id'].map(aggregation_map)

    if select_sample is not None:
        adata = adata[adata.obs['sample_id'].str.startswith(select_sample)]

    if is_spatial:
        spatial_dir = os.path.join(adata_path, 'spatial')
        if is_aggr:
            adata.uns["spatial"] = {
                library_id: {
                    "images": {res: imread(os.path.join(spatial_dir, library_id, f"tissue_{res}_image.png")) for res in ["hires", "lowres"]},
                    "scalefactors": json.loads(Path(os.path.join(spatial_dir, library_id, "scalefactors_json.json")).read_bytes())
                }
                for library_id in set(adata.obs["sample_id"].values)
                if os.path.isdir(os.path.join(spatial_dir, library_id))
            }
        else:
            adata.uns["spatial"] = {
                "images": {res: imread(os.path.join(spatial_dir, f"tissue_{res}_image.png")) for res in ["hires", "lowres"]},
                "scalefactors": json.loads(Path(os.path.join(spatial_dir, "scalefactors_json.json")).read_bytes())
            }

        # Read tissue positions
        tissue_pos_file = os.path.join(
            adata_path, 'aggr_tissue_positions_list.csv' if is_aggr else 'spatial/tissue_positions_list.csv'
        )
        positions = pd.read_csv(tissue_pos_file, header=None, index_col=0)
        positions.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_col_in_fullres', 'pxl_row_in_fullres']
        adata.obs = adata.obs.join(positions, how="left")
        adata.obsm['spatial'] = adata.obs[['pxl_row_in_fullres', 'pxl_col_in_fullres']].to_numpy().astype(float)
        adata.obs.drop(columns=['pxl_row_in_fullres', 'pxl_col_in_fullres'], inplace=True)

    if is_xenograft:
        gem_file_path = os.path.join(adata_path, 'analysis', 'gem_classification.csv')
        if os.path.exists(gem_file_path):
            gem_df = pd.read_csv(gem_file_path, index_col=0)
            gem_df["admix"] = gem_df['mm10'] / (gem_df['GRCh38'] + gem_df['mm10'])
            adata.obs = adata.obs.join(gem_df, how="left")
        else:
            print(f"File not found: {gem_file_path}. Ensure the GEM classification file is available.")

    adata.var_names_make_unique()
    adata.X = adata.X.toarray()
    adata.obs["X"] = adata.obsm["spatial"][:, 0]
    adata.obs["Y"] = adata.obsm["spatial"][:, 1]
    return adata

def save_list(l, file):
    """Save a list to a CSV file.

    Args:
        l (list): List to save.
        file (str): Output file path.
    """
    df = pd.DataFrame(l)
    df.to_csv(file, index=False)

def read_gmt_file(gmt_type_file):
    """Read a .gmt file and return it as a DataFrame in gene set format.

    Args:
        gmt_type_file (str): Path to the .gmt file.

    Returns:
        pandas.DataFrame: DataFrame with gene set names as columns.
    """
    data = []
    with open(gmt_type_file, 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            data.append(line)
    df_cell_types = pd.DataFrame(data).T
    df_cell_types.columns = df_cell_types.iloc[0, :]
    df_cell_types = df_cell_types.drop(0)
    return df_cell_types

def save_ranked_genes(result_df, results_file, top_genes=-1):
    """Save ranked genes per program/celltype to a CSV file.

    For each column of ``result_df``, genes are ranked by descending score and
    written as a column of gene names, padded to equal length.

    Args:
        result_df (pandas.DataFrame): DataFrame of gene scores.
        results_file (str): Output CSV file path.
        top_genes (int): If greater than 0, only the top N genes are saved. Defaults to -1.
    """
    ranked_data = {
        col: list(result_df[col].dropna().sort_values(ascending=False).index)
        for col in result_df.columns
    }
    max_len = max(len(lst) for lst in ranked_data.values())
    for col in ranked_data:
        ranked_data[col] += [None] * (max_len - len(ranked_data[col]))
    ranked_indices = pd.DataFrame(ranked_data)
    if top_genes > 0:
        ranked_indices = ranked_indices.head(top_genes)
    ranked_indices.to_csv(results_file, index=False)

def load_experiment_result(results_dir_path, sample_name, exp_name, mode, is_annotated=False):
    """Load an experiment result table by mode and sample.

    Args:
        results_dir_path (str): Directory where results are saved.
        sample_name (str): Sample identifier.
        exp_name (str): Experiment identifier.
        mode (str): Result mode, one of 'spots', 'genes', 'genescores',
            'r_genes', 'r_genescores', 'hvar', or 'annotation'.
        is_annotated (bool): If True, remap columns to annotation celltypes. Defaults to False.

    Returns:
        pandas.DataFrame: Result DataFrame.

    Raises:
        ValueError: If ``mode`` is not a recognized result mode.
        FileNotFoundError: If the corresponding result file does not exist.
    """
    mode_key_map = {
        "spots": "topics_per_spot",
        "genes": "genes_per_topic",
        "genescores": "genescores_per_topic",
        "r_genes": "ranked_genes",
        "r_genescores": "ranked_genescores",
        "hvar": "top_genes",
        "annotation": "annotation",
    }
    if mode not in mode_key_map:
        raise ValueError(f"Invalid mode: {mode}")
    mode_key = mode_key_map[mode]
    header_col = 0 if mode != "hvar" else None

    results_file = os.path.join(results_dir_path, f"{mode_key}_{exp_name}_{sample_name}.csv")
    if not os.path.exists(results_file):
        raise FileNotFoundError(f"{results_file} not found")

    results_df = pd.read_csv(results_file, index_col=0, header=header_col)
    results_df.index = results_df.index.map(str)

    if is_annotated:
        annotation = load_experiment_result(results_dir_path, sample_name, exp_name, mode="annotation")
        results_df = results_df.rename(columns=dict(zip(annotation['program'], annotation['celltype'])))
        # Remove columns that start with exp_name (possible leftovers)
        results_df = results_df.drop(columns=[col for col in results_df if col.startswith(exp_name)], errors='ignore')
    return results_df

def get_ground_truth(adata_spatial, mode="genes"):
    """Retrieve the ground truth matrix (genes x celltypes or spots) for a dataset.

    Args:
        adata_spatial (anndata.AnnData): Annotated spatial dataset.
        mode (str): Mode for the ground truth table, one of 'genes' or 'spots'.
            Defaults to "genes".

    Returns:
        pandas.DataFrame: Ground truth table.
    """
    base_dir = "Z:\\MorrissyLab Dropbox\\Visium_profiling\\benchmark"
    data_dir = os.path.join(base_dir, "data")
    if mode == "spots":
        ground_truth = adata_spatial.uns["ground_truth"]
    elif mode == "genes":
        ground_truth = pd.read_csv(os.path.join(data_dir, adata_spatial.uns["dataset_name"], 'cell_type_gene_df.csv'), index_col=0)
    ground_truth.index = ground_truth.index.map(str)
    return ground_truth
