import os
import gc
import math

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import distinctipy
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm

from .io import check_dir, load_experiment_result

# --------- Helper variables and sample orders ---------
printable_rotate_dict = {
    'BT143_x1_28d_i_0um_ctrl': 'x',
    'BT143_x1_28d_i_30um_ctrl': 'x',
    'BT143_x3_76d_i_0um_ctrl': '0',
    'BT143_x3_76d_i_40um_ctrl': 'xy',
    'BT143_x2_48d_i_0um_ctrl': 'x',
    'BT143_x2_48d_i_30um_ctrl': 'xy',
    'BT143_x4_76d_i_0um_ctrl': 'y',
    'BT143_x4_76d_i_30um_ctrl': 'xy'
}


samples_bt53 = [
    "BT53_m1_48h_100um_ctrl",
    "BT53_m1_48h_0um_RT",
    "BT53_m1_48h_100um_TMZ",
    "BT53_m1_48h_0um_RT-TMZ",
    "BT53_m2_48h_0um_ctrl",
    "BT53_m2_48h_0um_RT",
    "BT53_m2_48h_0um_TMZ",
    "BT53_m2_48h_0um_RT-TMZ",
    "BT53_m1_2w_0um_ctrl",
    "BT53_m1_2w_0um_RT",
    "BT53_m1_2w_100um_TMZ",
    "BT53_m1_2w_0um_RT-TMZ",
    "BT53_m2_2w_0um_ctrl",
    "BT53_m2_2w_0um_RT",
    "BT53_m2_2w_0um_TMZ",
    "BT53_m2_2w_0um_RT-TMZ",
    "NRG_m1_48h_0um_ctrl",
    "NRG_m1_48h_0um_RT",
    "NRG_m1_48h_0um_TMZ",
    "NRG_m1_48h_0um_RT-TMZ",
    "NRG_m2_48h_0um_ctrl",
    "NRG_m2_48h_0um_RT",
    "NRG_m2_48h_0um_TMZ",
    "NRG_m2_48h_0um_RT-TMZ",
    "NRG_m1_2w_0um_ctrl",
    "NRG_m1_2w_0um_RT",
    "NRG_m1_2w_0um_TMZ",
    "NRG_m1_2w_0um_RT-TMZ",
    "NRG_m2_2w_0um_ctrl",
    "NRG_m2_2w_0um_RT",
    "NRG_m2_2w_0um_TMZ",
    "NRG_m2_2w_0um_RT-TMZ"
]

# Function to create a dataframe of (x, y) locations for a grid of area X
def generate_random_locations(N):
    def get_square_dimensions(X):
        sqrt_X = math.sqrt(X)
        x = math.floor(sqrt_X)
        y = math.ceil(X / x)
        return x, y

    x, y = get_square_dimensions(N)
    locations = [(i % x, i // x) for i in range(N)]
    df = pd.DataFrame(locations, columns=['X', 'Y'])
    return df

# Define a custom colormap for spatial plotting
spatial_color_map = LinearSegmentedColormap.from_list("", ['cornsilk', 'magenta', 'navy'])


def plot_df_heatmap(df, title_name, x_name, y_name, results_dir_path, cmap = "Blues", is_cluster = False):
    """
    Plot a DataFrame as a heatmap or clustered heatmap and save as PDF.
    """
    fig_size = [10 + df.shape[1]/4, 10 + df.shape[0]/6]
    if(is_cluster):
        fig_cluster = sns.clustermap(df.astype(float).fillna(0), cmap=cmap, col_cluster=True, row_cluster=True, figsize=fig_size, xticklabels=True, yticklabels=True)
        fig_cluster.ax_heatmap.set_title(title_name)
        fig_cluster.ax_heatmap.set_xlabel(x_name)
        fig_cluster.ax_heatmap.set_ylabel(y_name)
        fig_cluster.ax_heatmap.tick_params(axis='y', labelrotation=0, labelright=True, labelleft=False)
        fig_cluster.ax_heatmap.tick_params(axis='x', labelrotation=90)
        for _, spine in fig_cluster.ax_heatmap.spines.items():
            spine.set_visible(True)
            spine.set_color('#aaaaaa')
        fig_cluster.cax.set_visible(False)
        fig_cluster.savefig(os.path.join(results_dir_path, f"clustermap_{title_name}.pdf") )
    else:
        fig, ax_plot = plt.subplots(figsize=fig_size, layout="constrained")
        sns.heatmap(df.astype(float), cmap=cmap, ax=ax_plot, xticklabels=True, yticklabels=True)
        ax_plot.set_title(title_name)
        ax_plot.set_xlabel(x_name)
        ax_plot.set_ylabel(y_name)
        ax_plot.tick_params(axis='y', labelrotation=0, labelright=True, labelleft=False)
        ax_plot.tick_params(axis='x', labelrotation=90)
        for _, spine in ax_plot.spines.items():
            spine.set_visible(True)
            spine.set_color('#aaaaaa')

        fig.savefig(os.path.join(results_dir_path, f"heatmap_{title_name}.pdf") )
    

def plot_spatial_topic(adata_spatial, topic, axe, use_scanpy=False):
    """
    Plot a single spatial topic on a given matplotlib axis.

    Parameters
    ----------
    adata_spatial : AnnData
        AnnData object with spatial info.
    topic : str
        Column name in adata_spatial.obs to plot.
    axe : matplotlib.axes.Axes
        Axis to plot on.
    use_scanpy : bool
        Use scanpy.pl.spatial (True) or direct scatter (False).
    """
    if use_scanpy:
        # Use Scanpy's built-in spatial plotting function
        scatter = sc.pl.spatial(
            adata_spatial,
            basis="spatial",
            color=topic,
            spot_size=150,
            color_map=spatial_color_map,
            ax=axe,
            show=False,
            return_fig=True
        )
    else:
        # Custom scatter plot for spatial data
        scatter = axe.scatter(
            x=adata_spatial.obs["X"],
            y=adata_spatial.obs["Y"],
            c=adata_spatial.obs[topic].values,
            cmap=spatial_color_map,
            vmin=0,
            vmax=1,
            s=10,
            marker='o'
        )

    axe.set_title(f"{topic}", fontsize=6)
    return scatter

def plot_spatial_all_topics(adata_spatial, rf_usages, results_dir_path, title_name="Spatial Topics Plot", fig_width = None, fig_height = None, is_show=False):
    """
    Plot multiple spatial topics on a grid of subplots and save the figure as a PDF.

    Parameters:
    - adata_spatial: AnnData object with spatial coordinates in .obs.
    - rf_usages: DataFrame, each column is a topic to plot, normalized within the function.
    - results_dir_path: str, directory to save the output PDF.
    - title_name: str, title for the entire figure.
    - is_show: bool, if True displays the plot; otherwise, only saves it.
    """
    params = adata_spatial.uns.get("params", {})

    plt.ioff()  # Turn off interactive plotting for faster processing

    # Normalize the rf_usages data
    # rf_usages = (rf_usages - np.min(rf_usages)) / (np.max(rf_usages) or 1)
    topics_list = rf_usages.columns

    # Join topics data with spatial AnnData and fill any NaN values with 0
    adata_spatial.obs = adata_spatial.obs.join(rf_usages).fillna(0)

    n_topics = len(topics_list)

    if params.get("generate_locations", False):
        df_locations =  generate_random_locations(N=len(adata_spatial.obs))
        adata_spatial.obs["X"] = df_locations["X"].values
        adata_spatial.obs["Y"] = df_locations["Y"].values

    # Set figure dimensions if not provided in params
    fig_width = params.get("fig_width", None)
    fig_height = params.get("fig_height", None)


    # width, height = len(adata_spatial.obs["X"].unique()), len(adata_spatial.obs["Y"].unique())
    # fig_width_single, fig_height_single = max(width / 5, 2), max(height / 5, 2)
    fig_width_single, fig_height_single = 8.5, 11
    
    ncols = int(math.ceil(math.sqrt(n_topics)))
    nrows = int(math.ceil(n_topics / ncols))

    fig_width = fig_width or (fig_width_single * ncols)
    fig_height = fig_height or (fig_height_single * nrows)

    # Create a figure and a grid of axes
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(fig_width, fig_height), constrained_layout=True)
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    axes = axes.flatten() if n_topics > 1 else [axes]

    # Plot each topic on its respective axis
    for topic, axe in zip(topics_list, axes[:n_topics]):
        scatter = plot_spatial_topic(adata_spatial, topic=topic, axe=axe, use_scanpy=False)

    # Turn off axes
    for axe in axes:
        axe.axis('off')

    cbar = fig.colorbar(scatter, ax=axes, location='right', fraction=0.01, pad=0.05)
    cbar.set_label('Intensity', rotation=270, labelpad=10)

    fig.suptitle(f'{title_name}', fontsize=10)

    # Save the figure
    os.makedirs(results_dir_path, exist_ok=True)
    plt.savefig(os.path.join(results_dir_path, f"topics_plot_{title_name}.pdf"), dpi=300)
    if is_show:
        plt.show()
    plt.close('all')


def plot_spatial_all_topics_aggr_manuscript(adata, rf_usages, results_dir_path, title_name="Spatial Topics Plot",
                                            same_legend = False, plot_lots=True, COLS=None, ROWS=None, filter_th=1):
    os.makedirs(results_dir_path, exist_ok=True)

    if(filter_th < 1):
        print(f"Filtered to {filter_th} percentile")
        percentiles = rf_usages.quantile(filter_th)
        rf_usages = rf_usages.apply(lambda col: col.where(col >= percentiles[col.name], 0))
        title_name += f"_filtered_{filter_th}_compressed"

    output_file = os.path.join(results_dir_path, f"topics_plot_{title_name}_printable.pdf")

    common_cols = adata.obs.columns.intersection(rf_usages.columns)
    adata.obs.drop(columns=common_cols, inplace=True)
    adata.obs = adata.obs.join(rf_usages.fillna(0), how="left").fillna(0)
    samples_list = adata.obs['sample_id'].unique()
    color_map = LinearSegmentedColormap.from_list("",['cornsilk','magenta','navy'])

    plt.ioff()
    with PdfPages(output_file) as pdf:
        COLS = COLS or math.ceil(math.sqrt(len(samples_list)))
        ROWS = ROWS or math.ceil(len(samples_list) / COLS)
        max_rows = ROWS*COLS

        for topic in tqdm(rf_usages.columns):
            if (filter_th < 1):
                adata_spatial = adata[adata.obs[topic] > 0,:].copy()
                


            if(max_rows == ROWS*COLS): 
                figure, axes_hist = plt.subplots(nrows=ROWS, ncols=COLS, figsize=(8.5, 11))
                if isinstance(axes_hist, np.ndarray):
                    axes_hist = axes_hist.flatten().tolist()
                else:
                    axes_hist = [axes_hist]
    
                iter_ax = 0
            

            top_third_value = sorted(list(adata_spatial.obs[topic].values), reverse = True)[2]
            VMIN = 0
            VMAX = 0.0001 if top_third_value == 0 else top_third_value
            for sample_plot in samples_list:
                selected_sample = adata_spatial[adata_spatial.obs["sample_id"]== sample_plot].copy()

                _ = sc.pl.spatial(selected_sample, img_key="hires", color_map=color_map,
                                    legend_loc=None, colorbar_loc=None,
                                    vmin = VMIN, vmax = VMAX, library_id = sample_plot, color=topic, size=1,
                                    show=False, return_fig=False, save= False, ax = axes_hist[iter_ax]) 
                
                axes_hist[iter_ax].get_xaxis().set_visible(False)
                axes_hist[iter_ax].get_yaxis().set_visible(False)
                axes_hist[iter_ax].set_title(f"{topic}", fontsize=4, pad=-10)
                axes_hist[iter_ax].set_aspect('equal', adjustable='box')  # makes axis square
                
                iter_ax+=1
                max_rows-=1
                del selected_sample

            # Adjust layout to fit plots and colorbar nicely
            sm = cm.ScalarMappable(cmap=plt.get_cmap(color_map), norm=plt.Normalize(vmin=VMIN, vmax=VMAX))
            sm.set_array([])  # Required for colorbar
            cbar = figure.colorbar(sm, ax=axes_hist[iter_ax-len(samples_list):iter_ax], orientation='vertical', pad=0.01)
            cbar.set_ticks([VMIN, VMAX])
            cbar.ax.set_yticklabels([f"{VMIN:.2f}", f"{VMAX:.1g}"], fontsize=4)
    
            if(max_rows == 0):
                for ax in axes_hist:
                    for spine in ax.spines.values():
                        spine.set_visible(False)
                    ax.xaxis.set_ticks_position('none')
                    ax.yaxis.set_ticks_position('none')
                    ax.set_xticklabels([])
                    ax.set_yticklabels([])
                pdf.savefig()
                plt.close(figure)
                plt.close('all')
                gc.collect()  # Collect garbage after each iteration to free memory
                max_rows = ROWS*COLS


        # Save the last figure if it exists and was not saved
        # if(max_rows !=0):
        if 'figure' in locals():
            for ax in axes_hist:
                for spine in ax.spines.values():
                    spine.set_visible(False)
                ax.xaxis.set_ticks_position('none')
                ax.yaxis.set_ticks_position('none')
                ax.set_xticklabels([])
                ax.set_yticklabels([])
            pdf.savefig(figure)
            plt.close(figure)
            plt.close('all')
            gc.collect()

    print(f"Plot {title_name} is completed")



def plot_spatial_all_topics_aggr(adata, rf_usages, results_dir_path, title_name="Spatial Topics Plot",
                                  same_legend=False, plot_topic=True, COLS=None, ROWS=None, image_key=None, filter_th = 1.0):
    
    os.makedirs(results_dir_path, exist_ok=True)

    if(filter_th < 1):
        print(f"Filtered to {filter_th} percentile")
        percentiles = rf_usages.quantile(filter_th)
        rf_usages = rf_usages.apply(lambda col: col.where(col >= percentiles[col.name], 0))
        title_name += f"_filtered_{filter_th}"


    output_file = os.path.join(results_dir_path, f"topics_plot_{title_name}.pdf")

    # Optimize merging step
    rf_usages = rf_usages.fillna(0)
    adata.obs = adata.obs.drop(columns=adata.obs.columns.intersection(rf_usages.columns), errors='ignore')
    adata.obs = adata.obs.join(rf_usages, how="left").fillna(0)

    color_map = LinearSegmentedColormap.from_list("", ['cornsilk', 'magenta', 'navy'])
    samples_list = adata.obs['sample_id'].unique()

    if "BT53" in results_dir_path:
        samples_list = sorted(samples_list, key=lambda x: samples_bt53.index(x) if x in samples_bt53 else len(samples_bt53))
        COLS = 8
    elif "pdx_merge" in results_dir_path:
        COLS = 4

    COLS = COLS or math.ceil(math.sqrt(len(samples_list)))
    ROWS = ROWS or math.ceil(len(samples_list) / COLS)
    print(f"plot cols={COLS} rows={ROWS}")
    plt.ioff()

    with PdfPages(output_file) as pdf:
        for topic in tqdm(rf_usages.columns, desc="Plotting topics"):
            topic_max_value = None
            if same_legend:
                topic_vals = adata.obs[topic].values
                if np.count_nonzero(topic_vals) >= 3:
                    topic_max_value = np.partition(topic_vals, -3)[-3]
                else:
                    topic_max_value = np.max(topic_vals)

            fig, axes = plt.subplots(nrows=ROWS, ncols=COLS, figsize=(COLS * 4, ROWS * 4))
            if isinstance(axes, np.ndarray):
                axes = axes.flatten().tolist()
            else:
                axes = [axes]

            for idx, sample_id in enumerate(samples_list):
                if idx >= len(axes): break  # safety

                ax = axes[idx]
                adata_sub = adata[adata.obs['sample_id'] == sample_id].copy()

                if plot_topic:
                    values = adata_sub.obs[topic].values
                    top_val = topic_max_value if same_legend else np.partition(values, -3)[-3] if len(values) >= 3 else np.max(values)
                    vmax = 0.0001 if top_val == 0 else top_val

                    sc.pl.spatial(adata_sub, img_key=image_key, color_map=color_map, legend_loc=None,
                                  colorbar_loc=None, library_id=sample_id, color=topic, size=1.3,
                                  vmin=0, vmax=vmax, show=False, return_fig=False, save=False, ax=ax)

                    cbar = plt.colorbar(cm.ScalarMappable(cmap=color_map), ax=ax, ticks=[0, 1], shrink=0.8, pad=0.03)
                    cbar.ax.set_yticklabels(["0.00", f"{vmax:.1g}"], fontsize=10)
                else:
                    sc.pl.spatial(adata_sub, img_key=image_key, library_id=sample_id, color=topic, size=1.3,
                                  show=False, return_fig=False, save=False, ax=ax)

                ax.set_title(f"{sample_id}_{topic}", fontsize=10)
                ax.axis("off")

                rotation = printable_rotate_dict.get(sample_id, "")
                if 'x' in rotation:
                    ax.invert_xaxis()
                if 'y' in rotation:
                    ax.invert_yaxis()

                del adata_sub

            for ax in axes[len(samples_list):]:
                ax.axis("off")

            fig.suptitle(f"{title_name}_{topic}")
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            gc.collect()

def plot_benchmark_methods_topics(results_dir, adata_spatial, experiments_list, fig_width = None, fig_height = None, use_scanpy=False, is_show=False):

    plt.ioff()  # Turn off interactive plotting for faster processing

    sample_name = adata_spatial.uns["dataset_name"]
    
    params = adata_spatial.uns["params"]

    results_dir_path = check_dir(os.path.join(results_dir, sample_name, "analysis"))
    print(results_dir_path)

    df_gt = adata_spatial.uns["ground_truth"]
    normalized_gt = df_gt.div(df_gt.sum(axis=1), axis=0)
    cell_types = df_gt.columns

    n_cell_types = len(cell_types)
    n_methods = len(experiments_list)

    if params.get("generate_locations"):
        df_locations =  generate_random_locations(N=len(adata_spatial.obs))
        adata_spatial.obs["X"] = df_locations["X"].values
        adata_spatial.obs["Y"] = df_locations["Y"].values

    width, height = len(adata_spatial.obs["X"].unique()), len(adata_spatial.obs["Y"].unique())
    fig_width_single, fig_height_single = max(width / 5, 2), max(height / 5, 2)

    fig_width = fig_width or min((fig_width_single * (n_methods + 1)), 50)
    fig_height = fig_height or min((fig_height_single * n_cell_types), 50)
    print(fig_width, fig_height)

    fig, axes = plt.subplots(ncols=n_methods + 1, nrows=n_cell_types, figsize=(fig_width, fig_height), constrained_layout=True)
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    axes = axes.flatten() if n_cell_types > 1 else [axes]

    plot_index = 0
    for cell_type in cell_types:
        adata_spatial.obs[f'GT_{cell_type}'] = normalized_gt[cell_type].values
        scatter = plot_spatial_topic(adata_spatial, topic=f'GT_{cell_type}', axe=axes[plot_index], use_scanpy=use_scanpy)
        plot_index+=1

        for exp_name in experiments_list:    
            spots_df = load_experiment_result(results_dir_path=os.path.join(results_dir, sample_name, f"{exp_name}_{sample_name}"), sample_name=sample_name, exp_name=exp_name, mode="spots", is_annotated=True)
            # spots_df = (spots_df - spots_df.min()) / (spots_df.max() - spots_df.min())
            adata_spatial.obs[f'{exp_name}_{cell_type}'] = spots_df[cell_type]
            scatter = plot_spatial_topic(adata_spatial, topic=f'{exp_name}_{cell_type}', axe=axes[plot_index], use_scanpy=use_scanpy)
            plot_index+=1

    for axe in axes:
        axe.axis('off')

    cbar = fig.colorbar(scatter, ax=axes, location='right', fraction=0.01, pad=0.05)
    cbar.set_label('Intensity', rotation=270, labelpad=10)

    plt.savefig(os.path.join(results_dir_path, f"methods_topics_plot_{sample_name}_.pdf"), dpi=300)
    if is_show:
        plt.show()
    plt.close('all')

