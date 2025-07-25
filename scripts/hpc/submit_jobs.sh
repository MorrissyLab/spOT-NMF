#!/bin/bash

# Define paths and other parameters
ADATA_PATH="/work/morrissy_lab/spatial_data/samples_aggregate/pdx_merge_all/outs/"
RESULTS_DIR="/work/morrissy_lab/spatial_data/wassgard_xenograft"


sbatch spotnmf_job.sh mm10 90 "$ADATA_PATH" "$RESULTS_DIR" cell2location
sbatch spotnmf_xenograft.sh mm10 90 "$ADATA_PATH" "$RESULTS_DIR" pdx_merge_all 0.001 0.05 0.01 0.01
sbatch spotnmf_visiumhd.sh GRCh38 50 "$ADATA_PATH" "$RESULTS_DIR" P5_CRC 8 0.001 0.05 0.01 0.01

echo "All jobs submitted!"