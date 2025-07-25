#!/bin/bash
#SBATCH --job-name=spotnmf
#SBATCH --output=logs/spotnmf_%j.out  # Log file
#SBATCH --error=logs/spotnmf_%j.err   # Error log file
#SBATCH --time=16:00:00                 # Runtime
#SBATCH --partition=gpu-a100,gpu-v100            # GPU partition
#SBATCH --gres=gpu:1                   # Request 1 GPU
#SBATCH --mem=256G                     # 128GB RAM
#SBATCH --cpus-per-task=16             # 16 CPU cores

# Activate Conda environment
source ~/miniconda/etc/profile.d/conda.sh
conda activate spotnmf

# Read input parameters from command-line arguments
GENOME=$1
K=$2
ADATA_PATH=$3
RESULTS_DIR=$4
SAMPLE_NAME=$5

echo "Running spotnmf for genome=${GENOME}, K=${K}"

## One Step
python cli.py spotnmf --sample_name "${SAMPLE_NAME}_${GENOME}_${K}" --adata_path "$ADATA_PATH" --results_dir "$RESULTS_DIR" --k "$K" --genome "$GENOME" --is_aggr --is_xeno

# Run the Python script with the provided arguments
## Three steps
# python cli.py deconvolve --sample_name "${SAMPLE_NAME}_${GENOME}_${K}" --adata_path "$ADATA_PATH" --results_dir "$RESULTS_DIR" --k "$K" --genome "$GENOME" --is_aggr --is_xeno
# python cli.py annotate --sample_name "${SAMPLE_NAME}_${GENOME}_${K}" --results_dir "$RESULTS_DIR" --k "$K" --genome "$GENOME" --is_aggr --is_xeno
# python cli.py plot --sample_name "${SAMPLE_NAME}_${GENOME}_${K}" --adata_path "$ADATA_PATH" --results_dir "$RESULTS_DIR"  --genome "$GENOME" --is_aggr --is_xeno


echo "Job completed for ${GENOME} with K=${K}"
