#!/bin/bash
#SBATCH --job-name=spotnmf
#SBATCH --time=20:00:00                 # Runtime
#SBATCH --partition=gpu-a100,gpu-v100            # GPU partition
#SBATCH --gres=gpu:1                   # Request 1 GPU
#SBATCH --mem=150G                     # 128GB RAM
#SBATCH --cpus-per-task=16            # 16 CPU cores


# Activate Conda environment
source ~/miniconda/etc/profile.d/conda.sh
conda activate spotnmf

# Input arguments
GENOME=$1
K=$2
ADATA_PATH=$3
RESULTS_DIR=$4
SAMPLE_NAME=$5
lr=$6
eps=$7
h=$8
w=$9


# Full sample name used in output
FULL_SAMPLE_NAME="${SAMPLE_NAME}_${GENOME}_${K}"
LOG_DIR="${RESULTS_DIR}/${FULL_SAMPLE_NAME}/logs"
mkdir -p "$LOG_DIR"

# Redirect output and error logs
exec > >(tee -a "${LOG_DIR}/spotnmf.out")
exec 2> >(tee -a "${LOG_DIR}/spotnmf.err" >&2)


echo "========== Starting spotnmf =========="
echo "Sample: $FULL_SAMPLE_NAME"
echo "Genome: $GENOME | K: $K | LR: $lr | EPS: $eps | H: $h | W: $w"
echo "======================================="

# Set HVG file path based on genome
if [ "$GENOME" == "mm10" ]; then
  HVG_FILE="/home/aly.abdelkareem1/spotnmf/varsha_odg/Mm_all.overdispersed_genes.txt"
elif [ "$GENOME" == "GRCh38" ]; then
  HVG_FILE="/home/aly.abdelkareem1/spotnmf/varsha_odg/Hs_all.overdispersed_genes.txt"
else
  echo "Unknown genome: $GENOME"
  exit 1
fi


python cli.py spotnmf \
  --sample_name "$FULL_SAMPLE_NAME" \
  --adata_path "$ADATA_PATH" \
  --results_dir "$RESULTS_DIR" \
  --k "$K" \
  --genome "$GENOME" \
  --is_aggr \
  --is_xeno \
  --lr "$lr" \
  --eps "$eps" \
  --h "$h" \
  --w "$w" \
  --hvg_file "$HVG_FILE" 
  # --normalize_rows

echo "========== Job completed for $FULL_SAMPLE_NAME =========="
