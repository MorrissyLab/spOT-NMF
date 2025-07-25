#!/bin/bash
#SBATCH --job-name=spotnmf
#SBATCH --time=23:30:00
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --mem=256G
#SBATCH --cpus-per-task=20

# Activate Conda environment
source ~/miniconda/etc/profile.d/conda.sh
conda activate spotnmf

# Input arguments
GENOME=$1
K=$2
ADATA_PATH=$3
RESULTS_DIR=$4
SAMPLE_NAME=$5
BIN=$6
LR=$7
EPS=$8
H=$9
W=${10}

# Full sample name used in output
FULL_SAMPLE_NAME="${SAMPLE_NAME}_${GENOME}_${K}_${BIN}um"
LOG_DIR="${RESULTS_DIR}/${FULL_SAMPLE_NAME}/logs"
mkdir -p "$LOG_DIR"

# Redirect output and error logs
exec > >(tee -a "${LOG_DIR}/spotnmf.out")
exec 2> >(tee -a "${LOG_DIR}/spotnmf.err" >&2)

echo "========== Starting spotnmf =========="
echo "Sample: $FULL_SAMPLE_NAME"
echo "Genome: $GENOME | K: $K | Bin: ${BIN}um | LR: $LR | EPS: $EPS | H: $H | W: $W"
echo "======================================="

# Run the main script
python cli.py spotnmf \
  --sample_name "$FULL_SAMPLE_NAME" \
  --adata_path "$ADATA_PATH" \
  --results_dir "$RESULTS_DIR" \
  --k "$K" \
  --genome "$GENOME" \
  --data_mode "visium_hd" \
  --bin_size "$BIN" \
  --lr "$LR" \
  --eps "$EPS" \
  --h "$H" \
  --w "$W"  # --is_aggr  # --is_xeno \


echo "========== Job completed for $FULL_SAMPLE_NAME =========="
