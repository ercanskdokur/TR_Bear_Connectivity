#!/bin/bash
#SBATCH --job-name=trbear_02
#SBATCH --partition=short
#SBATCH --account=users
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --exclude=dy02,dy03
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
echo "[INFO] host=$(hostname) | date=$(date)"

PROJECT_ROOT="${PROJECT_ROOT:-/path/to/tr_bear_connect}"
SIF=${PROJECT_ROOT}/programs/trbear.sif
SCRIPT=${PROJECT_ROOT}/scripts/02_predictors_present_to_tif.R

module purge
module load apptainer/1.4.1

apptainer exec \
  --pwd "${PROJECT_ROOT}/scripts" \
  --bind "${TB_DATA_ROOT:-${PROJECT_ROOT}/data}:${TB_DATA_ROOT:-${PROJECT_ROOT}/data}" \
  --bind "${PROJECT_ROOT}:${PROJECT_ROOT}" \
  "${SIF}" Rscript "${SCRIPT}"

echo "[INFO] Done. date=$(date)"
