#!/bin/bash
#SBATCH --job-name=08_present
#SBATCH --partition=short
#SBATCH --account=users
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
echo "[INFO] host=$(hostname) | date=$(date) | job=$SLURM_JOB_ID"

PROJECT_ROOT="${PROJECT_ROOT:-/path/to/tr_bear_connect}"
SIF=${PROJECT_ROOT}/programs/trbear.sif
SCRIPT=${PROJECT_ROOT}/scripts/08_postprocess_present.R

module purge
module load apptainer/1.4.1

apptainer exec \
  --pwd "${PROJECT_ROOT}/scripts" \
  --bind "${TB_DATA_ROOT:-${PROJECT_ROOT}/data}:${TB_DATA_ROOT:-${PROJECT_ROOT}/data}" \
  --bind "${PROJECT_ROOT}:${PROJECT_ROOT}" \
  "${SIF}" Rscript "${SCRIPT}"

echo "[INFO] Done. date=$(date)"
