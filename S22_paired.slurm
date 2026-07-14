#!/bin/bash
#SBATCH --job-name=S22_paired
#SBATCH --partition=short
#SBATCH --account=users
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -uo pipefail
echo "[INFO] host=$(hostname) | date=$(date) | job=${SLURM_JOB_ID:-NA}"
PROJECT_ROOT="${PROJECT_ROOT:-/path/to/tr_bear_connect}"
SIF=${PROJECT_ROOT}/programs/trbear.sif
module purge; module load apptainer/1.4.1
apptainer exec --pwd "${PROJECT_ROOT}/scripts" \
  --bind "${TB_DATA_ROOT:-${PROJECT_ROOT}/data}:${TB_DATA_ROOT:-${PROJECT_ROOT}/data}" \
  --bind "${PROJECT_ROOT}:${PROJECT_ROOT}" \
  "${SIF}" Rscript "${PROJECT_ROOT}/scripts/S22_paired.R" \
  && echo "[OK] S22_paired.R" || echo "[FAIL] S22_paired.R"
echo "[INFO] Done. date=$(date)"
