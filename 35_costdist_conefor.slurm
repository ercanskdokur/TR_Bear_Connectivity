#!/bin/bash
#SBATCH --job-name=35_costdist
#SBATCH --partition=mid
#SBATCH --account=users
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -uo pipefail
echo "[INFO] host=$(hostname) | date=$(date) | job=${SLURM_JOB_ID:-NA}"
PROJECT_ROOT="${PROJECT_ROOT:-/path/to/tr_bear_connect}"
SIF=${PROJECT_ROOT}/programs/trbear.sif
module purge
module load apptainer/1.4.1

apptainer exec --pwd "${PROJECT_ROOT}/scripts" \
  --bind "${TB_DATA_ROOT:-${PROJECT_ROOT}/data}:${TB_DATA_ROOT:-${PROJECT_ROOT}/data}" \
  --bind "${PROJECT_ROOT}:${PROJECT_ROOT}" \
  "${SIF}" Rscript "${PROJECT_ROOT}/scripts/35_costdist_conefor.R" \
  && echo "[OK] 35_costdist_conefor.R" || echo "[FAIL] 35_costdist_conefor.R"
echo "[INFO] Done. date=$(date)"
