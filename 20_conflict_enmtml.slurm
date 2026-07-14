#!/bin/bash
#SBATCH --job-name=20_conflict
#SBATCH --partition=long
#SBATCH --account=users
#SBATCH --time=48:00:00
#SBATCH --requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --exclude=dy02,dy03
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
echo "[INFO] host=$(hostname) | date=$(date) | job=$SLURM_JOB_ID"

PROJECT_ROOT="${PROJECT_ROOT:-/path/to/tr_bear_connect}"
SIF=${PROJECT_ROOT}/programs/trbear.sif
SCRIPT=${PROJECT_ROOT}/scripts/20_conflict_enmtml.R

module purge
module load apptainer/1.4.1

apptainer exec \
  --pwd "${PROJECT_ROOT}/scripts" \
  --bind "${TB_DATA_ROOT:-${PROJECT_ROOT}/data}:${TB_DATA_ROOT:-${PROJECT_ROOT}/data}" \
  --bind "${PROJECT_ROOT}:${PROJECT_ROOT}" \
  "${SIF}" Rscript "${SCRIPT}"

echo "[INFO] Done. date=$(date)"
