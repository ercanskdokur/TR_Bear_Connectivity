#!/bin/bash
#SBATCH --job-name=36_37_sens
#SBATCH --partition=mid
#SBATCH --account=users
#SBATCH --time=12:00:00
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

run_one () {
  echo "=== RUN $1 | $(date) ==="
  apptainer exec --pwd "${PROJECT_ROOT}/scripts" \
    --bind "${TB_DATA_ROOT:-${PROJECT_ROOT}/data}:${TB_DATA_ROOT:-${PROJECT_ROOT}/data}" \
    --bind "${PROJECT_ROOT}:${PROJECT_ROOT}" \
    "${SIF}" Rscript "${PROJECT_ROOT}/scripts/$1" \
    && echo "[OK] $1" || echo "[FAIL] $1"
}
run_one 36_c_sensitivity.R
run_one 37_patch_threshold.R
echo "[INFO] Done. date=$(date)"
