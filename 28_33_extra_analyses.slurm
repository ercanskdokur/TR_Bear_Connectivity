#!/bin/bash
#SBATCH --job-name=28_33_extra
#SBATCH --partition=short
#SBATCH --account=users
#SBATCH --time=02:00:00
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
  local script="$1"
  echo "==================================================================="
  echo "[INFO] RUN ${script} | $(date)"
  echo "==================================================================="
  apptainer exec \
    --pwd "${PROJECT_ROOT}/scripts" \
    --bind "${TB_DATA_ROOT:-${PROJECT_ROOT}/data}:${TB_DATA_ROOT:-${PROJECT_ROOT}/data}" \
    --bind "${PROJECT_ROOT}:${PROJECT_ROOT}" \
    "${SIF}" Rscript "${PROJECT_ROOT}/scripts/${script}" \
    && echo "[OK] ${script}" || echo "[FAIL] ${script} (continuing)"
}

run_one 28_dpc_rankshift.R
run_one 29_network_robustness.R
run_one 30_pa_network_gap.R
run_one 31_corridor_centroid_shift.R
run_one 32_corridor_validation.R
run_one 33_pinch_road_priority.R

echo "[INFO] All done. date=$(date)"
