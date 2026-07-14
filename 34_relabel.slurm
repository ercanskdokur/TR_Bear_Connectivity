#!/bin/bash
#SBATCH --job-name=34_relabel
#SBATCH --partition=short
#SBATCH --account=users
#SBATCH --time=01:30:00
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

run_one 34_core_crosswalk.R
run_one 27_conefor_components.R
run_one 28_dpc_rankshift.R
run_one 30_pa_network_gap.R

echo "[INFO] Done. date=$(date)"
