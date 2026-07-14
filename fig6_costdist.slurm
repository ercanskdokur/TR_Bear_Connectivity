#!/bin/bash
#SBATCH --job-name=fig6_cd
#SBATCH --partition=short
#SBATCH --account=users
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -uo pipefail
echo "[INFO] host=$(hostname) | date=$(date) | job=${SLURM_JOB_ID:-NA}"
PR=${PROJECT_ROOT}
SIF=${PR}/programs/trbear.sif
module purge; module load apptainer/1.4.1

apptainer exec --pwd "${PR}/scripts" \
  --bind "${TB_DATA_ROOT:-${PROJECT_ROOT}/data}:${TB_DATA_ROOT:-${PROJECT_ROOT}/data}" \
  --bind "${PROJECT_ROOT}:${PROJECT_ROOT}" \
  "${SIF}" Rscript "${PR}/scripts/fig6_costdist.R" \
  && echo "[OK] fig6_costdist.R" || echo "[FAIL] fig6_costdist.R"
echo "[INFO] Done. date=$(date)"
