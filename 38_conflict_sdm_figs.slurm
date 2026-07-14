#!/bin/bash
#SBATCH --job-name=38_conflict
#SBATCH --partition=short
#SBATCH --account=users
#SBATCH --time=00:40:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
module purge; module load apptainer/1.4.1
PR=${PROJECT_ROOT}
apptainer exec --pwd "$PR/scripts" --bind "${TB_DATA_ROOT:-${PROJECT_ROOT}/data}:${TB_DATA_ROOT:-${PROJECT_ROOT}/data}" \
  --bind "${PROJECT_ROOT}:${PROJECT_ROOT}" \
  "$PR/programs/trbear.sif" Rscript "$PR/scripts/38_conflict_sdm_figs.R" \
  && echo "[OK] 38" || echo "[FAIL] 38"
