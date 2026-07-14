#!/bin/bash
#SBATCH --job-name=14_unicor
#SBATCH --partition=long
#SBATCH --account=users
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --array=0-6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --exclude=dy02,dy03
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

SCENARIOS=(present 2041_2070_ssp126 2041_2070_ssp370 2041_2070_ssp585 \
           2071_2100_ssp126 2071_2100_ssp370 2071_2100_ssp585)
SCEN="${SCENARIOS[$SLURM_ARRAY_TASK_ID]}"

echo "[INFO] host=$(hostname) | date=$(date) | job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} | scen=${SCEN}"

PROJECT_ROOT="${PROJECT_ROOT:-/path/to/tr_bear_connect}"
SIF=${PROJECT_ROOT}/programs/trbear.sif
UNICOR_SRC=${PROJECT_ROOT}/programs/UNICOR/unicor
UTIL_SRC=${PROJECT_ROOT}/programs/UNICOR/utilities
SCEN_DIR=${PROJECT_ROOT}/outputs/derived/unicor/${SCEN}

if [[ ! -d "${SCEN_DIR}" ]]; then
  echo "[ERROR] scenario dir missing: ${SCEN_DIR}"; exit 1
fi

# Per-scenario UNICOR working tree so ../utilities/ resolves correctly
# and the 7 array tasks don't collide with each other.
WORK=${SCEN_DIR}/unicor
rm -rf "${WORK}"
mkdir -p "${WORK}"
cp -r "${UNICOR_SRC}/." "${WORK}/"
ln -sfn "${UTIL_SRC}" "${SCEN_DIR}/utilities"

# Move scenario inputs into the working tree (basenames in the .rip)
cp -f "${SCEN_DIR}/${SCEN}.rsg" "${WORK}/${SCEN}.rsg"
cp -f "${SCEN_DIR}/${SCEN}.rip" "${WORK}/${SCEN}.rip"
cp -f "${SCEN_DIR}/sources.xy" "${WORK}/sources.xy"

module purge
module load apptainer/1.4.1

cd "${WORK}"
echo "[INFO] CWD=$(pwd) | files:"; ls -la *.rsg *.rip *.xy

apptainer exec \
  --pwd "${WORK}" \
  --bind "${TB_DATA_ROOT:-${PROJECT_ROOT}/data}:${TB_DATA_ROOT:-${PROJECT_ROOT}/data}" \
  --bind "${PROJECT_ROOT}:${PROJECT_ROOT}" \
  "${SIF}" python3 UNICOR.py "${SCEN}.rip"

# Collect outputs into a tidy results/ dir
RESULTS=${SCEN_DIR}/results
mkdir -p "${RESULTS}"
shopt -s nullglob
for pat in "${SCEN}"*.csv "${SCEN}"*.txt "${SCEN}"*.tif "${SCEN}"*.asc "${SCEN}"*.log; do
  if [[ -e "${WORK}/${pat}" ]]; then
    mv "${WORK}/${pat}" "${RESULTS}/"
  fi
done
# Also catch the pickle / shelf files UNICOR writes
for f in "${WORK}/${SCEN}"* ; do
  [[ -e "$f" ]] || continue
  bn=$(basename "$f")
  [[ "$bn" == "${SCEN}.rsg" || "$bn" == "${SCEN}.rip" ]] && continue
  mv "$f" "${RESULTS}/" 2>/dev/null || true
done

echo "[INFO] Results in: ${RESULTS}"
ls -la "${RESULTS}/" | head -40

echo "[INFO] Done. date=$(date)"
