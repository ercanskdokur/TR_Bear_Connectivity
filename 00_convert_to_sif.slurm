#!/bin/bash
#SBATCH --job-name=trbear_sif
#SBATCH --partition=mid
#SBATCH --account=users
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --exclude=dy02,dy03
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
echo "[INFO] host=$(hostname) | date=$(date)"

PROJECT_ROOT="${PROJECT_ROOT:-/path/to/tr_bear_connect}"
TAR_GZ=${STAGE_DIR:-${PROJECT_ROOT}}/trbear.tar.gz
TAR_PATH=${STAGE_DIR:-${PROJECT_ROOT}}/trbear.tar
SIF_OUT=${PROJECT_ROOT}/programs/trbear.sif

module purge
module load apptainer/1.4.1
apptainer --version

export APPTAINER_TMPDIR=${APPTAINER_TMPDIR:-${PROJECT_ROOT}/apptainer_tmp}
export APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-${PROJECT_ROOT}/apptainer_cache}
mkdir -p "${APPTAINER_TMPDIR}" "${APPTAINER_CACHEDIR}"

if [ ! -f "${TAR_PATH}" ]; then
    echo "[INFO] decompressing ${TAR_GZ} -> ${TAR_PATH}"
    gunzip -c "${TAR_GZ}" > "${TAR_PATH}"
fi
ls -lh "${TAR_PATH}"

echo "[INFO] Building SIF from docker-archive (overwriting if exists)..."
apptainer build --force "${SIF_OUT}" "docker-archive://${TAR_PATH}"

echo "[INFO] SIF built. Sanity-check ENMTML pipeline deps..."
ls -lh "${SIF_OUT}"
apptainer exec "${SIF_OUT}" R --version | head -2
apptainer exec "${SIF_OUT}" python --version
apptainer exec "${SIF_OUT}" R --no-save -e "
  for (p in c('ENMTML','terra','sf','dismo','gam','kernlab','randomForest',
              'maxnet','gbm','spThin','landscapemetrics','rgdal','rgeos'))
    cat(sprintf('%-20s %s\n', p, as.character(packageVersion(p))))
"
apptainer exec "${SIF_OUT}" python -c "import numpy, rasterio, geopandas, networkx; print('py deps OK')"
apptainer exec "${SIF_OUT}" ls /opt/programs/UNICOR | head -10
apptainer exec "${SIF_OUT}" ls /opt/programs/conefor_cli 2>/dev/null || echo "[WARN] conefor_cli missing"

echo "[INFO] Cleaning up uncompressed tar..."
rm -f "${TAR_PATH}"

echo "[INFO] Done. date=$(date)"
