#!/usr/bin/env bash
set -euo pipefail

ARCH=$(uname -m)
OS=$OSTYPE
APP_NAME="napari-dea-${OS}-${ARCH}"
ENV_DIR="${APP_NAME}/env"
BIN_DIR="${APP_NAME}/bin"
PACK_DIR="env.tar.gz"


echo "[1/5] Cleaning previous build..."
rm -rf "${APP_NAME}" "${BIN_DIR}" "env.tar.gz" "${APP_NAME}.tar.gz" "tmp_env"

echo "[1/5] Creating directories..."
mkdir -p "${APP_NAME}" "${BIN_DIR}"

echo "[2/5] Building the conda environment (prefix)..."
mamba create -y -p ./tmp_env -c conda-forge \
  python=3.11

echo "[3/5] Installing spatialdata and your plugin..."
./tmp_env/bin/pip install \
  "spatialdata >= 0.6" \
  "napari-spatialdata[all,bermuda] >= 0.6.0" \
  "spatialdata-io >= 0.5" \
  "spatialdata-plot >= 0.2.13" \
  git+https://github.com/Yasas1994/napari-spatialdata-dea.git@main

echo "[4/5] Packing environment for redistribution..."
conda-pack -p ./tmp_env -o env.tar.gz

echo "[4.1/5] Extracting packed env into ${ENV_DIR}..."
mkdir -p "${ENV_DIR}"
tar -xzf env.tar.gz -C "${ENV_DIR}"

# Patch absolute paths inside env
echo "[4.2/5] Fixing environment relocation..."
"${ENV_DIR}/bin/conda-unpack"

echo "[4.5/5] Creating launcher..."
cat >"${BIN_DIR}/napari-dea" <<'EOF'
#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ENV_DIR="${APP_DIR}/env"

# Activate portable env
source "${ENV_DIR}/bin/activate"

exec python -m napari "$@"
EOF
chmod +x "${BIN_DIR}/napari-dea"

echo "[5/5] Creating README..."
cat >"${APP_NAME}/README.txt" <<EOF
napari-dea for ${OS} ${ARCH} $(date)

This directory contains:
- A relocatable conda environment in ./env
- A launcher script in ./bin/napari-dea

To launch napari:
    cd ${APP_NAME}
    chmod + ./bin/napari-dea
    ./bin/napari-dea <path to .zarr archive>

    You do **NOT** need conda on your system to run this :).
EOF

echo "[5.1/5] Creating distributable tarball..."
tar -czf "${APP_NAME}.tar.gz" "${APP_NAME}"

echo "[5.1/5] Removing tmp files..."
rm -r "tmp_env" "env.tar.gz"

echo "Done!"
echo "Output created: ${APP_NAME}.tar.gz"