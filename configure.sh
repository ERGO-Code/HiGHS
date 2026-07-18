#!/usr/bin/env bash
# Adaptive config: enable GPU if CUDA present, otherwise CPU-only. Both build and run.
#
# Usage:
#   ./configure.sh                       # adaptive: probe nvcc for GPU/CPU
#   ./configure.sh --no-cuda             # force CPU-only (even if nvcc present)
#   ./configure.sh --cuda                # force GPU (error if no nvcc)
#   ./configure.sh --build mybuild ...   # specify build dir + pass-through extra cmake args
#
# Environment variables:
#   CUDA_HOME          CUDA toolkit root (default: derived from nvcc path)
#   CMAKE_PREFIX_PATH_EXTRA  extra cmake prefix path (e.g. conda env / venv)
set -euo pipefail
cd "$(dirname "$0")"

# ---- Parse args ----
FORCE_GPU=""        # "" adaptive / "on" force on / "off" force off
BUILD_DIR="build"
EXTRA_CMAKE_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --no-cuda) FORCE_GPU="off"; shift;;
    --cuda)    FORCE_GPU="on";  shift;;
    --build)   BUILD_DIR="$2"; shift 2;;
    --build=*) BUILD_DIR="${1#--build=}"; shift;;
    *)         EXTRA_CMAKE_ARGS+=("$1"); shift;;
  esac
done

# ---- 1. GPU adaptive probe ----
NVCC_PATH=""
if [[ "$FORCE_GPU" != "off" ]]; then
  if command -v nvcc >/dev/null 2>&1; then
    NVCC_PATH="$(command -v nvcc)"
  elif [[ -n "${CUDA_HOME:-}" && -x "${CUDA_HOME}/bin/nvcc" ]]; then
    NVCC_PATH="${CUDA_HOME}/bin/nvcc"
  fi
fi

if [[ "$FORCE_GPU" == "on" && -z "$NVCC_PATH" ]]; then
  echo "[configure] error: --cuda specified but nvcc not found" >&2
  exit 1
fi

if [[ -n "${NVCC_PATH}" ]]; then
  echo "[configure] nvcc found: ${NVCC_PATH}"
  "${NVCC_PATH}" --version | tail -2
  EXTRA_CMAKE_ARGS+=(-DCUPDLP_GPU=ON)
  EXTRA_CMAKE_ARGS+=(-DCMAKE_CUDA_ARCHITECTURES=native)
  if [[ -z "${CUDA_HOME:-}" ]]; then
    CUDA_HOME="$(dirname "$(dirname "${NVCC_PATH}")")"
    export CUDA_HOME
  fi
  echo "[configure] CUDA_HOME=${CUDA_HOME} → CUPDLP_GPU=ON"
elif [[ "$FORCE_GPU" == "off" ]]; then
  echo "[configure] --no-cuda specified → CUPDLP_GPU=OFF (forced CPU-only)"
  EXTRA_CMAKE_ARGS+=(-DCUPDLP_GPU=OFF)
else
  echo "[configure] nvcc not found → CUPDLP_GPU=OFF (CPU-only, PDLP still runs on CPU)"
  EXTRA_CMAKE_ARGS+=(-DCUPDLP_GPU=OFF)
fi

# ---- 2. Force-enable gRPC server subproject ----
EXTRA_CMAKE_ARGS+=(-DHIGHS_BUILD_GRPC_SERVER=ON)
EXTRA_CMAKE_ARGS+=(-DCMAKE_BUILD_TYPE=Release)

# ---- 3. Adaptive CMAKE_PREFIX_PATH (conda env / venv / custom) ----
PREFIX_PATHS=()
[[ -n "${CONDA_PREFIX:-}" ]] && PREFIX_PATHS+=("${CONDA_PREFIX}")
[[ -n "${VIRTUAL_ENV:-}" ]]  && PREFIX_PATHS+=("${VIRTUAL_ENV}")
[[ -n "${CMAKE_PREFIX_PATH_EXTRA:-}" ]] && PREFIX_PATHS+=("${CMAKE_PREFIX_PATH_EXTRA}")
if [[ ${#PREFIX_PATHS[@]} -gt 0 ]]; then
  _ifs="${IFS}"; IFS=';'
  _joined="${PREFIX_PATHS[*]}"
  IFS="${_ifs}"
  EXTRA_CMAKE_ARGS+=(-DCMAKE_PREFIX_PATH="${_joined}")
  echo "[configure] CMAKE_PREFIX_PATH=${_joined}"
fi

echo "[configure] cmake -S. -B${BUILD_DIR} ${EXTRA_CMAKE_ARGS[*]}"
cmake -S. -B"${BUILD_DIR}" "${EXTRA_CMAKE_ARGS[@]}"
echo "[configure] done. Next: cmake --build ${BUILD_DIR} --parallel --target highs_grpc_server"
