#!/usr/bin/env bash
# Adaptive build script for highs-server (external interface to HiGHS).
#
# Two-step build:
#   Step 1: build & install HiGHS (with CUPDLP_GPU=ON if CUDA available)
#   Step 2: build highs-server against the installed HiGHS
#
# Usage:
#   ./configure.sh                          # adaptive: probe nvcc for GPU/CPU
#   ./configure.sh --no-cuda                # force CPU-only (even if nvcc present)
#   ./configure.sh --cuda                   # force GPU (error if no nvcc)
#   ./configure.sh --build-dir <dir>        # custom build dir (default: build)
#   ./configure.sh --install-prefix <path>  # custom HiGHS install prefix
#   ./configure.sh -- <extra cmake args>    # pass-through to HiGHS cmake
#
# Environment variables:
#   CUDA_HOME          CUDA toolkit root (default: derived from nvcc path)
#   CMAKE_PREFIX_PATH_EXTRA  extra cmake prefix path (e.g. conda env / venv)
set -euo pipefail
cd "$(dirname "$0")"

# ---- Parse args ----
FORCE_GPU=""
HIGHS_BUILD_DIR="build-highs"
SERVER_BUILD_DIR="build-server"
INSTALL_PREFIX="${PWD}/build-highs-install"
EXTRA_CMAKE_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --no-cuda) FORCE_GPU="off"; shift;;
    --cuda)    FORCE_GPU="on";  shift;;
    --build-dir)        SERVER_BUILD_DIR="$2"; shift 2;;
    --install-prefix)   INSTALL_PREFIX="$2"; shift 2;;
    --) shift; EXTRA_CMAKE_ARGS+=("$@"); break;;
    *) EXTRA_CMAKE_ARGS+=("$1"); shift;;
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

HIGHS_CMAKE_ARGS=()
# Use system gcc-11/g++-11 to avoid conda gcc-12 libstdc++ ABI mismatch
if command -v gcc-11 >/dev/null 2>&1 && command -v g++-11 >/dev/null 2>&1; then
  HIGHS_CMAKE_ARGS+=(-DCMAKE_C_COMPILER=gcc-11 -DCMAKE_CXX_COMPILER=g++-11)
fi
if [[ -n "${NVCC_PATH}" ]]; then
  echo "[configure] nvcc found: ${NVCC_PATH}"
  HIGHS_CMAKE_ARGS+=(-DCUPDLP_GPU=ON -DCMAKE_CUDA_ARCHITECTURES=native)
  if [[ -z "${CUDA_HOME:-}" ]]; then
    CUDA_HOME="$(dirname "$(dirname "${NVCC_PATH}")")"
    export CUDA_HOME
  fi
  # Extract CUDA version (e.g. 12.1) for compat patch decision
  CUDART_VERSION=$("${NVCC_PATH}" --version 2>/dev/null | grep -oP 'release \K[0-9]+\.[0-9]+' || echo "")
  echo "[configure] CUDA_HOME=${CUDA_HOME} CUDA=${CUDART_VERSION} -> CUPDLP_GPU=ON"
else
  echo "[configure] nvcc not found -> CUPDLP_GPU=OFF (CPU-only)"
  HIGHS_CMAKE_ARGS+=(-DCUPDLP_GPU=OFF)
fi

# ---- 2. Adaptive CMAKE_PREFIX_PATH (conda env / venv / custom) ----
PREFIX_PATHS=()
[[ -n "${CONDA_PREFIX:-}" ]] && PREFIX_PATHS+=("${CONDA_PREFIX}")
[[ -n "${VIRTUAL_ENV:-}" ]]  && PREFIX_PATHS+=("${VIRTUAL_ENV}")
[[ -n "${CMAKE_PREFIX_PATH_EXTRA:-}" ]] && PREFIX_PATHS+=("${CMAKE_PREFIX_PATH_EXTRA}")
JOIN_PREFIX=""
if [[ ${#PREFIX_PATHS[@]} -gt 0 ]]; then
  _ifs="${IFS}"; IFS=';'
  JOIN_PREFIX="${PREFIX_PATHS[*]}"
  IFS="${_ifs}"
  HIGHS_CMAKE_ARGS+=(-DCMAKE_PREFIX_PATH="${JOIN_PREFIX}")
  echo "[configure] CMAKE_PREFIX_PATH=${JOIN_PREFIX}"
fi

# ---- 3. CUDA <12.4 compatibility patch (applied to source tree, not committed) ----
# cusparseSpMV_preprocess only exists in CUDA 12.4+. For 12.1-12.3, comment out
# the two call sites in pdhg.cc. This patches the working tree but is never
# committed (zero upstream modifications in the repo).
if [[ -n "${NVCC_PATH}" ]] && [[ "${CUDART_VERSION:-}" == "12."[0-3]* ]]; then
  _pdhg="highs/pdlp/hipdlp/pdhg.cc"
  if grep -q 'cusparseSpMV_preprocess' "${_pdhg}" 2>/dev/null; then
    echo "[configure] CUDA <12.4 detected: patching ${_pdhg} (cusparseSpMV_preprocess)"
    sed -i '/CUSPARSE_CHECK(cusparseSpMV_preprocess(/,/));/s|^|//|' "${_pdhg}"
    # Mark so we can unpatch after build
    touch .patched-pdhg
  fi
fi

# ---- 4. Step 1: build & install HiGHS ----
HIGHS_CMAKE_ARGS+=(
  -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}"
  -DCMAKE_BUILD_TYPE=Release
  "${EXTRA_CMAKE_ARGS[@]}"
)
echo "[configure] Step 1: configure + build + install HiGHS"
echo "[configure]   cmake -S. -B${HIGHS_BUILD_DIR} ${HIGHS_CMAKE_ARGS[*]}"
cmake -S. -B"${HIGHS_BUILD_DIR}" "${HIGHS_CMAKE_ARGS[@]}"
cmake --build "${HIGHS_BUILD_DIR}" --parallel --target install

# Unpatch pdhg.cc after build (restore working tree to clean state)
if [[ -f .patched-pdhg ]]; then
  git checkout -- highs/pdlp/hipdlp/pdhg.cc 2>/dev/null || true
  rm -f .patched-pdhg
  echo "[configure] restored pdhg.cc to clean state"
fi

# ---- 4. Step 2: build highs-server against installed HiGHS ----
SERVER_CMAKE_ARGS=(
  -DCMAKE_PREFIX_PATH="${INSTALL_PREFIX}${JOIN_PREFIX:+;${JOIN_PREFIX}}"
  -DCMAKE_BUILD_TYPE=Release
)
echo "[configure] Step 2: configure + build highs-server"
echo "[configure]   cmake -Sserver -B${SERVER_BUILD_DIR} ${SERVER_CMAKE_ARGS[*]}"
cmake -Sserver -B"${SERVER_BUILD_DIR}" "${SERVER_CMAKE_ARGS[@]}"
cmake --build "${SERVER_BUILD_DIR}" --parallel --target highs_grpc_server

echo "[configure] done. Binary: ${SERVER_BUILD_DIR}/highs_grpc_server"
echo "[configure] next: ./${SERVER_BUILD_DIR}/highs_grpc_server --bind 127.0.0.1:50051"
