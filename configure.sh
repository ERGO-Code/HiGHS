#!/usr/bin/env bash
# 自适应配置：有 CUDA 就开 GPU，没有就 CPU-only，二者皆可编译运行
#
# 用法:
#   ./configure.sh                       # 自适应：探测 nvcc 决定 GPU/CPU
#   ./configure.sh --no-cuda             # 强制 CPU-only（即便有 nvcc）
#   ./configure.sh --cuda                # 强制 GPU（无 nvcc 则报错退出）
#   ./configure.sh --build mybuild ...   # 指定 build 目录 + 透传额外 cmake 参数
#
# 环境变量:
#   CUDA_HOME          CUDA toolkit 根目录（默认从 nvcc 路径推导）
#   CMAKE_PREFIX_PATH_EXTRA  额外的 cmake 前缀路径（如 conda env / venv）
set -euo pipefail
cd "$(dirname "$0")"

# ---- 解析参数 ----
FORCE_GPU=""        # "" 自适应 / "on" 强制开 / "off" 强制关
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

# ---- 1. GPU 自适应探测 ----
NVCC_PATH=""
if [[ "$FORCE_GPU" != "off" ]]; then
  if command -v nvcc >/dev/null 2>&1; then
    NVCC_PATH="$(command -v nvcc)"
  elif [[ -n "${CUDA_HOME:-}" && -x "${CUDA_HOME}/bin/nvcc" ]]; then
    NVCC_PATH="${CUDA_HOME}/bin/nvcc"
  fi
fi

if [[ "$FORCE_GPU" == "on" && -z "$NVCC_PATH" ]]; then
  echo "[configure] 错误：--cuda 指定但未找到 nvcc" >&2
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
  echo "[configure] --no-cuda 指定 → CUPDLP_GPU=OFF（强制 CPU-only）"
  EXTRA_CMAKE_ARGS+=(-DCUPDLP_GPU=OFF)
else
  echo "[configure] nvcc 未找到 → CUPDLP_GPU=OFF（CPU-only，PDLP 仍可跑在 CPU）"
  EXTRA_CMAKE_ARGS+=(-DCUPDLP_GPU=OFF)
fi

# ---- 2. 强制开启 gRPC server 子项目 ----
EXTRA_CMAKE_ARGS+=(-DHIGHS_BUILD_GRPC_SERVER=ON)
EXTRA_CMAKE_ARGS+=(-DCMAKE_BUILD_TYPE=Release)

# ---- 3. 自适应 CMAKE_PREFIX_PATH（conda env / venv / 自定义）----
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
echo "[configure] 完成。下一步: cmake --build ${BUILD_DIR} --parallel --target highs_grpc_server"
