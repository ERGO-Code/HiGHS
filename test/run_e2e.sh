#!/usr/bin/env bash
# E2E 测试脚本：修复 v1 sleep 2 不可靠、kill 无 trap 兜底
set -euo pipefail
cd "$(dirname "$0")/.."

# 自动探测二进制位置（CMake 可能放 build/ 或 build/bin/）
if [[ -z "${SERVER_BIN:-}" ]]; then
  for cand in ./build/bin/highs_grpc_server ./build/highs_grpc_server; do
    [[ -x "$cand" ]] && SERVER_BIN="$cand" && break
  done
fi
SERVER_BIN="${SERVER_BIN:?未找到 highs_grpc_server 二进制}"
BIND="${BIND:-127.0.0.1:50051}"

echo "[e2e] starting server: ${SERVER_BIN} --bind ${BIND}"
"${SERVER_BIN}" --bind "${BIND}" --max-concurrent 1 &
SERVER_PID=$!
trap 'kill ${SERVER_PID} 2>/dev/null || true; wait ${SERVER_PID} 2>/dev/null || true' EXIT INT TERM

# 轮询就绪（最多 30s），替代脆弱的 sleep 2
echo "[e2e] waiting for server ready..."
for i in $(seq 1 60); do
  if python3 -c "import grpc; grpc.insecure_channel('${BIND}').channel_ready_future().result(timeout=1)" 2>/dev/null; then
    echo "[e2e] server ready"
    break
  fi
  sleep 0.5
done

# 生成 Python 桩
cd test
echo "[e2e] generating python stubs..."
python3 -m grpc_tools.protoc -I../server/protos --python_out=. --grpc_python_out=. ../server/protos/solver.proto

echo "[e2e] running tests..."
python3 test_client.py
echo "[e2e] ALL TESTS PASSED"
