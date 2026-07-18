#!/usr/bin/env bash
# E2E test script: fixes v1 unreliable sleep 2 and missing trap cleanup
set -euo pipefail
cd "$(dirname "$0")/../.."

# Auto-detect binary location (CMake may put in build/ or build/bin/)
if [[ -z "${SERVER_BIN:-}" ]]; then
  for cand in ./build/bin/highs_grpc_server ./build/highs_grpc_server; do
    [[ -x "$cand" ]] && SERVER_BIN="$cand" && break
  done
fi
SERVER_BIN="${SERVER_BIN:?highs_grpc_server binary not found}"
BIND="${BIND:-127.0.0.1:50051}"

echo "[e2e] starting server: ${SERVER_BIN} --bind ${BIND}"
"${SERVER_BIN}" --bind "${BIND}" --max-concurrent 1 &
SERVER_PID=$!
trap 'kill ${SERVER_PID} 2>/dev/null || true; wait ${SERVER_PID} 2>/dev/null || true' EXIT INT TERM

# Poll for readiness (up to 30s), replaces fragile sleep 2
echo "[e2e] waiting for server ready..."
for i in $(seq 1 60); do
  if python3 -c "import grpc; grpc.insecure_channel('${BIND}').channel_ready_future().result(timeout=1)" 2>/dev/null; then
    echo "[e2e] server ready"
    break
  fi
  sleep 0.5
done

# Generate Python stubs
cd server/test
echo "[e2e] generating python stubs..."
python3 -m grpc_tools.protoc -I../protos --python_out=. --grpc_python_out=. ../protos/solver.proto

echo "[e2e] running tests..."
python3 test_client.py
echo "[e2e] ALL TESTS PASSED"
