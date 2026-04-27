#!/usr/bin/env bash
# End-to-end PyPSA SciGRID-DE benchmark for HiGHS.
#
# Builds the MPS via build_scigrid_mps.py if missing, runs the highs CLI with
# a time limit, and emits a summary line that can be diffed across HiGHS
# builds. Total wall time, simplex iterations, and the objective value go to
# stdout; the full HiGHS log goes to <mps>.log next to the MPS.
#
# Usage:
#   check/perf/scripts/run_endtoend.sh [HOURS] [HIGHS_BIN] [TIME_LIMIT_SEC]
#
# Defaults: 168h (~50s, fast iteration). Use 504h (~5-6 min) for the proper
# real-workload signal.

set -euo pipefail

HOURS="${1:-168}"
HIGHS_BIN="${2:-build/bin/highs}"
TIME_LIMIT="${3:-1200}"

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT"

MPS="/tmp/scigrid_${HOURS}h.mps"
LOG="/tmp/scigrid_${HOURS}h.log"

if [[ ! -f "$MPS" ]]; then
  echo "[run_endtoend] building $MPS via pypsa..." >&2
  python3 "check/perf/scripts/build_scigrid_mps.py" \
    --hours "$HOURS" --out "$MPS"
fi

if [[ ! -x "$HIGHS_BIN" ]]; then
  echo "[run_endtoend] highs binary not found at $HIGHS_BIN" >&2
  echo "[run_endtoend]   build with: cmake --build build --target highs-bin" >&2
  exit 1
fi

echo "[run_endtoend] solving $MPS with $HIGHS_BIN (time_limit=${TIME_LIMIT}s)" >&2

# Use bash's TIMEFORMAT so we get a clean elapsed wall-time figure that's easy
# to grep for in the output.
TIMEFORMAT='[run_endtoend] wall_time_sec=%R'
{ time "$HIGHS_BIN" --time_limit="$TIME_LIMIT" "$MPS" > "$LOG" 2>&1 ; } 2> /tmp/_runtime

# Parse a few headline numbers from the HiGHS log.
ITERS=$(grep -E "Simplex   iterations:" "$LOG" | awk '{print $NF}' || echo "?")
OBJ=$(grep -E "Objective value" "$LOG" | awk '{print $NF}' || echo "?")
HIGHS_TIME=$(grep -E "HiGHS run time" "$LOG" | awk '{print $NF}' || echo "?")
STATUS=$(grep -E "Model status" "$LOG" | awk -F': *' '{print $2}' || echo "?")
WALL=$(grep -E "wall_time_sec" /tmp/_runtime | awk -F= '{print $2}')

echo "instance=scigrid_${HOURS}h status=${STATUS} iters=${ITERS} obj=${OBJ} highs_time_s=${HIGHS_TIME} wall_time_s=${WALL}"
