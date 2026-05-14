#!/usr/bin/env python3
import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT / "build"
HIPDLP_BENCH_BIN = BUILD_DIR / "bin" / "hipdlp_pslp_bench"
CUPDLPX_BENCH_BIN = BUILD_DIR / "bin" / "highs_cupdlpx_bench"
LOG_DIR = BUILD_DIR / "highs_cupdlpx_logs"
RAW_CSV = BUILD_DIR / "highs_cupdlpx_bench_raw.csv"
COMPARISON_CSV = BUILD_DIR / "highs_cupdlpx_bench_comparison.csv"

INSTANCES = [
    "shell",
    "stair",
    "fit2p",
    "fit2d",
    "25fv47",
    "adlittle",
    "neso-2005-filtered",
    "neso-2245-filtered",
    "L1_sixm1000obs",
]

HIPDLP_FIELDS = [
    "instance",
    "mode",
    "orig_rows",
    "orig_cols",
    "orig_nnz",
    "reduced_rows",
    "reduced_cols",
    "reduced_nnz",
    "presolve_time",
    "iterations",
    "solve_time",
    "total_time",
    "objective",
    "model_status",
]

CUPDLPX_FIELDS = [
    "instance",
    "orig_rows",
    "orig_cols",
    "orig_nnz",
    "reduced_rows",
    "reduced_cols",
    "reduced_nnz",
    "presolve_time",
    "rescaling_time",
    "iterations",
    "solve_time",
    "total_time",
    "objective",
    "model_status",
    "termination_reason",
]

RAW_FIELDS = [
    "instance",
    "input_file",
    "solver",
    "orig_rows",
    "orig_cols",
    "orig_nnz",
    "reduced_rows",
    "reduced_cols",
    "reduced_nnz",
    "presolve_time",
    "rescaling_time",
    "iterations",
    "solve_time",
    "total_time",
    "objective",
    "model_status",
    "termination_reason",
]

COMPARISON_FIELDS = [
    "instance",
    "input_file",
    "availability",
    "hipdlp_model_status",
    "hipdlp_iterations",
    "hipdlp_presolve_time",
    "hipdlp_solve_time",
    "hipdlp_total_time",
    "hipdlp_objective",
    "cupdlpx_model_status",
    "cupdlpx_termination_reason",
    "cupdlpx_iterations",
    "cupdlpx_presolve_time",
    "cupdlpx_rescaling_time",
    "cupdlpx_solve_time",
    "cupdlpx_total_time",
    "cupdlpx_objective",
    "orig_rows",
    "orig_cols",
    "orig_nnz",
    "reduced_rows",
    "reduced_cols",
    "reduced_nnz",
    "nnz_reduction",
    "nnz_reduction_ratio",
    "objective_delta_cupdlpx_vs_hipdlp",
    "objective_relative_delta_cupdlpx_vs_hipdlp",
    "iteration_delta_cupdlpx_vs_hipdlp",
    "iteration_ratio_cupdlpx_vs_hipdlp",
    "total_time_delta_cupdlpx_vs_hipdlp",
    "total_time_ratio_cupdlpx_vs_hipdlp",
]


def find_instance(name: str) -> Path | None:
    for base_dir in (ROOT / "check" / "instances", Path("/srv/mps_da")):
        for suffix in (".mps.gz", ".mps", ".lp", ".gz"):
            candidate = base_dir / f"{name}{suffix}"
            if candidate.exists():
                return candidate
    return None


def parse_csv_line(output: str, fields: list[str]) -> dict[str, str]:
    for line in reversed(output.splitlines()):
        stripped = line.strip()
        parts = stripped.split(",")
        if len(parts) != len(fields):
            continue
        return dict(zip(fields, parts))
    raise ValueError("Could not find benchmark CSV line in command output")


def decode_subprocess_output(output: str | bytes | None) -> str:
    if output is None:
        return ""
    if isinstance(output, bytes):
        return output.decode("utf-8", errors="replace")
    return output


def make_incomplete_row(
    instance: str, solver: str, model_status: str, total_time: float | None = None
) -> dict[str, str]:
    row = {field: "" for field in RAW_FIELDS}
    row["instance"] = instance
    row["solver"] = solver
    row["model_status"] = model_status
    if total_time is not None:
        row["total_time"] = str(total_time)
    return row


def normalize_hipdlp_run(run: dict[str, str]) -> dict[str, str]:
    row = {field: "" for field in RAW_FIELDS}
    row["instance"] = run["instance"]
    row["solver"] = "hipdlp_highs_presolve"
    row["orig_rows"] = run["orig_rows"]
    row["orig_cols"] = run["orig_cols"]
    row["orig_nnz"] = run["orig_nnz"]
    row["reduced_rows"] = run["reduced_rows"]
    row["reduced_cols"] = run["reduced_cols"]
    row["reduced_nnz"] = run["reduced_nnz"]
    row["presolve_time"] = run["presolve_time"]
    row["iterations"] = run["iterations"]
    row["solve_time"] = run["solve_time"]
    row["total_time"] = run["total_time"]
    row["objective"] = run["objective"]
    row["model_status"] = run["model_status"]
    return row


def normalize_cupdlpx_run(run: dict[str, str]) -> dict[str, str]:
    row = {field: "" for field in RAW_FIELDS}
    row["instance"] = run["instance"]
    row["solver"] = "cupdlpx_highs_presolve"
    for field in (
        "orig_rows",
        "orig_cols",
        "orig_nnz",
        "reduced_rows",
        "reduced_cols",
        "reduced_nnz",
        "presolve_time",
        "rescaling_time",
        "iterations",
        "solve_time",
        "total_time",
        "objective",
        "model_status",
        "termination_reason",
    ):
        row[field] = run[field]
    return row


def to_float(value: str) -> float:
    return float(value) if value else math.nan


def to_int(value: str) -> int:
    return int(value) if value else 0


def safe_ratio(numerator: float, denominator: float) -> float:
    if denominator == 0 or math.isnan(denominator):
        return math.nan
    return numerator / denominator


def has_reliable_objective(run: dict[str, str] | None) -> bool:
    if not run:
        return False
    if run["model_status"] not in {"Optimal", "Model empty"}:
        return False
    try:
        return math.isfinite(to_float(run["objective"]))
    except ValueError:
        return False


def run_benchmark(
    cmd: list[str], fields: list[str], instance: str, solver: str, timeout_seconds: float | None
) -> tuple[dict[str, str], str]:
    try:
        completed = subprocess.run(
            cmd,
            cwd=BUILD_DIR,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as exc:
        stdout = decode_subprocess_output(exc.stdout)
        return (
            make_incomplete_row(instance, solver, "Timeout", timeout_seconds),
            stdout + f"\nTIMEOUT after {timeout_seconds} seconds\n",
        )
    if completed.returncode != 0:
        return (
            make_incomplete_row(instance, solver, f"ExitCode{completed.returncode}"),
            completed.stdout,
        )
    parsed = parse_csv_line(completed.stdout, fields)
    return parsed, completed.stdout


def build_comparison_row(
    instance: str,
    input_file: Path | None,
    hipdlp: dict[str, str] | None,
    cupdlpx: dict[str, str] | None,
    availability: str,
) -> dict[str, str]:
    row = {field: "" for field in COMPARISON_FIELDS}
    row["instance"] = instance
    row["input_file"] = "" if input_file is None else str(input_file)
    row["availability"] = availability

    if hipdlp:
        row["hipdlp_model_status"] = hipdlp["model_status"]
        row["hipdlp_iterations"] = hipdlp["iterations"]
        row["hipdlp_presolve_time"] = hipdlp["presolve_time"]
        row["hipdlp_solve_time"] = hipdlp["solve_time"]
        row["hipdlp_total_time"] = hipdlp["total_time"]
        row["hipdlp_objective"] = hipdlp["objective"]

    if cupdlpx:
        row["cupdlpx_model_status"] = cupdlpx["model_status"]
        row["cupdlpx_termination_reason"] = cupdlpx["termination_reason"]
        row["cupdlpx_iterations"] = cupdlpx["iterations"]
        row["cupdlpx_presolve_time"] = cupdlpx["presolve_time"]
        row["cupdlpx_rescaling_time"] = cupdlpx["rescaling_time"]
        row["cupdlpx_solve_time"] = cupdlpx["solve_time"]
        row["cupdlpx_total_time"] = cupdlpx["total_time"]
        row["cupdlpx_objective"] = cupdlpx["objective"]

    baseline = cupdlpx or hipdlp
    if baseline:
        row["orig_rows"] = baseline["orig_rows"]
        row["orig_cols"] = baseline["orig_cols"]
        row["orig_nnz"] = baseline["orig_nnz"]
        row["reduced_rows"] = baseline["reduced_rows"]
        row["reduced_cols"] = baseline["reduced_cols"]
        row["reduced_nnz"] = baseline["reduced_nnz"]
        orig_nnz = to_int(baseline["orig_nnz"])
        reduced_nnz = to_int(baseline["reduced_nnz"])
        row["nnz_reduction"] = str(orig_nnz - reduced_nnz)
        row["nnz_reduction_ratio"] = str(safe_ratio(orig_nnz - reduced_nnz, orig_nnz))

    if hipdlp and cupdlpx:
        hipdlp_iterations = to_int(hipdlp["iterations"])
        cupdlpx_iterations = to_int(cupdlpx["iterations"])
        hipdlp_total_time = to_float(hipdlp["total_time"])
        cupdlpx_total_time = to_float(cupdlpx["total_time"])
        row["iteration_delta_cupdlpx_vs_hipdlp"] = str(
            cupdlpx_iterations - hipdlp_iterations
        )
        row["iteration_ratio_cupdlpx_vs_hipdlp"] = str(
            safe_ratio(cupdlpx_iterations, hipdlp_iterations)
        )
        row["total_time_delta_cupdlpx_vs_hipdlp"] = str(
            cupdlpx_total_time - hipdlp_total_time
        )
        row["total_time_ratio_cupdlpx_vs_hipdlp"] = str(
            safe_ratio(cupdlpx_total_time, hipdlp_total_time)
        )

    if has_reliable_objective(hipdlp) and has_reliable_objective(cupdlpx):
        hipdlp_objective = to_float(hipdlp["objective"])
        cupdlpx_objective = to_float(cupdlpx["objective"])
        delta = cupdlpx_objective - hipdlp_objective
        scale = max(1.0, abs(hipdlp_objective), abs(cupdlpx_objective))
        row["objective_delta_cupdlpx_vs_hipdlp"] = str(delta)
        row["objective_relative_delta_cupdlpx_vs_hipdlp"] = str(delta / scale)

    return row


def write_log(log_path: Path, output: str) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(output)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--timeout",
        type=float,
        default=None,
        help="Per-run timeout in seconds.",
    )
    parser.add_argument(
        "--instance",
        action="append",
        dest="instances",
        help="Run only the named instance. Repeat to run multiple instances.",
    )
    parser.add_argument(
        "--solver",
        choices=("both", "hipdlp", "cupdlpx"),
        default="both",
        help="Run both solvers or just one of them.",
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop immediately when one benchmark command fails.",
    )
    args = parser.parse_args()

    required_bins = []
    if args.solver in ("both", "hipdlp"):
        required_bins.append(HIPDLP_BENCH_BIN)
    if args.solver in ("both", "cupdlpx"):
        required_bins.append(CUPDLPX_BENCH_BIN)

    missing = [str(path) for path in required_bins if not path.exists()]
    if missing:
        print("Missing benchmark binaries:", file=sys.stderr)
        for path in missing:
            print(f"  {path}", file=sys.stderr)
        return 1

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    instances = args.instances if args.instances else INSTANCES

    with RAW_CSV.open("w", newline="") as raw_file, COMPARISON_CSV.open(
        "w", newline=""
    ) as comparison_file:
        raw_writer = csv.DictWriter(raw_file, fieldnames=RAW_FIELDS)
        comparison_writer = csv.DictWriter(comparison_file, fieldnames=COMPARISON_FIELDS)
        raw_writer.writeheader()
        comparison_writer.writeheader()

        for instance in instances:
            input_file = find_instance(instance)
            if input_file is None:
                comparison_writer.writerow(
                    build_comparison_row(instance, None, None, None, "missing_input")
                )
                continue

            hipdlp_cmd = [str(HIPDLP_BENCH_BIN), str(input_file), "--mode=highs"]
            cupdlpx_cmd = [str(CUPDLPX_BENCH_BIN), str(input_file)]

            hipdlp_run = None
            cupdlpx_run = None

            if args.solver in ("both", "hipdlp"):
                print(f"Running {instance} with HiPDLP ...", file=sys.stderr)
                hipdlp_run, hipdlp_output = run_benchmark(
                    hipdlp_cmd,
                    HIPDLP_FIELDS,
                    instance,
                    "hipdlp_highs_presolve",
                    args.timeout,
                )
                hipdlp_log_path = LOG_DIR / f"{instance}_hipdlp.log"
                write_log(hipdlp_log_path, hipdlp_output)
                hipdlp_row = normalize_hipdlp_run(hipdlp_run)
                hipdlp_row["input_file"] = str(input_file)
                raw_writer.writerow(hipdlp_row)
                if "mode" not in hipdlp_run:
                    print(
                        f"HiPDLP failed on {instance}; see {hipdlp_log_path}",
                        file=sys.stderr,
                    )
                    if args.fail_fast:
                        comparison_writer.writerow(
                            build_comparison_row(
                                instance, input_file, None, None, "hipdlp_failed"
                            )
                        )
                        raw_file.flush()
                        comparison_file.flush()
                        return 1

            if args.solver in ("both", "cupdlpx"):
                print(f"Running {instance} with cuPDLPx ...", file=sys.stderr)
                cupdlpx_run, cupdlpx_output = run_benchmark(
                    cupdlpx_cmd,
                    CUPDLPX_FIELDS,
                    instance,
                    "cupdlpx_highs_presolve",
                    args.timeout,
                )
                cupdlpx_log_path = LOG_DIR / f"{instance}_cupdlpx.log"
                write_log(cupdlpx_log_path, cupdlpx_output)
                cupdlpx_row = normalize_cupdlpx_run(cupdlpx_run)
                cupdlpx_row["input_file"] = str(input_file)
                raw_writer.writerow(cupdlpx_row)
                if "termination_reason" not in cupdlpx_run:
                    print(
                        f"cuPDLPx failed on {instance}; see {cupdlpx_log_path}",
                        file=sys.stderr,
                    )
                    if args.fail_fast:
                        comparison_writer.writerow(
                            build_comparison_row(
                                instance, input_file, hipdlp_run, None, "cupdlpx_failed"
                            )
                        )
                        raw_file.flush()
                        comparison_file.flush()
                        return 1

            hipdlp_complete = hipdlp_run is not None and "mode" in hipdlp_run
            cupdlpx_complete = (
                cupdlpx_run is not None and "termination_reason" in cupdlpx_run
            )
            if hipdlp_complete and cupdlpx_complete:
                availability = "ok"
            elif hipdlp_complete:
                availability = "cupdlpx_failed"
            elif cupdlpx_complete:
                availability = "hipdlp_failed"
            else:
                availability = "both_failed"

            comparison_writer.writerow(
                build_comparison_row(
                    instance,
                    input_file,
                    hipdlp_run if hipdlp_complete else None,
                    cupdlpx_run if cupdlpx_complete else None,
                    availability,
                )
            )

            raw_file.flush()
            comparison_file.flush()

    print(COMPARISON_CSV)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
