#!/usr/bin/env python3
import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT / "build"
BENCH_BIN = BUILD_DIR / "bin" / "hipdlp_pslp_bench"
LOG_DIR = BUILD_DIR / "hipdlp_pslp_logs"
RAW_CSV = BUILD_DIR / "hipdlp_pslp_bench_raw.csv"
COMPARISON_CSV = BUILD_DIR / "hipdlp_pslp_bench_comparison.csv"

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

MODES = ["none", "highs", "pslp"]

CSV_FIELDS = [
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

COMPARISON_FIELDS = [
    "instance",
    "input_file",
    "availability",
    "none_model_status",
    "none_iterations",
    "none_presolve_time",
    "none_solve_time",
    "none_total_time",
    "none_objective",
    "highs_model_status",
    "highs_iterations",
    "highs_presolve_time",
    "highs_solve_time",
    "highs_total_time",
    "highs_objective",
    "pslp_model_status",
    "pslp_iterations",
    "pslp_presolve_time",
    "pslp_solve_time",
    "pslp_total_time",
    "pslp_objective",
    "orig_rows",
    "orig_cols",
    "orig_nnz",
    "highs_reduced_rows",
    "highs_reduced_cols",
    "highs_reduced_nnz",
    "reduced_rows",
    "reduced_cols",
    "reduced_nnz",
    "highs_nnz_reduction",
    "highs_nnz_reduction_ratio",
    "nnz_reduction",
    "nnz_reduction_ratio",
    "highs_objective_delta",
    "highs_objective_relative_delta",
    "pslp_objective_delta",
    "pslp_objective_relative_delta",
    "highs_iteration_delta",
    "highs_iteration_ratio",
    "pslp_iteration_delta",
    "pslp_iteration_ratio",
    "highs_total_time_delta",
    "highs_total_time_ratio",
    "pslp_total_time_delta",
    "pslp_total_time_ratio",
]


def find_instance(name: str) -> Path | None:
    for base_dir in (ROOT / "check" / "instances", Path("/srv/mps_da")):
        for suffix in (".mps.gz", ".mps", ".lp", ".gz"):
            candidate = base_dir / f"{name}{suffix}"
            if candidate.exists():
                return candidate
    return None


def parse_bench_csv_line(output: str) -> dict[str, str]:
    for line in reversed(output.splitlines()):
        stripped = line.strip()
        parts = stripped.split(",")
        if len(parts) != len(CSV_FIELDS):
            continue
        if parts[1] not in {"none", "plain", "highs", "pslp"}:
            continue
        return dict(zip(CSV_FIELDS, parts))
    raise ValueError("Could not find bench CSV line in command output")


def make_incomplete_row(instance: str, mode: str, model_status: str,
                        total_time: float | None = None) -> dict[str, str]:
    row = {field: "" for field in CSV_FIELDS}
    row["instance"] = instance
    row["mode"] = mode
    row["model_status"] = model_status
    if total_time is not None:
        row["total_time"] = str(total_time)
    return row


def decode_subprocess_output(output: str | bytes | None) -> str:
    if output is None:
        return ""
    if isinstance(output, bytes):
        return output.decode("utf-8", errors="replace")
    return output


def to_float(value: str) -> float:
    return float(value) if value else math.nan


def to_int(value: str) -> int:
    return int(value) if value else 0


def safe_ratio(numerator: float, denominator: float) -> float:
    if denominator == 0 or math.isnan(denominator):
        return math.nan
    return numerator / denominator


def has_reliable_objective(run: dict[str, str]) -> bool:
    if run["model_status"] not in {"Optimal", "Model empty"}:
        return False
    try:
        return math.isfinite(to_float(run["objective"]))
    except ValueError:
        return False


def is_complete_run(run: dict[str, str] | None) -> bool:
    if not run:
        return False
    if run["model_status"].startswith("ExitCode"):
        return False
    return run["model_status"] != "Timeout"


def run_mode(instance: str, input_file: Path, mode: str,
             timeout_seconds: float | None) -> tuple[dict[str, str], str]:
    cmd = [str(BENCH_BIN), str(input_file), f"--mode={mode}"]
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
            make_incomplete_row(instance, mode, "Timeout", timeout_seconds),
            stdout + f"\nTIMEOUT after {timeout_seconds} seconds\n",
        )
    if completed.returncode != 0:
        return (
            make_incomplete_row(instance, mode, f"ExitCode{completed.returncode}"),
            completed.stdout,
        )
    parsed = parse_bench_csv_line(completed.stdout)
    return parsed, completed.stdout


def build_comparison_row(instance: str, input_file: Path | None,
                         none: dict[str, str] | None,
                         highs: dict[str, str] | None,
                         pslp: dict[str, str] | None,
                         availability: str) -> dict[str, str]:
    row = {field: "" for field in COMPARISON_FIELDS}
    row["instance"] = instance
    row["input_file"] = "" if input_file is None else str(input_file)
    row["availability"] = availability

    mode_runs = {"none": none, "highs": highs, "pslp": pslp}
    for mode, run in mode_runs.items():
        if not run:
            continue
        row[f"{mode}_model_status"] = run["model_status"]
        row[f"{mode}_iterations"] = run["iterations"]
        row[f"{mode}_presolve_time"] = run["presolve_time"]
        row[f"{mode}_solve_time"] = run["solve_time"]
        row[f"{mode}_total_time"] = run["total_time"]
        row[f"{mode}_objective"] = run["objective"]

    baseline = next((run for run in (none, highs, pslp) if run), None)
    if baseline:
        row["orig_rows"] = baseline["orig_rows"]
        row["orig_cols"] = baseline["orig_cols"]
        row["orig_nnz"] = baseline["orig_nnz"]
    if highs:
        row["highs_reduced_rows"] = highs["reduced_rows"]
        row["highs_reduced_cols"] = highs["reduced_cols"]
        row["highs_reduced_nnz"] = highs["reduced_nnz"]
    if pslp:
        row["reduced_rows"] = pslp["reduced_rows"]
        row["reduced_cols"] = pslp["reduced_cols"]
        row["reduced_nnz"] = pslp["reduced_nnz"]

    if baseline and highs and baseline["orig_nnz"] and highs["reduced_nnz"]:
        orig_nnz = to_int(baseline["orig_nnz"])
        highs_reduced_nnz = to_int(highs["reduced_nnz"])
        row["highs_nnz_reduction"] = str(orig_nnz - highs_reduced_nnz)
        row["highs_nnz_reduction_ratio"] = str(
            safe_ratio(orig_nnz - highs_reduced_nnz, orig_nnz)
        )
    if baseline and pslp and baseline["orig_nnz"] and pslp["reduced_nnz"]:
        orig_nnz = to_int(baseline["orig_nnz"])
        reduced_nnz = to_int(pslp["reduced_nnz"])
        row["nnz_reduction"] = str(orig_nnz - reduced_nnz)
        row["nnz_reduction_ratio"] = str(safe_ratio(orig_nnz - reduced_nnz, orig_nnz))

    if is_complete_run(none) and is_complete_run(highs):
        none_total_time = to_float(none["total_time"])
        highs_total_time = to_float(highs["total_time"])
        none_iterations = to_int(none["iterations"])
        highs_iterations = to_int(highs["iterations"])
        row["highs_iteration_delta"] = str(highs_iterations - none_iterations)
        row["highs_iteration_ratio"] = str(
            safe_ratio(highs_iterations, none_iterations)
        )
        row["highs_total_time_delta"] = str(highs_total_time - none_total_time)
        row["highs_total_time_ratio"] = str(
            safe_ratio(highs_total_time, none_total_time)
        )
    if is_complete_run(none) and is_complete_run(pslp):
        none_total_time = to_float(none["total_time"])
        pslp_total_time = to_float(pslp["total_time"])
        none_iterations = to_int(none["iterations"])
        pslp_iterations = to_int(pslp["iterations"])
        row["pslp_iteration_delta"] = str(pslp_iterations - none_iterations)
        row["pslp_iteration_ratio"] = str(
            safe_ratio(pslp_iterations, none_iterations)
        )
        row["pslp_total_time_delta"] = str(pslp_total_time - none_total_time)
        row["pslp_total_time_ratio"] = str(
            safe_ratio(pslp_total_time, none_total_time)
        )

    if is_complete_run(none) and is_complete_run(highs) and has_reliable_objective(none) and has_reliable_objective(highs):
        none_objective = to_float(none["objective"])
        highs_objective = to_float(highs["objective"])
        row["highs_objective_delta"] = str(highs_objective - none_objective)
        row["highs_objective_relative_delta"] = str(
            safe_ratio(highs_objective - none_objective, max(1.0, abs(none_objective)))
        )
    if is_complete_run(none) and is_complete_run(pslp) and has_reliable_objective(none) and has_reliable_objective(pslp):
        none_objective = to_float(none["objective"])
        pslp_objective = to_float(pslp["objective"])
        row["pslp_objective_delta"] = str(pslp_objective - none_objective)
        row["pslp_objective_relative_delta"] = str(
            safe_ratio(pslp_objective - none_objective, max(1.0, abs(none_objective)))
        )
    return row


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run experimental HiPDLP benchmarks for none/highs/pslp presolve modes and write CSV files."
    )
    parser.add_argument(
        "--instances",
        nargs="*",
        default=INSTANCES,
        help="Instance names to run. Defaults to the built-in experiment set.",
    )
    parser.add_argument(
        "--raw-csv",
        type=Path,
        default=RAW_CSV,
        help="Path for row-per-run CSV output.",
    )
    parser.add_argument(
        "--comparison-csv",
        type=Path,
        default=COMPARISON_CSV,
        help="Path for per-instance comparison CSV output.",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=LOG_DIR,
        help="Directory for per-instance raw stdout logs.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=300.0,
        help="Per-mode timeout in seconds. Timed-out runs are recorded and the batch continues.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not BENCH_BIN.exists():
        print(f"Missing binary: {BENCH_BIN}", file=sys.stderr)
        return 1

    args.raw_csv.parent.mkdir(parents=True, exist_ok=True)
    args.comparison_csv.parent.mkdir(parents=True, exist_ok=True)
    args.log_dir.mkdir(parents=True, exist_ok=True)

    with args.raw_csv.open("w", newline="") as raw_file, args.comparison_csv.open(
        "w", newline=""
    ) as comparison_file:
        raw_writer = csv.DictWriter(raw_file, fieldnames=CSV_FIELDS)
        comparison_writer = csv.DictWriter(comparison_file, fieldnames=COMPARISON_FIELDS)
        raw_writer.writeheader()
        comparison_writer.writeheader()

        for instance in args.instances:
            input_file = find_instance(instance)
            if input_file is None:
                print(f"Skipping {instance}: no input file found", file=sys.stderr)
                comparison_writer.writerow(
                    build_comparison_row(
                        instance, None, None, None, None, "missing_input"
                    )
                )
                comparison_file.flush()
                continue

            print(f"Running {instance} none ...", file=sys.stderr)
            none_row, none_stdout = run_mode(
                instance, input_file, "none", args.timeout_seconds
            )
            print(f"Running {instance} highs ...", file=sys.stderr)
            highs_row, highs_stdout = run_mode(
                instance, input_file, "highs", args.timeout_seconds
            )
            print(f"Running {instance} pslp ...", file=sys.stderr)
            pslp_row, pslp_stdout = run_mode(
                instance, input_file, "pslp", args.timeout_seconds
            )

            (args.log_dir / f"{instance}_none.log").write_text(none_stdout)
            (args.log_dir / f"{instance}_highs.log").write_text(highs_stdout)
            (args.log_dir / f"{instance}_pslp.log").write_text(pslp_stdout)

            raw_writer.writerow(none_row)
            raw_writer.writerow(highs_row)
            raw_writer.writerow(pslp_row)
            raw_file.flush()

            comparison_writer.writerow(
                build_comparison_row(
                    instance,
                    input_file,
                    none_row,
                    highs_row,
                    pslp_row,
                    "ok"
                    if all(is_complete_run(row) for row in (none_row, highs_row, pslp_row))
                    else "partial",
                )
            )
            comparison_file.flush()

    print(args.raw_csv)
    print(args.comparison_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
