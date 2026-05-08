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

CSV_FIELDS = [
    "instance",
    "mode",
    "orig_rows",
    "orig_cols",
    "orig_nnz",
    "reduced_rows",
    "reduced_cols",
    "reduced_nnz",
    "pslp_time",
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
    "plain_model_status",
    "plain_iterations",
    "plain_solve_time",
    "plain_total_time",
    "plain_objective",
    "pslp_model_status",
    "pslp_iterations",
    "pslp_pslp_time",
    "pslp_solve_time",
    "pslp_total_time",
    "pslp_objective",
    "orig_rows",
    "orig_cols",
    "orig_nnz",
    "reduced_rows",
    "reduced_cols",
    "reduced_nnz",
    "nnz_reduction",
    "nnz_reduction_ratio",
    "objective_delta",
    "objective_relative_delta",
    "iteration_delta",
    "iteration_ratio",
    "total_time_delta",
    "total_time_ratio",
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
        if parts[1] not in {"plain", "pslp"}:
            continue
        return dict(zip(CSV_FIELDS, parts))
    raise ValueError("Could not find bench CSV line in command output")


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


def run_mode(instance: str, input_file: Path, mode: str) -> tuple[dict[str, str], str]:
    cmd = [str(BENCH_BIN), str(input_file), f"--mode={mode}"]
    completed = subprocess.run(
        cmd,
        cwd=BUILD_DIR,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"{BENCH_BIN.name} failed for {instance} mode={mode} with code {completed.returncode}\n"
            f"{completed.stdout}"
        )
    parsed = parse_bench_csv_line(completed.stdout)
    return parsed, completed.stdout


def build_comparison_row(instance: str, input_file: Path | None,
                         plain: dict[str, str] | None,
                         pslp: dict[str, str] | None,
                         availability: str) -> dict[str, str]:
    row = {field: "" for field in COMPARISON_FIELDS}
    row["instance"] = instance
    row["input_file"] = "" if input_file is None else str(input_file)
    row["availability"] = availability
    if not plain or not pslp:
        return row

    plain_total_time = to_float(plain["total_time"])
    pslp_total_time = to_float(pslp["total_time"])
    plain_objective = to_float(plain["objective"])
    pslp_objective = to_float(pslp["objective"])
    plain_iterations = to_int(plain["iterations"])
    pslp_iterations = to_int(pslp["iterations"])
    orig_nnz = to_int(plain["orig_nnz"])
    reduced_nnz = to_int(pslp["reduced_nnz"])

    row.update(
        {
            "plain_model_status": plain["model_status"],
            "plain_iterations": plain["iterations"],
            "plain_solve_time": plain["solve_time"],
            "plain_total_time": plain["total_time"],
            "plain_objective": plain["objective"],
            "pslp_model_status": pslp["model_status"],
            "pslp_iterations": pslp["iterations"],
            "pslp_pslp_time": pslp["pslp_time"],
            "pslp_solve_time": pslp["solve_time"],
            "pslp_total_time": pslp["total_time"],
            "pslp_objective": pslp["objective"],
            "orig_rows": plain["orig_rows"],
            "orig_cols": plain["orig_cols"],
            "orig_nnz": plain["orig_nnz"],
            "reduced_rows": pslp["reduced_rows"],
            "reduced_cols": pslp["reduced_cols"],
            "reduced_nnz": pslp["reduced_nnz"],
            "nnz_reduction": str(orig_nnz - reduced_nnz),
            "nnz_reduction_ratio": str(safe_ratio(orig_nnz - reduced_nnz, orig_nnz)),
            "iteration_delta": str(pslp_iterations - plain_iterations),
            "iteration_ratio": str(safe_ratio(pslp_iterations, plain_iterations)),
            "total_time_delta": str(pslp_total_time - plain_total_time),
            "total_time_ratio": str(safe_ratio(pslp_total_time, plain_total_time)),
        }
    )
    if has_reliable_objective(plain) and has_reliable_objective(pslp):
        row["objective_delta"] = str(pslp_objective - plain_objective)
        row["objective_relative_delta"] = str(
            safe_ratio(pslp_objective - plain_objective, max(1.0, abs(plain_objective)))
        )
    return row


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run experimental HiPDLP vs PSLP+HiPDLP benchmarks and write CSV files."
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
                    build_comparison_row(instance, None, None, None, "missing_input")
                )
                comparison_file.flush()
                continue

            print(f"Running {instance} plain ...", file=sys.stderr)
            plain_row, plain_stdout = run_mode(instance, input_file, "plain")
            print(f"Running {instance} pslp ...", file=sys.stderr)
            pslp_row, pslp_stdout = run_mode(instance, input_file, "pslp")

            (args.log_dir / f"{instance}_plain.log").write_text(plain_stdout)
            (args.log_dir / f"{instance}_pslp.log").write_text(pslp_stdout)

            raw_writer.writerow(plain_row)
            raw_writer.writerow(pslp_row)
            raw_file.flush()

            comparison_writer.writerow(
                build_comparison_row(instance, input_file, plain_row, pslp_row, "ok")
            )
            comparison_file.flush()

    print(args.raw_csv)
    print(args.comparison_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
