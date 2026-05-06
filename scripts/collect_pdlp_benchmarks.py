#!/usr/bin/env python3
import csv
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT / "build"
HIGHS_BIN = BUILD_DIR / "bin" / "highs"
OPTIONS_FILE = BUILD_DIR / "option.txt"
OUTPUT_CSV = BUILD_DIR / "pdlp_benchmarks.csv"

INSTANCES = [
    "shell",
    "stair",
    "fit2p",
    "fit2d",
    "25fv47",
    "adlittle",
    "neso-2005-filtered",
    "neso-2245-filtered",
    "L1_sixm1000obs"
]

SUMMARY_PATTERNS = {
    "termination_status": re.compile(r"Termination Status:\s*(.+)"),
    "model_status": re.compile(r"Model status\s*:\s*(.+)"),
    "total_iterations": re.compile(r"Total Iterations:\s*([0-9]+)"),
    "total_time_sec": re.compile(r"Total Time:\s*([0-9.]+)"),
    "final_primal_objective": re.compile(r"Final Primal Objective:\s*([-+0-9.eE]+)"),
    "final_duality_gap": re.compile(r"Final Duality Gap:\s*([-+0-9.eE]+)"),
    "final_primal_feasibility": re.compile(r"Final Primal Feasibility:\s*([-+0-9.eE]+)"),
    "final_dual_feasibility": re.compile(r"Final Dual Feasibility:\s*([-+0-9.eE]+)"),
    "pdlp_iterations": re.compile(r"PDLP\s+iterations:\s*([0-9]+)"),
    "objective_value": re.compile(r"Objective value\s*:\s*([-+0-9.eE]+)"),
    "pd_objective_error": re.compile(r"P-D objective error\s*:\s*([-+0-9.eE]+)"),
    "highs_run_time_sec": re.compile(r"HiGHS run time\s*:\s*([0-9.]+)"),
}


def find_instance(name: str) -> Path:
    for suffix in (".mps.gz", ".mps", ".lp", ".gz"):
        candidate = Path("/srv/mps_da") / f"{name}{suffix}"
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"No input file found for {name}")


def parse_summary(text: str) -> dict:
    result = {key: "" for key in SUMMARY_PATTERNS}
    for line in text.splitlines():
        for key, pattern in SUMMARY_PATTERNS.items():
            match = pattern.search(line)
            if match:
                result[key] = match.group(1).strip()
    return result


def has_useful_summary(summary: dict) -> bool:
    # Consider the run usable if at least one canonical end-of-run metric exists.
    return bool(
        summary.get("pdlp_iterations")
        or summary.get("objective_value")
        or summary.get("highs_run_time_sec")
    )


def main() -> int:
    if not HIGHS_BIN.exists():
        print(f"Missing binary: {HIGHS_BIN}", file=sys.stderr)
        return 1
    if not OPTIONS_FILE.exists():
        print(f"Missing options file: {OPTIONS_FILE}", file=sys.stderr)
        return 1

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "instance",
        "input_file",
        "highs_exit_code",
        "termination_status",
        "model_status",
        "total_iterations",
        "total_time_sec",
        "final_primal_objective",
        "final_duality_gap",
        "final_primal_feasibility",
        "final_dual_feasibility",
        "pdlp_iterations",
        "objective_value",
        "pd_objective_error",
        "highs_run_time_sec",
    ]
    with OUTPUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        f.flush()

        for instance in INSTANCES:
            input_file = find_instance(instance)
            cmd = [
                str(HIGHS_BIN),
                str(input_file),
                "--options_file",
                str(OPTIONS_FILE),
            ]
            print(f"Running {instance} ...", file=sys.stderr)
            completed = subprocess.run(
                cmd,
                cwd=BUILD_DIR,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
            summary = parse_summary(completed.stdout)
            if completed.returncode != 0 and not summary.get("model_status"):
                summary["model_status"] = "Unknown"
            if completed.returncode != 0 and not has_useful_summary(summary):
                print(completed.stdout, file=sys.stderr)
                raise RuntimeError(
                    f"Highs failed on {instance} with code {completed.returncode}"
                )
            if completed.returncode != 0:
                print(
                    (
                        f"Warning: HiGHS exited with code {completed.returncode} on {instance}, "
                        "but summary metrics were found; recording row."
                    ),
                    file=sys.stderr,
                )
            row = {"instance": instance, "input_file": str(input_file)}
            row["highs_exit_code"] = str(completed.returncode)
            row.update(summary)
            writer.writerow(row)
            f.flush()

    print(OUTPUT_CSV)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
