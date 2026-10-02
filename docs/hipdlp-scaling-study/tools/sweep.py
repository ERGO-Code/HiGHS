#!/usr/bin/env python3
"""Run the (scratch) migration adapter over instances x scaling configs, append CSV rows.

usage: sweep.py INSTANCE_LIST CSV_OUT TIME_LIMIT WORKERS [CONFIG,...] [BACKEND]
Rows already present in CSV_OUT are skipped, so the sweep can be resumed.
"""
import csv, json, os, subprocess, sys, threading, time
from concurrent.futures import ThreadPoolExecutor

SCRATCH = "/tmp/claude-1005/-home-yzhou-Github-hipdlp2/1e818617-21b5-47ff-b840-29674cc27652/scratchpad"
ADAPTER = os.environ.get("SWEEP_ADAPTER", f"{SCRATCH}/diag/adapter_x")
POOL = "/srv/mps_da"

CONFIGS = {
    "A":  ["--scaling", "auto"],                              # CR + Ruiz10 + PC  (current default)
    "L":  ["--scaling", "legacy"],                            # Ruiz10 + PC       (legacy scaling)
    "P":  ["--scaling", "legacy", "--scaling-mask", "4"],     # PC only
    "R":  ["--scaling", "legacy", "--scaling-mask", "1"],     # Ruiz only
    "C":  ["--scaling", "curtis_reid"],                       # CR only
    "O":  ["--scaling", "off"],
    "AP": ["--scaling", "auto", "--scaling-mask", "4"],       # CR + PC
    "AR": ["--scaling", "auto", "--scaling-mask", "1"],       # CR + Ruiz
}
FIELDS = ["instance", "config", "status", "iterations", "solver_wall_sec", "kkt_pass",
          "rel_primal", "rel_dual", "rel_gap", "objective"]

def main():
    inst_file, out, tl, workers = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4])
    cfgs = sys.argv[5].split(",") if len(sys.argv) > 5 and sys.argv[5] else list(CONFIGS)
    backend = sys.argv[6] if len(sys.argv) > 6 else "cuda"
    insts = [l.strip() for l in open(inst_file) if l.strip()]
    done = set()
    if os.path.exists(out):
        for r in csv.DictReader(open(out)):
            done.add((r["instance"], r["config"]))
    lock = threading.Lock()
    newfile = not os.path.exists(out)
    f = open(out, "a", newline="")
    w = csv.DictWriter(f, fieldnames=FIELDS)
    if newfile:
        w.writeheader(); f.flush()

    jobs = [(i, c) for i in insts for c in cfgs if (i, c) not in done]
    print(f"{len(jobs)} jobs", flush=True)

    def path_of(i):
        p = f"{SCRATCH}/mps/{i}.mps"
        return p if os.path.exists(p) else f"{POOL}/{i}.mps.gz"

    def run(job):
        i, c = job
        cmd = [ADAPTER, "--instance", path_of(i), "--implementation", "rewrite", "--backend", backend,
               "--kkt-tolerance", "1e-4", "--time-limit", str(tl)] + CONFIGS[c]
        row = {"instance": i, "config": c}
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=tl + 120)
            j = json.loads(p.stdout)
            row.update(status=j["status"], iterations=j["iterations"], solver_wall_sec=j["solver_wall_sec"],
                       kkt_pass=j["kkt_pass"], rel_primal=j["relative_primal_residual"],
                       rel_dual=j["relative_dual_residual"], rel_gap=j["relative_objective_gap"],
                       objective=j["objective"])
        except Exception as e:  # noqa
            row.update(status="RUNERR:" + type(e).__name__)
        with lock:
            w.writerow(row); f.flush()

    t0 = time.time()
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(run, jobs))
    print(f"done in {time.time()-t0:.0f}s", flush=True)

main()
