import os
import subprocess
import time
import re

# List of lightweight benchmark instances from /home/yunus/Masaüstü/Projects/HiGHS/mps
instances = [
    "markshare_4_0.mps.gz",
    "pk1.mps.gz",
    "neos5.mps.gz",
    "mas76.mps.gz",
    "assign1-5-8.mps.gz"
]

mps_dir = "/home/yunus/Masaüstü/Projects/HiGHS/mps"
build_dir = "/home/yunus/Masaüstü/Projects/HiGHS/build"
results_dir = "/home/yunus/Masaüstü/Projects/HiGHS/benchmark_results"

solvers = {
    "original": os.path.join(build_dir, "bin/highs_original"),
    "optimized": os.path.join(build_dir, "bin/highs")
}

def parse_highs_output(stdout):
    # Regexes to extract information from solver output
    status_match = re.search(r"Status\s+(\w+)", stdout)
    primal_match = re.search(r"Primal bound\s+([-\d\.e+]+)", stdout)
    dual_match = re.search(r"Dual bound\s+([-\d\.e+]+)", stdout)
    nodes_match = re.search(r"Nodes\s+(\d+)", stdout)
    iters_match = re.search(r"LP iterations\s+(\d+)", stdout)
    time_match = re.search(r"Timing\s+([\d\.]+)", stdout)

    status = status_match.group(1) if status_match else "Unknown"
    primal = primal_match.group(1) if primal_match else "N/A"
    dual = dual_match.group(1) if dual_match else "N/A"
    nodes = nodes_match.group(1) if nodes_match else "0"
    iters = iters_match.group(1) if iters_match else "0"
    timing = time_match.group(1) if time_match else "0.00"

    return {
        "status": status,
        "primal": primal,
        "dual": dual,
        "nodes": nodes,
        "iters": iters,
        "timing": timing
    }

results = {inst: {} for inst in instances}

print(f"Starting comparison benchmarking on {len(instances)} instances...")

for inst in instances:
    inst_path = os.path.join(mps_dir, inst)
    inst_name = inst.split(".")[0]
    print(f"\n==================== Solving {inst} ====================")
    
    for solver_name, solver_bin in solvers.items():
        print(f"Running {solver_name} solver...")
        
        # Paths for saving log and solution
        log_file_path = os.path.join(results_dir, f"{inst_name}_{solver_name}.log")
        sol_file_path = os.path.join(results_dir, f"{inst_name}_{solver_name}.sol")
        
        # Prepare command args
        # We output solution to sol_file_path, and set time limit to 60 seconds
        cmd = [
            solver_bin,
            inst_path,
            "--time_limit", "60",
            "--solution_file", sol_file_path
        ]
        
        start_time = time.time()
        try:
            process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=90)
            stdout = process.stdout
            stderr = process.stderr
        except subprocess.TimeoutExpired as e:
            stdout = e.stdout if e.stdout else "TIMEOUT"
            stderr = e.stderr if e.stderr else ""
            print(f"Solver {solver_name} timed out!")
            
        # Save output log
        with open(log_file_path, "w") as f:
            f.write(stdout)
            if stderr:
                f.write("\n\n=== STDERR ===\n")
                f.write(stderr)
                
        metrics = parse_highs_output(stdout)
        results[inst][solver_name] = metrics
        print(f"  Result: {metrics['status']}, Obj: {metrics['primal']}, Nodes: {metrics['nodes']}, Time: {metrics['timing']}s")

# Generate beautiful markdown report
report_path = os.path.join(results_dir, "comparison_report.md")
with open(report_path, "w") as f:
    f.write("# HiGHS Original vs Cache-Optimized Solver Benchmark Comparison\n\n")
    f.write("This report presents a head-to-head comparison between the original HiGHS MIP solver and the L1/L2 cache-optimized version (prefetch + loop unrolling) on lightweight MIPLIB instances (solved with a 60-second time limit).\n\n")
    
    f.write("## Comparison Table\n\n")
    f.write("| Instance | Solver Version | Status | Primal Bound (Obj) | Dual Bound | Nodes | LP Iters | Time (s) |\n")
    f.write("| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |\n")
    
    for inst in instances:
        inst_name = inst.split(".")[0]
        orig = results[inst]["original"]
        opt = results[inst]["optimized"]
        
        orig_time = float(orig['timing'])
        opt_time = float(opt['timing'])
        speedup = orig_time / opt_time if opt_time > 0 else 0
        
        f.write(f"| **{inst_name}** | Original | {orig['status']} | {orig['primal']} | {orig['dual']} | {orig['nodes']} | {orig['iters']} | {orig['timing']} |\n")
        f.write(f"| | Optimized | {opt['status']} | {opt['primal']} | {opt['dual']} | {opt['nodes']} | {opt['iters']} | {opt['timing']} |\n")
        f.write(f"| | **Speedup** | | | | | | **{speedup:.2f}x** |\n")
        f.write("| --- | --- | --- | --- | --- | --- | --- | --- |\n")
        
    f.write("\n## Saved Artifacts\n")
    f.write("The `.sol` (solution) and `.log` (run log) files have been saved to `/home/yunus/Masaüstü/Projects/HiGHS/benchmark_results/` for all runs.\n")

print(f"\nBenchmarking complete! Report generated at {report_path}")
