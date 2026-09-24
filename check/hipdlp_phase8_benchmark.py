"""Phase8 representative HiPDLP gate; planning and summarizing never solve LPs."""

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import random
import selectors
import statistics
import subprocess
import sys
import time


NAMES = ("afiro", "e226", "pypsa_medium", "pypsa-eur-elec_50_168h")
BACKENDS = ("cpu", "cuda")
ROLES = ("legacy", "preoptimization", "candidate")
FIXED_OPTIONS = {"kkt_tolerance": 1e-4, "pdlp_optimality_tolerance": 1e-4,
                 "presolve": "off", "threads": 1, "precision": "fp64"}
DEFAULT_PROTOCOL = Path("/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/"
                        "20260908T000926Z-gpu-resident/performance_protocol.json")


def file_hash(filename):
    with Path(filename).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(destination, payload):
    encoded = json.dumps(payload, allow_nan=False, indent=2, sort_keys=True) + "\n"
    with Path(destination).open("x", encoding="utf-8") as stream:
        stream.write(encoded)


def read_json(source):
    payload = json.loads(Path(source).read_text(encoding="utf-8"))
    json.dumps(payload, allow_nan=False)
    return payload


def load_harness(protocol):
    spec = importlib.util.spec_from_file_location("lp_benchmark_reference", protocol["benchmark_harness"])
    module = importlib.util.module_from_spec(spec)
    previous = sys.dont_write_bytecode
    try:
        sys.dont_write_bytecode = True
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = previous
    return module


def load_instances(protocol):
    harness = load_harness(protocol)
    manifests = {name: harness.load_manifest(Path(name))
                 for name in {item["manifest"] for item in protocol["instances"]}}
    with Path(protocol["reference_csv"]).open(newline="", encoding="utf-8") as stream:
        references = [row for row in csv.DictReader(stream)
                      if row["scope"] == "energy-full" and row["solver"] == "highs-ipm"]
    instances = []
    for item in protocol["instances"]:
        matches = [row for row in manifests[item["manifest"]] if row["name"] == item["name"]]
        if len(matches) != 1:
            raise ValueError(f"Missing or duplicate manifest entry: {item['name']}")
        row = dict(matches[0], weight=item["weight"], reference_source=item["manifest"])
        if not row["reference_objective"]:
            candidates = [reference for reference in references if reference["instance"] == item["name"]]
            if (len(candidates) != 1 or candidates[0]["accepted"].lower() != "true"
                    or candidates[0]["status"] != "OPTIMAL"):
                raise ValueError(f"No unique accepted Phase7 IPM reference: {item['name']}")
            row.update(reference_objective=candidates[0]["primal_objective"],
                       reference_source=protocol["reference_csv"],
                       reference_source_directory=candidates[0].get("source_directory", ""))
        if not math.isfinite(float(row["reference_objective"])):
            raise ValueError(f"Nonfinite reference objective: {item['name']}")
        row["instance_path"] = str((Path(item["data_root"]) / row["path"]).resolve())
        instances.append(row)
    return instances


def compile_command(protocol, solver):
    source, build = Path(solver["source"]), Path(solver["build"])
    command = [solver["compiler"], "-std=c++11", "-O3", "-g", "-DNDEBUG", "-pthread"]
    if solver["role"] == "legacy":
        command.append("-DHIPDLP_PERFORMANCE_LEGACY")
    command += ["-I" + str(include) for include in
                (build, source / "highs", source / "app", source / "highs/pdlp/cupdlp",
                 source / "highs/pdlp/hipdlp")]
    return command + [str(Path(solver["driver"]).with_name("HiPdlpPerformance.cpp")),
                      solver["library"], "-Wl,-rpath," + solver["library_dir"],
                      "-o", solver["driver"]]


def runtime_command(artifact):
    prefix = ([artifact["loader"], "--library-path", artifact["library_dir"]]
              if artifact.get("loader") else [])
    return prefix + [artifact.get("driver", artifact.get("binary"))]


def make_plan(protocol, output):
    validate_protocol(protocol)
    instances, jobs, blockers = load_instances(protocol), [], []
    for solver in protocol["solvers"]:
        for artifact in (solver["driver"], solver["library"],
                         str(Path(solver["build"]) / "HConfig.h"),
                         str(Path(solver["build"]) / "CMakeCache.txt")):
            if not Path(artifact).is_file():
                blockers.append(f"{solver['id']}: missing {artifact}")
        for instance in instances:
            directory = Path(output).resolve() / "runs" / f"{solver['id']}--{instance['name']}"
            command = runtime_command(solver) + [
                "--instance", instance["instance_path"], "--backend", solver["backend"],
                "--step-size", protocol["step_size"], "--time-limit", str(protocol["time_limit_sec"]),
                "--iteration-limit", str(protocol["iteration_limit"]), "--output-dir", str(directory),
                "--warmups", str(protocol["warmups"]), "--repetitions", str(protocol["repetitions"]),
                "--seed", str(protocol["seed"])]
            jobs.append({"solver": solver, "instance": instance, "directory": str(directory),
                         "command": command, "timeout_sec": protocol["time_limit_sec"] + 30})
    for instance in instances:
        if not Path(instance["instance_path"]).is_file():
            blockers.append("missing instance: " + instance["instance_path"])
    random.Random(protocol["seed"]).shuffle(jobs)
    return {"schema": "hipdlp-phase8-plan/1", "protocol": protocol, "instances": instances,
            "jobs": jobs, "blockers": blockers,
            "compile_commands": {solver["id"]: compile_command(protocol, solver)
                                 for solver in protocol["solvers"]}}


def run_driver(command, directory, timeout_sec, environment):
    started = time.perf_counter()
    timed_out, next_sequence, pending = False, 0, b""
    with (Path(directory) / "stdout.log").open("xb") as output, \
            (Path(directory) / "stderr.log").open("xb") as errors:
        with subprocess.Popen(command, cwd=directory, env=environment,
                              stdout=subprocess.PIPE, stderr=errors) as child:
            selector = selectors.DefaultSelector()
            selector.register(child.stdout, selectors.EVENT_READ)
            deadline = time.monotonic() + timeout_sec
            try:
                while selector.get_map():
                    if not selector.select(max(0, deadline - time.monotonic())):
                        raise subprocess.TimeoutExpired(command, timeout_sec)
                    chunk = os.read(child.stdout.fileno(), 65536)
                    if not chunk:
                        selector.unregister(child.stdout)
                        break
                    output.write(chunk)
                    output.flush()
                    pending += chunk
                    while b"\n" in pending:
                        line, pending = pending.split(b"\n", 1)
                        try:
                            event = json.loads(line)
                        except (UnicodeError, ValueError):
                            continue
                        if (isinstance(event, dict) and event.get("event") == "begin"
                                and type(event.get("sequence")) is int
                                and event["sequence"] == next_sequence):
                            deadline = time.monotonic() + timeout_sec
                            next_sequence += 1
                child.wait(timeout=max(0, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                timed_out = True
            finally:
                if child.poll() is None:
                    child.kill()
                child.wait()
                output.write(child.stdout.read())
                selector.close()
    return {"command": command, "return_code": child.returncode, "timed_out": timed_out,
            "process_wall_sec": time.perf_counter() - started, "pid": child.pid,
            "per_repetition_timeout_sec": timeout_sec}


def child_environment(artifact):
    environment = {key: value for key, value in os.environ.items()
                   if key not in ("LD_PRELOAD", "LD_AUDIT", "LD_DEBUG", "LD_PROFILE", "LD_LIBRARY_PATH")}
    environment.update(LD_LIBRARY_PATH=artifact["library_dir"], OMP_NUM_THREADS="1",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    return environment


def production_inputs(solver):
    build, source = Path(solver["build"]), Path(solver["source"])
    cache = (build / "CMakeCache.txt").read_text().splitlines()
    for required in ("BUILD_TESTING:BOOL=OFF", "CMAKE_BUILD_TYPE:STRING=Release",
                     "CMAKE_HOME_DIRECTORY:INTERNAL=" + str(source.resolve())):
        if required not in cache:
            raise ValueError(f"{solver['id']}: production requires {required}")
    header = (build / "HConfig.h").read_text()
    defined = {line.split()[1] for line in header.splitlines() if line.startswith("#define ")}
    cuda_macro = "CUPDLP_GPU" if solver["role"] == "legacy" else "HIGHS_HIPDLP_CUDA_DETECTED"
    if (cuda_macro in defined) != (solver["backend"] == "cuda"):
        raise ValueError("Build backend disagrees with requested backend: " + solver["id"])
    commands = (build / "compile_commands.json").read_text()
    if "-DHIGHS_HIPDLP_TESTING" in commands:
        raise ValueError("Testing instrumentation is not a production baseline")
    library = Path(solver["library"])
    build_library = build / "lib" / "libhighs.so.1.15"
    if file_hash(library) != file_hash(build_library):
        raise ValueError("Library does not match its configured build")
    files = [build / name for name in ("CMakeCache.txt", "HConfig.h", "compile_commands.json")]
    files += [source / "CMakeLists.txt", source / "Version.txt", library, build_library]
    for folder in ("highs", "app", "extern", "cmake"):
        files += [file for file in (source / folder).rglob("*") if file.is_file()]
    return {str(file): file_hash(file) for file in sorted(set(files))}


def driver_process_failure(record):
    code = record.get("process_return_code")
    if type(code) is not int or code != 0:
        return f"driver_process_return_code={code}"
    if record.get("process_timed_out") is not False:
        return "driver_process_timeout_or_missing_state"
    return ""


def validate_correctness_evidence(protocol, evidence):
    expected = {(f"candidate-{backend}", name, strategy)
                for backend in BACKENDS
                for strategy, names in (("pid", NAMES), ("adaptive", NAMES[:2]))
                for name in names}
    if (evidence.get("all_accepted") is not True
            or evidence.get("planned") != len(expected)
            or evidence.get("completed") != len(expected)):
        raise ValueError("correctness evidence is not a complete accepted campaign")
    rows = evidence.get("rows", [])
    keys = []
    for row in rows:
        key = (row.get("solver"), row.get("instance"), row.get("requested_step_size"))
        keys.append(key)
        effective = "fixed" if key[2] == "adaptive" else "pid"
        if (row.get("accepted") is not True or row.get("status") != "OPTIMAL"
                or row.get("kkt_pass") is not True or row.get("kkt_input_valid") is not True
                or driver_process_failure(row) or row.get("effective_step_size") != effective
                or any(row.get(field) != value for field, value in FIXED_OPTIONS.items())):
            raise ValueError("correctness evidence contains an invalid candidate result")
        for field, limit in (("relative_primal_residual", 1e-4),
                             ("relative_dual_residual", 1e-4),
                             ("relative_objective_gap", 1e-4),
                             ("reference_objective_relative_error", 1e-3)):
            value = row.get(field)
            if (isinstance(value, bool) or not isinstance(value, (int, float))
                    or not math.isfinite(value) or not 0 <= value <= limit):
                raise ValueError("correctness evidence has invalid " + field)
        for field in ("seed", "iteration_limit", "time_limit_sec"):
            if field in protocol and row.get(field) != protocol[field]:
                raise ValueError("correctness evidence configuration differs: " + field)
    if len(keys) != len(expected) or set(keys) != expected:
        raise ValueError("correctness evidence has missing or duplicate coverage")


def accept_result(protocol, instance, driver, verification):
    result = dict(driver, **{key: value for key, value in verification.items()
                             if key not in ("status", "return_status", "iterations", "solver_time_sec")})
    result.setdefault("return_code", -1)
    result.setdefault("timed_out", False)
    acceptance = {"statuses": ["OPTIMAL"], "max_reference_objective_relative_error": 1e-3}
    acceptance.update({"max_" + field: 1e-4 for field in
                       ("relative_primal_residual", "relative_dual_residual", "relative_objective_gap")})
    result = load_harness(protocol).evaluate_result(instance, result, {"acceptance": acceptance})
    duration = result.get("solver_time_sec")
    reasons = []
    if (isinstance(duration, bool) or not isinstance(duration, (int, float))
            or not math.isfinite(duration) or duration <= 0):
        reasons.append("invalid_solver_time")
    if result.get("primal_objective") is None:
        reasons.append("missing_primal_objective")
    process_reason = driver_process_failure(driver)
    if process_reason:
        reasons.append(process_reason)
    if reasons:
        result["accepted"] = False
        result["rejection_reason"] = ";".join(filter(None, [result["rejection_reason"], *reasons]))
    return result


def verify_job(protocol, job, process):
    directory, solver, instance = Path(job["directory"]), job["solver"], job["instance"]
    rows = []
    requested = protocol["step_size"]
    legacy = solver["role"] == "legacy"
    effective = ("pid" if legacy else "fixed") if requested == "adaptive" else requested
    for kind, count in (("warmup", protocol["warmups"]), ("measurement", protocol["repetitions"])):
        for index in range(count):
            stem = f"{kind}-{index:03d}"
            identity = {"solver": solver["id"], "instance": instance["name"], "kind": kind, "index": index}
            try:
                row = read_json(directory / f"{stem}.json")
                sequence = index + (protocol["warmups"] if kind == "measurement" else 0)
                expected = dict(FIXED_OPTIONS, instance=instance["instance_path"], backend=solver["backend"],
                                implementation="legacy" if legacy else "rewrite", pid=process["pid"],
                                kind=kind, index=index, sequence=sequence, cold=sequence == 0,
                                requested_step_size=requested, requested_step_size_enum={"fixed": 0, "adaptive": 1, "pid": 3}[requested],
                                effective_step_size=effective, scaling="legacy-native-default" if legacy else "auto",
                                solution_file=f"{stem}.sol", time_limit_sec=protocol["time_limit_sec"])
                expected.update({field: protocol[field] for field in ("warmups", "repetitions", "seed", "iteration_limit")})
                mismatches = [key for key, value in expected.items() if row.get(key) != value]
                if mismatches:
                    raise ValueError("driver metadata mismatch: " + ",".join(mismatches))
                verified_path = directory / f"verified-{stem}.json"
                if verified_path.exists():
                    raise FileExistsError(verified_path)
                command = runtime_command(protocol["verifier"]) + [
                    "--instance", instance["instance_path"], "--solution", str(directory / f"{stem}.sol"),
                    "--solver", solver["id"], "--status", row["status"],
                    "--solver-time", str(row["solver_time_sec"]), "--iterations", str(row["iterations"]),
                    "--tolerance", "0.0001", "--output", str(verified_path)]
                started = time.perf_counter()
                with (directory / f"verifier-{stem}.stdout.log").open("xb") as output, \
                        (directory / f"verifier-{stem}.stderr.log").open("xb") as errors:
                    verified = subprocess.run(command, cwd=directory, env=child_environment(protocol["verifier"]),
                                              stdout=output, stderr=errors, timeout=protocol["time_limit_sec"] + 30)
                verification = read_json(verified_path) if verified_path.is_file() else {}
                verification.update(return_code=verified.returncode, timed_out=False)
                row.update(process_return_code=process["return_code"],
                           process_timed_out=process["timed_out"])
                row = accept_result(protocol, instance, row, verification)
                row.update(verifier_time_sec=time.perf_counter() - started, verifier_command=command)
            except (OSError, ValueError, KeyError, TypeError, subprocess.TimeoutExpired) as error:
                row = {"accepted": False, "rejection_reason": f"{type(error).__name__}: {error}",
                       "solver_time_sec": None}
            row.update(identity, process_wall_sec=process["process_wall_sec"],
                       process_return_code=process["return_code"], process_timed_out=process["timed_out"])
            write_json(directory / f"accepted-{stem}.json", row)
            rows.append(row)
    return rows


def verify_provenance(protocol, protocol_file):
    files = [Path(protocol_file), Path(__file__), Path(protocol["driver_source"]),
             Path(protocol["benchmark_harness"]), Path(protocol["reference_csv"])]
    files += [Path(item["manifest"]) for item in protocol["instances"]]
    files += [Path(protocol["verifier"]["binary"])]
    files += [file for file in Path(protocol["verifier"]["library_dir"]).iterdir() if file.is_file()]
    hashes = {str(file): file_hash(file) for file in sorted(set(files))}
    for solver in protocol["solvers"]:
        proof_file = Path(solver["driver"]).with_name("build-provenance.json")
        proof = read_json(proof_file)
        current = production_inputs(solver)
        if proof["compile_command"] != compile_command(protocol, solver) or proof["production_inputs"] != current:
            raise ValueError("Source/build/config drift for " + solver["id"])
        for filename, expected in proof["file_sha256"].items():
            if file_hash(filename) != expected:
                raise ValueError("Artifact drift: " + filename)
        if proof["file_sha256"][solver["driver"]] != file_hash(solver["driver"]):
            raise ValueError("Driver fingerprint missing or mismatched")
        if file_hash(Path(solver["driver"]).with_name("HiPdlpPerformance.cpp")) != file_hash(protocol["driver_source"]):
            raise ValueError("Driver source changed after its compilation freeze")
        hashes.update(current)
        hashes.update(proof["file_sha256"])
        hashes[str(proof_file)] = file_hash(proof_file)
    return hashes


def validate_protocol(protocol):
    if (sorted(item["name"] for item in protocol["instances"]) != sorted(NAMES)
            or any(item["weight"] != 1.0 for item in protocol["instances"])):
        raise ValueError("The four predeclared equally weighted instances are required")
    expected = sorted((f"{role}-{backend}", role, backend)
                      for role in ROLES for backend in BACKENDS)
    actual = sorted((item["id"], item["role"], item["backend"])
                    for item in protocol["solvers"])
    if actual != expected:
        raise ValueError("Both backends of legacy, preoptimization and candidate are required")
    if protocol["warmups"] < 1 or protocol["repetitions"] < 5:
        raise ValueError("At least one warmup and five measurements are required")
    if protocol["bootstrap_samples"] != 10000 or protocol["gamma"] != 1.05:
        raise ValueError("The predeclared gate uses 10000 bootstraps and gamma=1.05")
    if "step_size" in protocol and (protocol["step_size"] not in ("fixed", "pid", "adaptive")
            or not math.isfinite(protocol["time_limit_sec"]) or protocol["time_limit_sec"] <= 0
            or protocol.get("timeout_grace_sec", 30) != 30 or protocol["iteration_limit"] < 0):
        raise ValueError("Invalid solver settings or timeout grace")


def ratio_statistics(ratios, seed):
    if not ratios:
        return {"gm": None, "ci95": None}
    logs = [math.log(ratio) for ratio in ratios]
    generator = random.Random(seed)
    samples = sorted(math.exp(statistics.mean(generator.choices(logs, k=len(logs))))
                     for unused in range(10000))
    bounds = []
    for quantile in (0.025, 0.975):
        location = quantile * (len(samples) - 1)
        lower = math.floor(location)
        bounds.append(samples[lower] + (location - lower)
                      * (samples[lower + 1] - samples[lower]))
    return {"gm": math.exp(statistics.mean(logs)), "ci95": bounds}


def summarize(protocol, rows):
    validate_protocol(protocol)
    indexed = {}
    for row in rows:
        key = tuple(row.get(field) for field in ("solver", "instance", "kind", "index"))
        indexed.setdefault(key, []).append(row)
    unavailable, medians, expected_keys = [], {}, set()
    for solver in protocol["solvers"]:
        for instance in protocol["instances"]:
            measurements = []
            complete = True
            for kind, count in (("warmup", protocol["warmups"]),
                                ("measurement", protocol["repetitions"])):
                for index in range(count):
                    key = (solver["id"], instance["name"], kind, index)
                    expected_keys.add(key)
                    matches = indexed.get(key, [])
                    reason = ""
                    if len(matches) != 1:
                        reason = "missing" if not matches else "duplicate"
                    else:
                        row = matches[0]
                        duration = row.get("solver_time_sec")
                        if row.get("accepted") is not True:
                            reason = row.get("rejection_reason") or "rejected"
                        elif driver_process_failure(row):
                            reason = driver_process_failure(row)
                        elif (isinstance(duration, bool)
                              or not isinstance(duration, (int, float))
                              or not math.isfinite(duration) or duration <= 0):
                            reason = "invalid_solver_time"
                    if reason:
                        complete = False
                        unavailable.append(dict(zip(
                            ("solver", "instance", "kind", "index", "reason"),
                            (*key, reason))))
                    elif kind == "measurement":
                        measurements.append(matches[0]["solver_time_sec"])
            if complete:
                medians[(solver["id"], instance["name"])] = statistics.median(measurements)
    for key in indexed.keys() - expected_keys:
        unavailable.append({"unexpected_record": list(key), "reason": "unplanned"})
    backends, regressions, all_ratios = {}, [], []
    for backend in BACKENDS:
        pairs, before_ratios = [], []
        for name in NAMES:
            legacy = medians.get((f"legacy-{backend}", name))
            candidate = medians.get((f"candidate-{backend}", name))
            before = medians.get((f"preoptimization-{backend}", name))
            if candidate is not None and before is not None:
                before_ratios.append(candidate / before)
            if legacy is None or candidate is None:
                continue
            ratio = candidate / legacy
            pairs.append({"instance": name, "legacy_median_sec": legacy,
                          "candidate_median_sec": candidate, "ratio": ratio,
                          "legacy_cold_sec": indexed[(f"legacy-{backend}", name, "warmup", 0)][0]["solver_time_sec"],
                          "candidate_cold_sec": indexed[(f"candidate-{backend}", name, "warmup", 0)][0]["solver_time_sec"]})
            if ratio > 1.2:
                regressions.append({"backend": backend, "instance": name, "ratio": ratio})
        ratios = [pair["ratio"] for pair in pairs]
        all_ratios.extend(ratios)
        result = ratio_statistics(ratios, protocol["seed"])
        result.update(planned_pairs=4, jointly_accepted_pairs=len(pairs), pairs=pairs,
                      candidate_over_preoptimization=ratio_statistics(before_ratios, protocol["seed"]))
        result["verdict"] = "PASS" if (len(pairs) == 4 and result["gm"] <= protocol["gamma"]
                            and not any(item["backend"] == backend for item in regressions)) else "NOT_PASS"
        backends[backend] = result
    return {
        "scope": "provisional_representative_gate",
        "release_coverage": "BLOCKED_NOT_EVALUATED",
        "verdict": "PASS" if not unavailable and all(
            result["verdict"] == "PASS" for result in backends.values()) else "NOT_PASS",
        "incomplete": bool(unavailable), "unavailable": unavailable,
        "planned_pairs": 8, "jointly_accepted_pairs": len(all_ratios),
        "diagnostic_gm": (math.exp(statistics.mean(map(math.log, all_ratios)))
                          if all_ratios else None),
        "backends": backends, "regressions_over_20pct": regressions,
        "bootstrap_samples": 10000, "seed": protocol["seed"], "gamma": protocol["gamma"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("plan", "run", "summarize"))
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--correctness-evidence", type=Path)
    arguments = parser.parse_args()
    if arguments.action == "run" and (not arguments.execute or not arguments.correctness_evidence
                                       or not arguments.correctness_evidence.is_file()):
        parser.error("run requires --execute and --correctness-evidence FILE after main's serialized correctness approval")
    output = arguments.output_dir.resolve()
    if arguments.action == "summarize":
        plan = read_json(output / "plan.json")
        rows_file = output / "records.jsonl"
        rows = [json.loads(line) for line in rows_file.read_text().splitlines()] if rows_file.exists() else []
        summary = summarize(plan["protocol"], rows)
        completion = read_json(output / "completion.json") if (output / "completion.json").exists() else {}
        if completion.get("provenance_ok") is not True:
            summary.update(verdict="NOT_PASS", incomplete=True, provenance_ok=False)
        print(json.dumps(summary, allow_nan=False, indent=2))
        return 0 if summary["verdict"] == "PASS" else 2
    protocol = read_json(arguments.protocol)
    if arguments.action == "run":
        validate_correctness_evidence(protocol, read_json(arguments.correctness_evidence))
    plan = make_plan(protocol, output)
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "plan.json", plan)
    if arguments.action == "plan":
        print(json.dumps({"plan": str(output / "plan.json"), "planned_jobs": len(plan["jobs"]),
                          "blockers": plan["blockers"], "execution": "NOT_RUN"}, indent=2))
        return 0
    rows, hashes, blockers = [], {}, list(plan["blockers"])
    if not blockers:
        try:
            hashes = verify_provenance(protocol, arguments.protocol)
            hashes[str(arguments.correctness_evidence.resolve())] = file_hash(arguments.correctness_evidence)
            for instance in plan["instances"]:
                if file_hash(instance["instance_path"]) != instance["sha256"]:
                    raise ValueError("Model SHA256 mismatch: " + instance["name"])
                hashes[instance["instance_path"]] = instance["sha256"]
        except (OSError, ValueError, KeyError) as error:
            blockers.append(str(error))
    write_json(output / "input-sha256.json", hashes)
    if not blockers:
        with (output / "records.jsonl").open("x", encoding="utf-8") as stream:
            for job in plan["jobs"]:
                directory = Path(job["directory"])
                directory.mkdir(parents=True, exist_ok=False)
                environment = child_environment(job["solver"])
                process = run_driver(job["command"], directory, job["timeout_sec"], environment)
                process["environment"] = {key: environment.get(key) for key in
                    ("LD_LIBRARY_PATH", "CUDA_VISIBLE_DEVICES", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
                write_json(directory / "process.json", process)
                for row in verify_job(protocol, job, process):
                    rows.append(row)
                    stream.write(json.dumps(row, allow_nan=False) + "\n")
                    stream.flush()
        blockers += ["post-run artifact drift: " + filename for filename, digest in hashes.items()
                     if not Path(filename).is_file() or file_hash(filename) != digest]
    write_json(output / "completion.json", {"provenance_ok": not blockers, "blockers": blockers})
    summary = summarize(protocol, rows)
    if blockers:
        summary.update(verdict="NOT_PASS", incomplete=True, provenance_blockers=blockers)
    write_json(output / "summary.json", summary)
    print(json.dumps(summary, allow_nan=False, indent=2))
    return 0 if summary["verdict"] == "PASS" else 2


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError) as error:
        sys.exit(str(error))