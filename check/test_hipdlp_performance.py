import copy
import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("hipdlp_phase8_benchmark.py")
NAMES = ("afiro", "e226", "pypsa_medium", "pypsa-eur-elec_50_168h")
RUNROOT = Path("/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/"
               "20260908T000926Z-gpu-resident")
AFIRO = Path("/home/yanyu/benchmark-data/lp/netlib/feasible/afiro.mps")


def synthetic_protocol():
    return {
        "instances": [{"name": name, "weight": 1.0} for name in NAMES],
        "solvers": [
            {"id": f"{role}-{backend}", "role": role, "backend": backend}
            for role in ("legacy", "preoptimization", "candidate")
            for backend in ("cpu", "cuda")
        ],
        "warmups": 1,
        "repetitions": 5,
        "seed": 20260908,
        "bootstrap_samples": 10000,
        "gamma": 1.05,
    }


def synthetic_rows(protocol):
    return [
        {
            "solver": solver["id"],
            "instance": instance["name"],
            "kind": kind,
            "index": index,
            "accepted": True,
            "process_return_code": 0,
            "process_timed_out": False,
            "solver_time_sec": 1.0 if solver["role"] == "candidate" else 2.0,
        }
        for solver in protocol["solvers"]
        for instance in protocol["instances"]
        for kind, count in (("warmup", protocol["warmups"]),
                            ("measurement", protocol["repetitions"]))
        for index in range(count)
    ]


@unittest.skipUnless(os.environ.get("HIPDLP_KKT_VERIFIER"),
                     "set HIPDLP_KKT_VERIFIER to test the external verifier")
class OriginalKktReportTests(unittest.TestCase):
    def test_report_uses_verified_offset_free_objective_accumulation(self):
        verifier = Path(os.environ["HIPDLP_KKT_VERIFIER"]).resolve()
        self.assertTrue(verifier.is_file())
        for sign in (1, -1):
            with self.subTest(sign=sign), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                model = root / "cancellation.mps"
                primal = root / "primal.txt"
                dual = root / "dual.txt"
                output = root / "verification.json"
                model.write_text(
                    "NAME CANCELLATION\nROWS\n N OBJ\nCOLUMNS\n"
                    f" X0 OBJ {sign}\n X1 OBJ {sign}\n X2 OBJ {-sign * 10**16}\n"
                    f"RHS\n RHS1 OBJ {-sign * 10**16}\nBOUNDS\n"
                    " FX BND X0 1\n FX BND X1 1\n FX BND X2 1\nENDATA\n")
                primal.write_text("1\n1\n1\n")
                dual.write_text("")
                result = subprocess.run([
                    str(verifier), "--instance", str(model),
                    "--primal-values", str(primal), "--row-dual-values", str(dual),
                    "--solver", "fixture", "--status", "OPTIMAL", "--solver-time", "0",
                    "--iterations", "0", "--tolerance", "0.0001", "--output", str(output)
                ], text=True, capture_output=True, timeout=10, check=False)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                record = json.loads(output.read_text())
                self.assertEqual(record["kkt_contract"], "hipdlp-original-kkt/2")
                self.assertTrue(record["kkt_pass"])
                self.assertTrue(record["kkt_complementarity_ok"])
                self.assertEqual(record["relative_objective_gap"], 0.0)
                self.assertEqual(record["primal_objective"], sign * 2.0)


class HarnessTest(unittest.TestCase):
    def setUp(self):
        self.assertTrue(SCRIPT.is_file(), "Phase8 benchmark implementation is missing")
        spec = importlib.util.spec_from_file_location("phase8_benchmark", SCRIPT)
        self.harness = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.harness)
        self.protocol = synthetic_protocol()
        self.rows = synthetic_rows(self.protocol)


class GateTests(HarnessTest):
    def test_missing_cuda_never_passes(self):
        rows = [row for row in self.rows if row["solver"].endswith("-cpu")]
        summary = self.harness.summarize(self.protocol, rows)
        self.assertEqual(summary["verdict"], "NOT_PASS")
        self.assertTrue(summary["incomplete"])
        self.assertEqual(summary["planned_pairs"], 8)
        self.assertEqual(summary["jointly_accepted_pairs"], 4)
        self.assertAlmostEqual(summary["diagnostic_gm"], 0.5)
        self.assertEqual(summary["release_coverage"], "BLOCKED_NOT_EVALUATED")
        self.assertEqual(len(summary["unavailable"]), 72)

    def test_complete_paired_medians_and_bootstrap(self):
        summary = self.harness.summarize(self.protocol, self.rows)
        self.assertEqual(summary["verdict"], "PASS")
        self.assertFalse(summary["incomplete"])
        self.assertEqual(summary["jointly_accepted_pairs"], 8)
        self.assertEqual(summary["scope"], "provisional_representative_gate")
        for backend in ("cpu", "cuda"):
            result = summary["backends"][backend]
            self.assertAlmostEqual(result["gm"], 0.5)
            self.assertEqual(result["ci95"], [0.5, 0.5])
            self.assertEqual(result["planned_pairs"], 4)
        self.assertEqual(summary["bootstrap_samples"], 10000)
        self.assertEqual(summary, self.harness.summarize(self.protocol, self.rows))

    def test_confidence_interval_is_reported_not_an_extra_gate(self):
        for row in self.rows:
            if row["solver"] == "candidate-cpu":
                row["solver_time_sec"] = (
                    1.8 if row["instance"] in ("afiro", "e226") else 2.2)
        summary = self.harness.summarize(self.protocol, self.rows)
        cpu = summary["backends"]["cpu"]
        self.assertLess(cpu["gm"], 1.05)
        self.assertGreater(cpu["ci95"][1], 1.05)
        self.assertEqual(cpu["verdict"], "PASS")
        self.assertEqual(summary["verdict"], "PASS")

    def test_rejected_or_missing_measurement_cannot_be_dropped(self):
        for reason in ("process_timeout", "kkt_pass=False", "missing"):
            with self.subTest(reason=reason):
                rows = copy.deepcopy(self.rows)
                target = next(row for row in rows if row["solver"] == "legacy-cuda"
                              and row["instance"] == "e226"
                              and row["kind"] == "measurement")
                if reason == "missing":
                    rows.remove(target)
                else:
                    target.update(accepted=False, rejection_reason=reason)
                summary = self.harness.summarize(self.protocol, rows)
                self.assertEqual(summary["verdict"], "NOT_PASS")
                self.assertTrue(summary["incomplete"])
                self.assertEqual(summary["jointly_accepted_pairs"], 7)
                self.assertEqual(len(summary["unavailable"]), 1)

    def test_one_backend_regression_cannot_hide_in_combined_gm(self):
        for row in self.rows:
            if row["solver"] == "candidate-cpu":
                row["solver_time_sec"] = 2.2
        summary = self.harness.summarize(self.protocol, self.rows)
        self.assertEqual(summary["verdict"], "NOT_PASS")
        self.assertFalse(summary["incomplete"])
        self.assertLess(summary["diagnostic_gm"], 1.05)
        self.assertAlmostEqual(summary["backends"]["cpu"]["gm"], 1.1)

    def test_completed_records_do_not_hide_failed_driver_process(self):
        for field, value in (("process_return_code", 2),
                             ("process_timed_out", True),
                             ("process_return_code", None)):
            with self.subTest(field=field, value=value):
                rows = copy.deepcopy(self.rows)
                for row in rows:
                    if row["solver"] == "candidate-cuda" and row["instance"] == "afiro":
                        row[field] = value
                summary = self.harness.summarize(self.protocol, rows)
                self.assertEqual(summary["verdict"], "NOT_PASS")
                self.assertEqual(len(summary["unavailable"]), 6)

    def test_twenty_percent_regressions_are_listed(self):
        for row in self.rows:
            if row["solver"] == "candidate-cpu" and row["instance"] == "afiro":
                row["solver_time_sec"] = 2.6
        summary = self.harness.summarize(self.protocol, self.rows)
        self.assertEqual(summary["regressions_over_20pct"], [
            {"backend": "cpu", "instance": "afiro", "ratio": 1.3}
        ])
        self.assertEqual(summary["verdict"], "NOT_PASS")

    def test_nonfinite_nonpositive_and_duplicate_rows_are_unavailable(self):
        for value in (math.nan, math.inf, 0.0, -1.0):
            with self.subTest(value=value):
                rows = copy.deepcopy(self.rows)
                rows[0]["solver_time_sec"] = value
                self.assertTrue(self.harness.summarize(self.protocol, rows)["incomplete"])
        self.rows.append(copy.deepcopy(self.rows[0]))
        self.assertTrue(self.harness.summarize(self.protocol, self.rows)["incomplete"])

    def test_warmups_are_excluded_from_medians_but_not_coverage(self):
        for row in self.rows:
            if row["kind"] == "warmup":
                row["solver_time_sec"] = 1000.0
        self.assertEqual(self.harness.summarize(self.protocol, self.rows)["diagnostic_gm"], 0.5)
        self.rows[0]["accepted"] = False
        self.assertTrue(self.harness.summarize(self.protocol, self.rows)["incomplete"])

    def test_outliers_do_not_replace_paired_medians(self):
        for row in self.rows:
            if row["solver"] == "candidate-cpu" and row["instance"] in ("afiro", "e226"):
                row["solver_time_sec"] = 2.0
            if row["solver"] == "candidate-cpu" and row["kind"] == "measurement" and row["index"] == 4:
                row["solver_time_sec"] = 1000.0
        summary = self.harness.summarize(self.protocol, self.rows)
        self.assertAlmostEqual(summary["backends"]["cpu"]["gm"], math.sqrt(0.5))
        self.assertEqual(summary["backends"]["cpu"]["ci95"], [0.5, 1.0])
        self.assertEqual(summary, self.harness.summarize(self.protocol, self.rows))

    def test_reduced_denominator_and_under_sampling_are_rejected(self):
        for field, value in (("repetitions", 4), ("warmups", 0),
                             ("instances", self.protocol["instances"][:-1]),
                             ("solvers", self.protocol["solvers"][:-1])):
            with self.subTest(field=field):
                protocol = copy.deepcopy(self.protocol)
                protocol[field] = value
                with self.assertRaises(ValueError):
                    self.harness.summarize(protocol, self.rows)


class PlannerTests(HarnessTest):
    def setUp(self):
        super().setUp()
        self.assertTrue(callable(getattr(self.harness, "make_plan", None)), "Phase8 planner is missing")
        scratch = RUNROOT / "performance-audit" / "scratch"
        scratch.mkdir(parents=True, exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(dir=scratch)
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        harness_root = Path("/home/yanyu/github/lp-benchmarks")
        self.protocol.update(
            benchmark_harness=str(harness_root / "run_benchmark.py"),
            reference_csv=str(RUNROOT.parent.parent / "phase7" / "summary.csv"),
            time_limit_sec=300, iteration_limit=2147483647, step_size="pid",
            timeout_grace_sec=30, driver_source=str(SCRIPT.with_name("HiPdlpPerformance.cpp")),
        )
        for instance in self.protocol["instances"]:
            suite = "netlib-feasible" if instance["name"] in ("afiro", "e226") else "energy-continuous"
            instance.update(manifest=str(harness_root / "manifests" / f"{suite}.csv"),
                            data_root="/home/yanyu/benchmark-data/lp/" +
                            ("netlib/feasible" if suite == "netlib-feasible" else "energy"))
        for solver in self.protocol["solvers"]:
            solver.update(source=str(RUNROOT / "frozen" / "source"),
                          build=str(self.root / "unbuilt" / solver["id"]),
                          library=str(self.root / "unbuilt" / solver["id"] / "lib" / "libhighs.so.1.15"),
                          driver=str(self.root / "drivers" / solver["id"] / "HiPdlpPerformance"),
                          compiler="/usr/bin/clang++", loader=None)
            solver["library_dir"] = str(Path(solver["library"]).parent)

    def test_planner_preserves_missing_production_jobs_without_execution(self):
        output = self.root / "not-created"
        plan = self.harness.make_plan(self.protocol, output)
        self.assertFalse(output.exists())
        self.assertEqual(len(plan["jobs"]), 24)
        self.assertTrue(plan["blockers"])
        pairs = [(job["solver"]["id"], job["instance"]["name"]) for job in plan["jobs"]]
        self.assertEqual(len(set(pairs)), 24)
        self.assertEqual(sum(solver.endswith("-cuda") for solver, name in pairs), 12)
        self.assertNotEqual(pairs, sorted(pairs))
        self.assertEqual(plan, self.harness.make_plan(self.protocol, output))
        other = copy.deepcopy(self.protocol)
        other["seed"] += 1
        self.assertNotEqual(plan["jobs"], self.harness.make_plan(other, output)["jobs"])
        for job in plan["jobs"]:
            command = job["command"]
            self.assertEqual(command[command.index("--warmups") + 1], "1")
            self.assertEqual(command[command.index("--repetitions") + 1], "5")
            self.assertEqual(command[command.index("--seed") + 1], "20260908")
            self.assertEqual(job["timeout_sec"], 330)

    def test_manifest_and_phase7_reference_values_are_loaded(self):
        instances = {row["name"]: row for row in self.harness.load_instances(self.protocol)}
        self.assertEqual(set(instances), set(NAMES))
        self.assertAlmostEqual(float(instances["afiro"]["reference_objective"]), -464.7531429)
        self.assertAlmostEqual(float(instances["pypsa_medium"]["reference_objective"]), 3147589.711596659)
        self.assertEqual(float(instances["pypsa-eur-elec_50_168h"]["reference_objective"]), 31872370698.000687)
        for row in instances.values():
            self.assertEqual(len(row["sha256"]), 64)
            self.assertTrue(Path(row["instance_path"]).is_file())

    def test_invalid_or_ambiguous_energy_reference_is_not_silently_used(self):
        fields = ("scope", "solver", "instance", "status", "accepted", "primal_objective")
        for defect in ("unaccepted", "nonfinite", "duplicate", "missing"):
            with self.subTest(defect=defect):
                reference = self.root / f"{defect}.csv"
                rows = [dict(zip(fields, ("energy-full", "highs-ipm", name, "OPTIMAL", "True", "12")))
                        for name in ("pypsa_medium", "pypsa-eur-elec_50_168h")]
                if defect == "unaccepted":
                    rows[0]["accepted"] = "False"
                elif defect == "nonfinite":
                    rows[0]["primal_objective"] = "nan"
                elif defect == "duplicate":
                    rows.append(dict(rows[0]))
                else:
                    rows.pop()
                with reference.open("x", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=fields)
                    writer.writeheader()
                    writer.writerows(rows)
                protocol = dict(self.protocol, reference_csv=str(reference))
                with self.assertRaises(ValueError):
                    self.harness.load_instances(protocol)

    def test_compile_commands_pair_headers_library_and_legacy_macro(self):
        for solver in self.protocol["solvers"]:
            command = self.harness.compile_command(self.protocol, solver)
            self.assertEqual(command[0], solver["compiler"])
            self.assertIn("-I" + solver["build"], command)
            self.assertIn("-I" + solver["source"] + "/highs", command)
            self.assertIn(solver["library"], command)
            self.assertEqual("-DHIPDLP_PERFORMANCE_LEGACY" in command, solver["role"] == "legacy")
            self.assertEqual(command[command.index("-o") + 1], solver["driver"])

    def test_json_artifact_collision_and_nonfinite_output_are_refused(self):
        output = self.root / "protected.json"
        self.harness.write_json(output, {"seed": 20260908})
        original = output.read_bytes()
        with self.assertRaises(FileExistsError):
            self.harness.write_json(output, {"seed": 0})
        self.assertEqual(output.read_bytes(), original)
        invalid = self.root / "nonfinite.json"
        with self.assertRaises(ValueError):
            self.harness.write_json(invalid, {"duration": math.nan})
        self.assertFalse(invalid.exists())

    def test_run_requires_explicit_authorization(self):
        protocol_file = self.root / "protocol.json"
        protocol_file.write_text(json.dumps(self.protocol))
        output = self.root / "must-not-run"
        result = subprocess.run([sys.executable, "-B", str(SCRIPT), "run",
                                 "--protocol", str(protocol_file), "--output-dir", str(output)],
                                capture_output=True, text=True, timeout=5)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--execute", result.stderr)
        self.assertFalse(output.exists())

    def test_failed_correctness_evidence_stops_before_creating_run(self):
        protocol_file = self.root / "protocol.json"
        protocol_file.write_text(json.dumps(self.protocol))
        evidence = self.root / "failed-correctness.json"
        evidence.write_text(json.dumps({"planned": 12, "completed": 12,
                                        "all_accepted": False, "rows": []}))
        output = self.root / "must-not-start"
        result = subprocess.run(
            [sys.executable, "-B", str(SCRIPT), "run", "--protocol", str(protocol_file),
             "--output-dir", str(output), "--execute", "--correctness-evidence", str(evidence)],
            capture_output=True, text=True, timeout=5)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("correctness evidence", result.stderr)
        self.assertFalse(output.exists())


class CorrectnessEvidenceTests(HarnessTest):
    def setUp(self):
        super().setUp()
        self.assertTrue(callable(getattr(self.harness, "validate_correctness_evidence", None)))
        rows = []
        for backend in ("cpu", "cuda"):
            for strategy, names in (("pid", NAMES), ("adaptive", NAMES[:2])):
                for name in names:
                    rows.append(dict(self.harness.FIXED_OPTIONS,
                        solver="candidate-" + backend, instance=name,
                        seed=self.protocol["seed"],
                        requested_step_size=strategy,
                        effective_step_size="fixed" if strategy == "adaptive" else "pid",
                        accepted=True, status="OPTIMAL", kkt_pass=True, kkt_input_valid=True,
                        relative_primal_residual=0.0, relative_dual_residual=0.0,
                        relative_objective_gap=0.0, reference_objective_relative_error=0.0,
                        process_return_code=0, process_timed_out=False))
        self.evidence = dict(planned=12, completed=12, all_accepted=True, rows=rows)

    def test_complete_evidence_is_valid(self):
        self.harness.validate_correctness_evidence(self.protocol, self.evidence)

    def test_invalid_record_or_coverage_is_not_a_prerequisite(self):
        for defect in ("missing", "duplicate", "rejected", "metrics", "process", "strategy", "tolerance"):
            with self.subTest(defect=defect):
                evidence = copy.deepcopy(self.evidence)
                if defect == "missing":
                    evidence["rows"].pop()
                elif defect == "duplicate":
                    evidence["rows"][-1] = dict(evidence["rows"][0])
                else:
                    field, value = {
                        "rejected": ("accepted", False),
                        "metrics": ("relative_dual_residual", 0.0002),
                        "process": ("process_return_code", 2),
                        "strategy": ("effective_step_size", "fixed"),
                        "tolerance": ("pdlp_optimality_tolerance", 0.0002),
                    }[defect]
                    evidence["rows"][0][field] = value
                with self.assertRaisesRegex(ValueError, "correctness evidence"):
                    self.harness.validate_correctness_evidence(self.protocol, evidence)


class AcceptanceTests(HarnessTest):
    def setUp(self):
        super().setUp()
        self.assertTrue(callable(getattr(self.harness, "accept_result", None)), "External acceptance is missing")
        self.protocol["benchmark_harness"] = "/home/yanyu/github/lp-benchmarks/run_benchmark.py"
        self.reference = {"reference_objective": "100"}
        self.driver = {"status": "OPTIMAL", "return_status": "ok", "iterations": 40,
                       "solver_time_sec": 1.0, "process_return_code": 0,
                       "process_timed_out": False}
        self.verification = {
            "status": "OPTIMAL", "return_status": "ok", "iterations": 40,
            "return_code": 0, "timed_out": False, "primal_objective": 100.0,
            "kkt_input_valid": True, "kkt_pass": True,
            "relative_primal_residual": 0.0, "relative_dual_residual": 0.0,
            "relative_objective_gap": 0.0,
        }

    def test_reference_error_uses_plan_one_per_mille_gate(self):
        verification = dict(self.verification, primal_objective=100.05)
        result = self.harness.accept_result(
            self.protocol, self.reference, self.driver, verification)
        self.assertTrue(result["accepted"], result["rejection_reason"])

    def test_verifier_success_does_not_hide_failed_driver(self):
        for field, value in (("process_return_code", 2),
                             ("process_timed_out", True),
                             ("process_return_code", None)):
            with self.subTest(field=field, value=value):
                driver = dict(self.driver, **{field: value})
                result = self.harness.accept_result(
                    self.protocol, self.reference, driver, self.verification)
                self.assertFalse(result["accepted"])
                self.assertIn("driver_process", result["rejection_reason"])

    def test_independent_kkt_and_reference_are_both_required(self):
        result = self.harness.accept_result(self.protocol, self.reference, self.driver, self.verification)
        self.assertTrue(result["accepted"])
        for field, value in (("kkt_pass", False), ("primal_objective", None),
                             ("primal_objective", 102.0), ("relative_dual_residual", 0.0002),
                             ("timed_out", True), ("return_code", 3)):
            with self.subTest(field=field, value=value):
                verification = dict(self.verification, **{field: value})
                result = self.harness.accept_result(self.protocol, self.reference, self.driver, verification)
                self.assertFalse(result["accepted"])
                self.assertTrue(result["rejection_reason"])

    def test_verifier_cannot_upgrade_solver_status_or_invalid_duration(self):
        for field, value in (("status", "TIME_LIMIT"), ("return_status", "error"),
                             ("solver_time_sec", 0), ("solver_time_sec", math.nan)):
            with self.subTest(field=field, value=value):
                driver = dict(self.driver, **{field: value})
                result = self.harness.accept_result(self.protocol, self.reference, driver, self.verification)
                self.assertFalse(result["accepted"])


class ProductionTests(HarnessTest):
    def test_testing_enabled_audit_is_not_a_production_baseline(self):
        self.assertTrue(callable(getattr(self.harness, "production_inputs", None)), "Production audit is missing")
        build = Path("/home/yanyu/builds/hipdlp2-phase8-audit-cpu-20260908T000926Z")
        solver = {"id": "candidate-cpu", "role": "candidate", "backend": "cpu",
                  "build": str(build), "source": "/home/yanyu/github/hipdlp2",
                  "library": str(build / "lib" / "libhighs.so.1.15")}
        with self.assertRaisesRegex(ValueError, "BUILD_TESTING"):
            self.harness.production_inputs(solver)


class ProcessTests(HarnessTest):
    def setUp(self):
        super().setUp()
        self.assertTrue(callable(getattr(self.harness, "run_driver", None)), "Phase8 timeout monitor is missing")
        scratch = RUNROOT / "performance-audit" / "scratch"
        scratch.mkdir(parents=True, exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(dir=scratch)
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def test_timeout_kills_only_owned_child_and_preserves_flushed_record(self):
        code = ("import json, signal; from pathlib import Path; "
                "print(json.dumps({'event':'begin','sequence':0}), flush=True); "
                "Path('warmup-000.json').write_text('{\"completed\":true}'); "
                "print(json.dumps({'event':'begin','sequence':1}), flush=True); signal.pause()")
        result = self.harness.run_driver([sys.executable, "-B", "-c", code], self.root,
                                         0.15, dict(os.environ))
        self.assertTrue(result["timed_out"])
        self.assertLess(result["return_code"], 0)
        self.assertGreater(result["process_wall_sec"], 0.1)
        self.assertEqual(json.loads((self.root / "warmup-000.json").read_text()), {"completed": True})
        self.assertFalse((self.root / "measurement-000.json").exists())
        self.assertIn('"sequence": 1', (self.root / "stdout.log").read_text())

    def test_successful_child_exit_and_stderr_are_preserved(self):
        code = "import sys; print('ordinary output'); print('diagnostic', file=sys.stderr)"
        result = self.harness.run_driver([sys.executable, "-B", "-c", code], self.root,
                                         1, dict(os.environ))
        self.assertFalse(result["timed_out"])
        self.assertEqual(result["return_code"], 0)
        self.assertGreater(result["process_wall_sec"], 0)
        self.assertIn("ordinary output", (self.root / "stdout.log").read_text())
        self.assertIn("diagnostic", (self.root / "stderr.log").read_text())


@unittest.skipUnless(os.environ.get("HIPDLP_PERFORMANCE_CPU_DRIVERS"),
                     "Set HIPDLP_PERFORMANCE_CPU_DRIVERS to opt into standalone CPU compilation")
class CompileHelperTests(HarnessTest):
    def test_compile_only_helper_freezes_and_checks_provenance(self):
        helper_file = RUNROOT / "performance_compile_drivers.py"
        self.assertTrue(helper_file.is_file(), "Compile-only provenance helper is missing")
        spec = importlib.util.spec_from_file_location("phase8_compile_helper", helper_file)
        helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(helper)
        scratch = RUNROOT / "performance-audit" / "scratch"
        scratch.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=scratch) as temporary:
            root = Path(temporary)
            frozen = RUNROOT / "frozen"
            solver = {"id": "legacy-cpu", "role": "legacy", "backend": "cpu",
                      "source": str(frozen / "legacy-source"),
                      "build": str(frozen / "legacy-build-cpu"),
                      "library": str(frozen / "legacy-runtime-cpu/lib/libhighs.so.1.15"),
                      "library_dir": str(frozen / "legacy-runtime-cpu/lib"),
                      "loader": str(frozen / "legacy-runtime-cpu/lib/ld-linux-aarch64.so.1"),
                      "compiler": "/usr/bin/clang++", "driver": str(root / "build" / "HiPdlpPerformance")}
            verifier = frozen / "legacy-verifier"
            protocol = dict(self.protocol, solvers=[solver], instances=[],
                            driver_source=str(SCRIPT.with_name("HiPdlpPerformance.cpp")),
                            harness_script=str(SCRIPT),
                            benchmark_harness="/home/yanyu/github/lp-benchmarks/run_benchmark.py",
                            reference_csv=str(RUNROOT.parent.parent / "phase7/summary.csv"),
                            verifier={"binary": str(verifier / "bin/external_lp_solution_verifier"),
                                      "library_dir": str(verifier / "lib"),
                                      "loader": str(verifier / "lib/ld-linux-aarch64.so.1")})
            protocol_file = root / "protocol.json"
            protocol_file.write_text(json.dumps(protocol))
            proof = helper.build_one(protocol, solver)
            self.assertTrue(Path(solver["driver"]).is_file())
            self.assertIn("-DHIPDLP_PERFORMANCE_LEGACY", proof["compile_command"])
            self.assertIn("libhighs.so.1.15", proof["ldd"])
            self.assertNotIn("libcuda", proof["ldd"])
            self.assertEqual(proof["file_sha256"][solver["driver"]], self.harness.file_hash(solver["driver"]))
            hashes = self.harness.verify_provenance(protocol, protocol_file)
            self.assertEqual(hashes[solver["driver"]], proof["file_sha256"][solver["driver"]])
            with self.assertRaises(FileExistsError):
                helper.build_one(protocol, solver)
            snapshot = Path(solver["driver"]).with_name("HiPdlpPerformance.cpp")
            with snapshot.open("a") as stream:
                stream.write("\n")
            with self.assertRaises(ValueError):
                self.harness.verify_provenance(protocol, protocol_file)


@unittest.skipUnless(os.environ.get("HIPDLP_PERFORMANCE_CPU_DRIVERS"),
                     "Set HIPDLP_PERFORMANCE_CPU_DRIVERS to opt into CPU afiro checks")
class DriverTests(HarnessTest):
    def setUp(self):
        super().setUp()
        self.drivers = [Path(value) for value in
                        os.environ["HIPDLP_PERFORMANCE_CPU_DRIVERS"].split(os.pathsep)]
        for driver in self.drivers:
            self.assertTrue(driver.is_file(), f"CPU audit driver is missing: {driver}")
        scratch = RUNROOT / "performance-audit" / "scratch"
        scratch.mkdir(parents=True, exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(dir=scratch)
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.environment = dict(os.environ)
        for key in ("LD_LIBRARY_PATH", "LD_PRELOAD", "LD_AUDIT"):
            self.environment.pop(key, None)
        self.environment["CUDA_VISIBLE_DEVICES"] = ""
        for driver in self.drivers:
            loader = subprocess.run(["ldd", str(driver)], env=self.environment,
                                    capture_output=True, text=True, check=True)
            self.assertNotIn("libcuda", loader.stdout.lower())
            self.assertNotIn("libcublas", loader.stdout.lower())
            self.assertNotIn("libcusparse", loader.stdout.lower())

    def command(self, driver, directory, instance=AFIRO, strategy="pid"):
        return [str(driver), "--instance", str(instance), "--backend", "cpu",
                "--step-size", strategy, "--time-limit", "5", "--iteration-limit", "250000",
                "--output-dir", str(directory), "--warmups", "1", "--repetitions", "2",
                "--seed", "20260908"]

    def invoke(self, command, directory):
        return subprocess.run(command, env=self.environment, cwd=directory,
                              capture_output=True, text=True, timeout=35)

    def test_cpu_same_process_results_and_raw_solution(self):
        for driver in self.drivers:
            with self.subTest(driver=driver):
                directory = self.root / driver.name
                directory.mkdir()
                result = self.invoke(self.command(driver, directory), directory)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                files = [directory / "warmup-000.json", directory / "measurement-000.json",
                         directory / "measurement-001.json"]
                rows = [json.loads(file.read_text()) for file in files]
                self.assertEqual([row["cold"] for row in rows], [True, False, False])
                self.assertEqual([row["sequence"] for row in rows], [0, 1, 2])
                self.assertEqual(len({row["pid"] for row in rows}), 1)
                self.assertEqual(len({row["read_time_sec"] for row in rows}), 1)
                self.assertEqual(len(list(directory.glob("*.sol"))), 3)
                for row in rows:
                    self.assertEqual(row["status"], "OPTIMAL")
                    self.assertEqual(row["backend"], "cpu")
                    self.assertEqual(row["effective_step_size"], "pid")
                    self.assertEqual(row["requested_step_size_enum"], 3)
                    self.assertEqual(row["kkt_tolerance"], 1e-4)
                    self.assertEqual(row["pdlp_optimality_tolerance"], 1e-4)
                    self.assertEqual(row["primal_feasibility_tolerance"], 1e-4)
                    self.assertEqual(row["dual_feasibility_tolerance"], 1e-4)
                    self.assertEqual(row["presolve"], "off")
                    self.assertEqual(row["threads"], 1)
                    self.assertGreater(row["solver_time_sec"], 0)
                    self.assertLessEqual(row["solver_time_sec"], row["repeat_wall_sec"])
                    solution = directory / row["solution_file"]
                    self.assertIn("# Dual solution values", solution.read_text())
                    verifier_root = RUNROOT / "frozen" / "legacy-verifier"
                    output = directory / f"verified-{row['sequence']}.json"
                    command = [str(verifier_root / "lib" / "ld-linux-aarch64.so.1"),
                               "--library-path", str(verifier_root / "lib"),
                               str(verifier_root / "bin" / "external_lp_solution_verifier"),
                               "--instance", str(AFIRO), "--solution", str(solution),
                               "--solver", driver.name, "--status", row["status"],
                               "--solver-time", str(row["solver_time_sec"]),
                               "--iterations", str(row["iterations"]), "--tolerance", "0.0001",
                               "--output", str(output)]
                    verified = self.invoke(command, directory)
                    self.assertEqual(verified.returncode, 0, verified.stdout + verified.stderr)
                    payload = json.loads(output.read_text())
                    self.assertTrue(payload["kkt_input_valid"])
                    self.assertTrue(payload["kkt_pass"], payload)
                    self.assertLess(abs(payload["primal_objective"] + 464.7531429) / 465.7531429, 1e-4)

    def test_orchestrator_verifies_cpu_afiro_records(self):
        self.assertTrue(callable(getattr(self.harness, "verify_job", None)), "Per-record verification is missing")
        verifier_root = RUNROOT / "frozen" / "legacy-verifier"
        self.protocol.update(
            benchmark_harness="/home/yanyu/github/lp-benchmarks/run_benchmark.py",
            repetitions=2, time_limit_sec=5, iteration_limit=250000, step_size="pid",
            verifier={"binary": str(verifier_root / "bin" / "external_lp_solution_verifier"),
                      "library_dir": str(verifier_root / "lib"),
                      "loader": str(verifier_root / "lib" / "ld-linux-aarch64.so.1")},
        )
        for driver in self.drivers:
            directory = self.root / driver.name
            directory.mkdir()
            command = self.command(driver, directory)
            process = self.harness.run_driver(command, directory, 35, self.environment)
            job = {"directory": str(directory), "solver": {
                "id": driver.name, "role": "legacy" if driver.name.startswith("legacy") else "candidate",
                "backend": "cpu"}, "instance": {
                    "name": "afiro", "instance_path": str(AFIRO), "reference_objective": "-464.7531429"}}
            rows = self.harness.verify_job(self.protocol, job, process)
            self.assertEqual(len(rows), 3)
            self.assertTrue(all(row["accepted"] for row in rows), rows)
            self.assertTrue(all(row["instance"] == "afiro" for row in rows))
            self.assertTrue(all(row["solver"] == driver.name for row in rows))
            self.assertEqual(len(list(directory.glob("verified-*.json"))), 3)

    def test_strategy_translation_and_long_escaped_model_path(self):
        model_directory = self.root / 'model "quoted"\\with\ncontrol'
        for index in range(12):
            model_directory /= f"{index:02d}-" + "long" * 22
        model_directory.mkdir(parents=True)
        model = model_directory / "afiro.mps"
        model.symlink_to(AFIRO)
        for driver in self.drivers:
            for strategy, enum in (("fixed", 0), ("adaptive", 1)):
                with self.subTest(driver=driver, strategy=strategy):
                    directory = self.root / f"{driver.name}-{strategy}"
                    directory.mkdir()
                    command = self.command(driver, directory, model, strategy)
                    command[command.index("--warmups") + 1] = "0"
                    command[command.index("--repetitions") + 1] = "1"
                    result = self.invoke(command, directory)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    row = json.loads((directory / "measurement-000.json").read_text())
                    expected = "pid" if row["implementation"] == "legacy" and strategy == "adaptive" else "fixed"
                    self.assertEqual(row["effective_step_size"], expected)
                    self.assertEqual(row["requested_step_size_enum"], enum)
                    self.assertEqual(row["instance"], str(model))
                    self.assertEqual(row["scaling"], "legacy-native-default" if
                                     row["implementation"] == "legacy" else "auto")

    def test_invalid_options_and_backend_are_rejected_without_solving(self):
        for driver in self.drivers:
            for flag, value in (("--backend", "cuda"), ("--backend", "bogus"),
                                ("--step-size", "bogus"), ("--time-limit", "nan"),
                                ("--iteration-limit", "-1"), ("--seed", "999999999999999999999")):
                with self.subTest(driver=driver, flag=flag, value=value):
                    with tempfile.TemporaryDirectory(dir=self.root) as temporary:
                        directory = Path(temporary)
                        command = self.command(driver, directory)
                        command[command.index(flag) + 1] = value
                        result = self.invoke(command, directory)
                        self.assertNotEqual(result.returncode, 0)
                        self.assertEqual(list(directory.glob("*.sol")), [])

    def test_output_collision_preserves_completed_files(self):
        for driver in self.drivers:
            directory = self.root / driver.name
            directory.mkdir()
            command = self.command(driver, directory)
            command[command.index("--warmups") + 1] = "0"
            command[command.index("--repetitions") + 1] = "1"
            result = self.invoke(command, directory)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            before = {file.name: file.read_bytes() for file in directory.iterdir()}
            result = self.invoke(command, directory)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(before, {file.name: file.read_bytes() for file in directory.iterdir()})


if __name__ == "__main__":
    unittest.main()