import json
from html.parser import HTMLParser
from pathlib import Path
import re
import unittest
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs" / "hipdlp-architecture.html"
SOLVER_HEADER = ROOT / "highs" / "pdlp" / "hipdlp" / "algorithm" / "solver.hpp"
SYSTEM_NODE_IDS = {
    "node-wrapper",
    "node-api",
    "node-model",
    "node-algorithm",
    "node-backend",
    "node-postsolve",
    "node-kkt",
}


class ArchitecturePageParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()
        self.tabs = []
        self.remote_assets = []
        self.links = []
        self.data_chunks = []
        self.in_data = False
        self.elements = {}
        self.lane_buttons = []
        self.backend_lane_cards = []
        self.shared_operation_cards = []
        self.style_chunks = []
        self.in_style = False
        self.current_node = None
        self.detail_node = None
        self.node_text = {node_id: [] for node_id in SYSTEM_NODE_IDS}
        self.node_detail_count = {node_id: 0 for node_id in SYSTEM_NODE_IDS}
        self.node_detail_text = {node_id: [] for node_id in SYSTEM_NODE_IDS}
        self.node_detail_links = {node_id: [] for node_id in SYSTEM_NODE_IDS}

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        element_id = attributes.get("id")
        if element_id:
            self.ids.add(element_id)
            self.elements[element_id] = {"tag": tag, "attributes": attributes}
        if tag == "article" and element_id in SYSTEM_NODE_IDS:
            self.current_node = element_id
        classes = set(attributes.get("class", "").split())
        if tag == "button" and "data-lane" in attributes:
            self.lane_buttons.append(attributes)
        if tag == "div" and "lane-card" in classes:
            self.backend_lane_cards.append(attributes)
        if tag == "li" and "operation-card" in classes:
            self.shared_operation_cards.append(attributes)
        if (
            tag == "details"
            and self.current_node is not None
            and "component-detail" in classes
        ):
            self.detail_node = self.current_node
            self.node_detail_count[self.current_node] += 1
        if tag == "a" and self.detail_node is not None:
            self.node_detail_links[self.detail_node].append(
                attributes.get("href", "")
            )
        if attributes.get("role") == "tab":
            self.tabs.append(attributes)
        for name in ("href", "src"):
            value = attributes.get(name, "")
            if value.startswith(("http://", "https://", "//")):
                self.remote_assets.append(value)
            if name == "href" and value:
                self.links.append(value)
        if tag == "script" and element_id == "architecture-data":
            self.in_data = True
        if tag == "style":
            self.in_style = True

    def handle_endtag(self, tag):
        if tag == "script" and self.in_data:
            self.in_data = False
        if tag == "style":
            self.in_style = False
        if tag == "details" and self.detail_node is not None:
            self.detail_node = None
        if tag == "article" and self.current_node is not None:
            self.current_node = None

    def handle_data(self, data):
        if self.in_data:
            self.data_chunks.append(data)
        if self.in_style:
            self.style_chunks.append(data)
        if self.current_node is not None:
            self.node_text[self.current_node].append(data)
        if self.detail_node is not None:
            self.node_detail_text[self.detail_node].append(data)


class ArchitecturePageTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(PAGE.is_file(), "Create docs/hipdlp-architecture.html")
        self.source = PAGE.read_text(encoding="utf-8")
        self.parser = ArchitecturePageParser()
        self.parser.feed(self.source)
        self.data = json.loads("".join(self.parser.data_chunks))

    @staticmethod
    def normalized(chunks):
        return " ".join("".join(chunks).split())

    def test_page_is_offline_and_has_five_views(self):
        self.assertEqual(self.parser.remote_assets, [])
        self.assertEqual(len(self.parser.tabs), 5)
        self.assertTrue({
            "system-map",
            "data-residency",
            "iteration-correctness",
            "bottlenecks",
            "learning-review",
        } <= self.parser.ids)
        self.assertEqual(self.data["snapshot"], {
            "branch": "feature/hipdlp-rewrite",
            "commit": "49b85f4f3",
            "generated": "2026-09-23",
            "workingTree": "dirty",
        })

    def test_system_map_has_the_real_solve_boundaries(self):
        self.assertTrue(SYSTEM_NODE_IDS <= self.parser.ids)
        for label in (
            "solveLpHiPdlp",
            "hipdlp::api::solve",
            "model::buildProblem",
            "algorithm::solve",
            "backend::createBackend",
            "model::postsolve",
            "verification::verifyKkt",
        ):
            self.assertIn(label, self.source)

    def test_data_model_matches_current_problem_and_iteration_storage(self):
        structures = {item["id"]: item for item in self.data["structures"]}
        self.assertEqual(
            set(structures),
            {"highs-lp", "original-model", "problem-data", "handles",
             "cpu-resources", "cuda-resources", "host-iterate",
             "highs-solution"},
        )
        inventory = self.data["vectorInventory"]
        self.assertEqual(inventory["columnSized"], 8)
        self.assertEqual(inventory["rowSized"], 8)
        self.assertEqual(inventory["workspaceColumnSized"], 1)
        self.assertEqual(inventory["workspaceRowSized"], 1)
        self.assertEqual(
            inventory["candidateDownloadDoubles"], "2*(rows+cols)"
        )
        self.assertEqual(inventory["columnVectors"], [
            "current_primal",
            "next_primal",
            "primal_anchor",
            "primal_candidate",
            "reflected_primal",
            "current_aty",
            "candidate_aty",
            "aty_delta",
        ])
        self.assertEqual(inventory["rowVectors"], [
            "current_dual",
            "next_dual",
            "dual_anchor",
            "dual_candidate",
            "reflected_dual",
            "current_ax",
            "reflected_ax",
            "candidate_ax",
        ])
        self.assertEqual(
            self.data["problemData"]["transformationMapFields"],
            [
                "original_to_transformed_row",
                "transformed_to_original_row",
                "original_row_sign",
                "slack_col_for_original_row",
            ],
        )
        self.assertEqual(
            self.data["cudaResourceCategories"],
            {
                "CudaProblemResource": {
                    "spmvScratch": ["spmv_buffer"],
                },
                "CudaBackend::Impl": {
                    "executionHandles": ["stream", "cublas", "cusparse"],
                    "timing": ["timing_events"],
                    "errorState": ["error_flag"],
                    "reductionScratch": [
                        "reduction_result",
                        "difference_scratch",
                    ],
                },
            },
        )
        for token in (
            "TransformationMetadata",
            "ScalingMetadata",
            "UnscaledNorms",
            "ProblemHandle",
            "VectorHandle",
            "WorkspaceHandle",
            "CommandBufferHandle",
            "PinnedDownloadBuffer",
        ):
            self.assertIn(token, self.source)

    def test_backend_card_states_broad_auto_fallback_as_a_review_risk(self):
        card_text = self.normalized(self.parser.node_text["node-backend"])
        detail_text = self.normalized(
            self.parser.node_detail_text["node-backend"]
        )
        for fact in (
            "源码事实",
            "createCudaBackend",
            "任何非 OK 错误",
            "显式 kCuda",
            "返回 createCudaBackend 的错误",
            "不进入 CPU fallback",
            "不代表期望行为",
        ):
            with self.subTest(fact=fact):
                self.assertIn(fact, card_text)
        review_text = detail_text.partition("Review：")[2]
        for risk in (
            "风险",
            "吞掉",
            "kBackendUnavailable",
            "kAllocationFailure",
            "kBackendFailure",
        ):
            with self.subTest(risk=risk):
                self.assertIn(risk, review_text)

    def test_flow_uses_control_and_correctness_visual_roles(self):
        for decision_id in ("decision-api-route", "decision-backend-kind"):
            with self.subTest(decision_id=decision_id):
                self.assertIn(decision_id, self.parser.elements)
                classes = set(
                    self.parser.elements[decision_id]["attributes"]
                    .get("class", "")
                    .split()
                )
                self.assertIn("control-fill", classes)
                self.assertNotIn("host-fill", classes)
                self.assertNotIn("correctness-fill", classes)

        self.assertIn("gate-kkt", self.parser.elements)
        kkt_classes = set(
            self.parser.elements["gate-kkt"]["attributes"]
            .get("class", "")
            .split()
        )
        self.assertIn("correctness-fill", kkt_classes)
        self.assertNotIn("host-fill", kkt_classes)

        styles = self.normalized(self.parser.style_chunks)
        self.assertRegex(
            styles,
            r"\.control-fill\s*\{[^}]*fill:\s*var\(--control-soft\)"
            r"[^}]*stroke:\s*var\(--control\)[^}]*\}",
        )
        self.assertRegex(
            styles,
            r"\.correctness-fill\s*\{[^}]*fill:\s*var\(--correctness-soft\)"
            r"[^}]*stroke:\s*var\(--correctness\)[^}]*\}",
        )

    def test_every_system_node_has_complete_review_details(self):
        required_labels = (
            "所有权",
            "输入",
            "输出",
            "执行位置",
            "失败",
            "源码",
            "Review：",
        )
        for node_id in sorted(SYSTEM_NODE_IDS):
            with self.subTest(node_id=node_id):
                self.assertEqual(self.parser.node_detail_count[node_id], 1)
                detail_text = self.normalized(
                    self.parser.node_detail_text[node_id]
                )
                for label in required_labels:
                    self.assertIn(label, detail_text)
                review_question = detail_text.partition("Review：")[2]
                self.assertGreaterEqual(len(review_question), 12)
                self.assertIn("？", review_question)
                links = self.parser.node_detail_links[node_id]
                self.assertEqual(len(links), 1)
                self.assertTrue(links[0].startswith("../highs/"))

    def test_all_local_source_links_resolve(self):
        source_links = [
            link for link in self.parser.links
            if link.startswith("../highs/") or link.startswith("../check/")
        ]
        self.assertGreaterEqual(len(source_links), 7)
        missing = []
        for link in source_links:
            path = unquote(urlsplit(link).path)
            if not (PAGE.parent / path).resolve().is_file():
                missing.append(link)
        self.assertEqual(missing, [])

    def test_iteration_and_cuda_constants_match_source(self):
        iteration = self.data["iteration"]
        self.assertEqual(iteration["spmvPerStep"], 2)
        self.assertEqual(iteration["convergenceCadence"], 40)
        self.assertEqual(iteration["recordedSteps"], 39)
        self.assertEqual(iteration["parameterDoubles"], 117)
        self.assertEqual(iteration["parameterBytes"], 936)
        self.assertEqual(iteration["checkPacketBytes"], 112)

        kernels = {item["name"]: item for item in self.data["kernels"]}
        self.assertEqual(kernels["projectedPrimalKernel"]["threads"], 256)
        self.assertEqual(kernels["projectedPrimalKernel"]["maxBlocks"], 65535)
        self.assertEqual(kernels["checkPartialsKernel"]["threads"], 256)
        self.assertEqual(kernels["checkPartialsKernel"]["maxBlocks"], 1024)
        self.assertEqual(kernels["checkPacketKernel"]["blocks"], 1)
        self.assertEqual(kernels["cusparseSpMV"]["grid"], "vendor-controlled")
        self.assertEqual(kernels["checkPacketKernel"]["dynamicSharedBytes"], 0)

    def test_replay_semantics_match_solver_scheduling(self):
        self.assertEqual(
            self.data["iteration"]["replaySemantics"],
            {
                "initialStep": "direct",
                "firstPostRestartStep": "direct",
                "shortTail": "direct",
                "submittedBlockSteps": 39,
                "hostTimeLimitChecks": "before-and-after-block",
                "hostInterruptible": False,
            },
        )
        for phrase in (
            "The initial step and the first post-restart step execute directly.",
            "Short tails with fewer than 39 remaining steps execute directly.",
            "a submitted 39-step block cannot be interrupted by host "
            "time-limit checks",
        ):
            with self.subTest(phrase=phrase):
                self.assertIn(phrase, self.source)

    def test_lane_controls_filter_only_backend_cards_and_default_to_compare(self):
        self.assertEqual(
            [
                (button["data-lane"], button["aria-pressed"])
                for button in self.parser.lane_buttons
            ],
            [("cpu", "false"), ("cuda", "false"), ("compare", "true")],
        )
        self.assertEqual(len(self.parser.backend_lane_cards), 10)
        self.assertEqual(
            {
                card.get("data-backend")
                for card in self.parser.backend_lane_cards
            },
            {"cpu", "cuda"},
        )
        self.assertTrue(
            all("hidden" not in card for card in self.parser.backend_lane_cards)
        )
        self.assertEqual(len(self.parser.shared_operation_cards), 5)
        self.assertTrue(
            all(
                "data-backend" not in card and "hidden" not in card
                for card in self.parser.shared_operation_cards
            )
        )
        for behavior in (
            'function activateLane(lane)',
            'document.querySelectorAll(".lane-card[data-backend]")',
            'button.setAttribute("aria-pressed", String(isActive))',
            'card.hidden = lane !== "compare" && card.dataset.backend !== lane',
            'button.addEventListener("click", () => activateLane(',
            'activateLane("compare")',
        ):
            with self.subTest(behavior=behavior):
                self.assertIn(behavior, self.source)

    def test_iteration_facts_stack_below_480px(self):
        styles = self.normalized(self.parser.style_chunks)
        self.assertIn(
            "@media (max-width: 480px) { .field-groups, "
            ".iteration-facts { grid-template-columns: 1fr; } }",
            styles,
        )
        self.assertRegex(
            styles,
            r"\.resource-card \{[^}]*min-width: 0;"
            r"[^}]*overflow-wrap: anywhere;",
        )
        self.assertRegex(
            styles,
            r"\.operation-card \{[^}]*min-width: 0;",
        )
        self.assertIn(
            ".operation-card code { overflow-wrap: anywhere; }",
            styles,
        )

    def test_iteration_view_names_real_backend_and_kernel_symbols(self):
        for token in (
            "enqueueHalpernPdhgStep",
            "fusedProjectedPrimalUpdateAndReflect",
            "multiplyAndProjectedDualUpdateAndReflect",
            "fusedVectorBlend",
            "multiplyTranspose",
            "executeSpmv",
            "cusparseSpMV",
            "projectedPrimalKernel&lt;true&gt;",
            "projectedDualKernel&lt;true&gt;",
            "vectorPairBlendKernel",
            "checkPartialsKernel",
            "checkPacketKernel",
            "cudaGraphLaunch",
        ):
            self.assertIn(token, self.source)

    def test_error_rail_covers_explicit_failure_surfaces(self):
        for phrase in (
            "Invalid option / argument",
            "Backend unavailable / not ready",
            "Allocation failure",
            "Invalid model",
            "Dimension / alias violation",
            "Non-finite arithmetic",
            "CUDA / cuBLAS / cuSPARSE failure",
            "Graph capture / instantiate / replay failure",
            "Invalid resident check",
            "Postsolve / original-model KKT failure",
        ):
            self.assertIn(phrase, self.source)

    def test_performance_evidence_retains_scope_and_limits(self):
        evidence = {item["id"]: item for item in self.data["evidence"]}
        self.assertEqual(evidence["phase8-gate"]["date"], "2026-09-08")
        self.assertEqual(
            evidence["cuda-fusion-replay"]["date"], "2026-09-22"
        )
        self.assertEqual(
            evidence["cuda-fusion-replay"]["reductionPercent"],
            [41.5, 52.4],
        )
        self.assertFalse(evidence["cuda-fusion-replay"]["releaseGate"])
        self.assertEqual(evidence["large20"]["acceptedPairs"], 0)
        for phrase in (
            "不能外推到全部模型",
            "不是 release acceptance",
            "每步成本",
            "迭代数量",
            "host oracle fallback",
        ):
            self.assertIn(phrase, self.source)

    def test_estimator_constants_and_reference_case(self):
        estimator = self.data["estimator"]
        self.assertEqual(estimator, {
            "doubleBytes": 8,
            "columnVectors": 8,
            "rowVectors": 8,
            "matrixCopies": 2,
            "checkPacketBytes": 112,
        })
        rows, cols, nnz, index_bytes = 1000, 500, 10000, 8
        matrix = (
            2 * nnz * (8 + index_bytes)
            + ((rows + 1) + (cols + 1)) * index_bytes
        )
        vectors = 8 * (rows + cols) * 8
        candidate = 2 * (rows + cols) * 8
        self.assertEqual(matrix, 332016)
        self.assertEqual(vectors, 96000)
        self.assertEqual(candidate, 24000)

    def test_performance_workbench_has_scoped_hotspots_and_report_evidence(
            self):
        self.assertTrue(
            {f"hotspot-{number}" for number in range(1, 7)}
            <= self.parser.ids
        )
        for heading in (
            "症状",
            "当前证据",
            "看什么指标",
            "下一步实验",
            "不能破坏",
        ):
            self.assertEqual(self.source.count(f"<h4>{heading}</h4>"), 6)

        report_paths = (
            "docs/hipdlp-phase8-cpu-diagnosis.md",
            "docs/hipdlp-phase8-gpu-performance-report.md",
            "docs/hipdlp-phase8-step-fusion.md",
            "docs/hipdlp-phase8-large20-screen.md",
        )
        for path in report_paths:
            with self.subTest(path=path):
                self.assertIn(f"<code>{path}</code>", self.source)
                self.assertNotIn(path, self.parser.links)

    def test_estimator_interface_validates_and_hides_stale_success(self):
        for element_id in (
            "memory-estimator",
            "rows",
            "cols",
            "nnz",
            "index-bytes",
            "estimator-error",
            "estimator-results",
            "matrix-bytes",
            "vector-bytes",
            "candidate-bytes",
            "packet-bytes",
        ):
            with self.subTest(element_id=element_id):
                self.assertIn(element_id, self.parser.elements)

        self.assertEqual(
            self.parser.elements["estimator-error"]["attributes"]["role"],
            "alert",
        )
        self.assertIn(
            "hidden",
            self.parser.elements["estimator-error"]["attributes"],
        )
        self.assertIn(
            "hidden",
            self.parser.elements["estimator-results"]["attributes"],
        )
        self.assertIn('<option value="4">', self.source)
        self.assertIn('<option value="8" selected>', self.source)

        for behavior in (
            "function parseCount(input, label)",
            "Number.isSafeInteger(value) || value < 0",
            "${label} 必须是非负安全整数",
            "function estimateMemory({ rows, cols, nnz, indexBytes })",
            "function formatBytes(bytes)",
            "function updateEstimator()",
            "估算结果超过 JavaScript 安全整数范围",
        ):
            with self.subTest(behavior=behavior):
                self.assertIn(behavior, self.source)

        update_start = self.source.index("function updateEstimator()")
        update_end = self.source.index(
            'estimatorForm.addEventListener("input"', update_start
        )
        update_source = self.source[update_start:update_end]
        self.assertLess(
            update_source.index("estimatorResults.hidden = true;"),
            update_source.index("try {"),
        )
        self.assertIn("estimatorError.textContent = error.message;", update_source)
        self.assertIn("estimatorError.hidden = false;", update_source)
        self.assertIn("estimatorResults.hidden = false;", update_source)

    def test_estimator_rejects_blank_before_number_conversion(self):
        parse_start = self.source.index("function parseCount(input, label)")
        parse_end = self.source.index(
            "function estimateMemory", parse_start
        )
        parse_source = self.source[parse_start:parse_end]
        blank_guard = 'input.value.trim() === ""'
        self.assertIn(blank_guard, parse_source)
        self.assertLess(
            parse_source.index(blank_guard),
            parse_source.index("Number(input.value)"),
        )

        update_start = self.source.index("function updateEstimator()")
        update_end = self.source.index(
            'estimatorForm.addEventListener("input"', update_start
        )
        update_source = self.source[update_start:update_end]
        self.assertLess(
            update_source.index("estimatorResults.hidden = true;"),
            update_source.index("try {"),
        )

    def test_estimator_restricts_index_width_before_estimating(self):
        update_start = self.source.index("function updateEstimator()")
        update_end = self.source.index(
            'estimatorForm.addEventListener("input"', update_start
        )
        update_source = self.source[update_start:update_end]
        parse_index = "const indexBytes = parseCount("
        width_guard = "if (indexBytes !== 4 && indexBytes !== 8)"
        estimate_call = "const estimates = estimateMemory({"
        self.assertIn(parse_index, update_source)
        self.assertIn(width_guard, update_source)
        self.assertIn(
            "index width 必须是 4 或 8 bytes",
            update_source,
        )
        self.assertLess(
            update_source.index(parse_index),
            update_source.index(width_guard),
        )
        self.assertLess(
            update_source.index(width_guard),
            update_source.index(estimate_call),
        )

    def test_diagnostics_map_connects_solve_fields_to_profilers(self):
        self.assertIn("diagnostics-map", self.parser.ids)
        for field in (
            "backend_setup_seconds",
            "iterations_seconds",
            "step_seconds",
            "boundary_seconds",
            "download_seconds",
            "host_verification_seconds",
            "candidate_downloads",
            "ordinary_check_d2h_bytes",
            "ordinary_check_transfers",
            "ordinary_check_synchronizations",
            "host_oracle_checks",
            "restarts",
            "command_buffer_replays",
        ):
            with self.subTest(field=field):
                self.assertIn(field, self.source)
        for tool in (
            "solver trace",
            "CPU profiler",
            "Nsight Systems",
            "Nsight Compute",
            "Compute Sanitizer",
        ):
            with self.subTest(tool=tool):
                self.assertIn(tool, self.source)

    def test_diagnostics_map_covers_every_solver_field(self):
        solver_source = SOLVER_HEADER.read_text(encoding="utf-8")
        diagnostics = re.search(
            r"struct SolveDiagnostics \{(?P<body>.*?)SolveDiagnostics\(\);",
            solver_source,
            re.DOTALL,
        )
        self.assertIsNotNone(diagnostics)
        solver_fields = set(re.findall(
            r"(?:std::size_t|double)\s+([a-z][a-z0-9_]*);",
            diagnostics.group("body"),
        ))

        table = re.search(
            r'<table id="diagnostics-map".*?</table>',
            self.source,
            re.DOTALL,
        )
        self.assertIsNotNone(table)
        mapped_fields = set(re.findall(
            r"<code>([a-z][a-z0-9_]*)</code>",
            table.group(),
        ))
        self.assertEqual(mapped_fields, solver_fields)

    def test_cuda_replay_card_states_release_acceptance_limit(self):
        card = re.search(
            r'<article class="hotspot-card" id="hotspot-4">'
            r"(?P<body>.*?)</article>",
            self.source,
            re.DOTALL,
        )
        self.assertIsNotNone(card)
        self.assertIn("不能外推到全部模型", card.group("body"))
        self.assertIn("不是 release acceptance", card.group("body"))


if __name__ == "__main__":
    unittest.main()
