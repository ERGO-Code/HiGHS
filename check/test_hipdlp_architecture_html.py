from collections import Counter
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
VIEW_IDS = {
    "system-map",
    "data-residency",
    "iteration-correctness",
    "bottlenecks",
    "learning-review",
}
COMPONENT_DRAWER_FIELDS = {
    "responsibility",
    "inputs",
    "outputs",
    "owner-lifetime",
    "residency",
    "errors",
    "invariants",
    "source-order",
    "tests",
    "good-signals",
    "danger-signals",
    "review-questions",
}
EVIDENCE_CATEGORIES = {
    "all",
    "fact",
    "measured",
    "hypothesis",
    "unavailable",
}
PRODUCTION_CUSTOM_KERNELS = {
    "fillKernel",
    "scaleKernel",
    "axpyKernel",
    "differenceKernel",
    "projectedPrimalKernel",
    "projectedDualKernel",
    "vectorBlendKernel",
    "vectorPairBlendKernel",
    "averageKernel",
    "finiteCheckKernel",
    "dotKernel",
    "normKernel",
    "differenceNormKernel",
    "checkPartialsKernel",
    "checkPacketKernel",
    "certificateRecoveryKernel",
    "certificateNormKernel",
    "certificateNormalizeKernel",
    "certificateConditionsKernel",
    "certificateFlagsKernel",
}
GENERIC_256_KERNELS = {
    "fillKernel",
    "scaleKernel",
    "axpyKernel",
    "differenceKernel",
    "projectedPrimalKernel",
    "projectedDualKernel",
    "vectorBlendKernel",
    "vectorPairBlendKernel",
    "averageKernel",
    "finiteCheckKernel",
}
GENERIC_KERNEL_GRIDS = {
    "fillKernel": "blockCount(size)",
    "scaleKernel": "blockCount(size)",
    "axpyKernel": "blockCount(size)",
    "differenceKernel": "blockCount(size)",
    "projectedPrimalKernel": "blockCount(num_cols)",
    "projectedDualKernel": "blockCount(num_rows)",
    "vectorBlendKernel": "blockCount(size)",
    "vectorPairBlendKernel": "blockCount(max(primal_size, dual_size))",
    "averageKernel": "blockCount(size)",
    "finiteCheckKernel": "blockCount(size)",
}
SERIAL_1X1_KERNELS = {
    "dotKernel",
    "normKernel",
    "differenceNormKernel",
}
PARTIAL_1024_KERNELS = {
    "checkPartialsKernel",
    "certificateRecoveryKernel",
    "certificateNormalizeKernel",
    "certificateConditionsKernel",
}
FINAL_1X256_KERNELS = {
    "checkPacketKernel",
    "certificateNormKernel",
    "certificateFlagsKernel",
}
CHECK_AND_CERTIFICATE_KERNELS = {
    "checkPartialsKernel",
    "checkPacketKernel",
    "certificateRecoveryKernel",
    "certificateNormKernel",
    "certificateNormalizeKernel",
    "certificateConditionsKernel",
    "certificateFlagsKernel",
}
STATIC_SHARED_BY_KERNEL = {
    "checkPartialsKernel": "CheckPartial[256]",
    "checkPacketKernel": "CheckPartial[256]",
    "certificateRecoveryKernel": "CertificateNormPartial[256]",
    "certificateNormKernel": "CertificateNormPartial[256]",
    "certificateNormalizeKernel": "CertificateNormPartial[256]",
    "certificateConditionsKernel": "ScreenPartial[256]",
    "certificateFlagsKernel": "ScreenPartial[256]",
}
CERTIFICATE_KERNEL_SEQUENCE = [
    "certificateRecoveryKernel",
    "certificateNormKernel(normalization)",
    "certificateNormalizeKernel",
    "certificateNormKernel(verifier)",
    "certificateConditionsKernel",
    "certificateFlagsKernel",
]
CUDA_LESSON_IDS = {
    f"cuda-lesson-{index}" for index in range(1, 9)
}
ANSWER_LABELS = {
    "结论",
    "代码证据",
    "Review 推理",
    "静态边界",
    "运行时工具",
}
INTEGRATED_EXERCISE_TOPICS = {
    "primal-candidate-lifecycle",
    "normal-check-d2h",
    "grid-stride-coverage",
    "graph-dependency-lifetime",
}
REVIEW_PASS_IDS = {
    "interface-ownership",
    "numerical-correctness",
    "concurrency-cuda-semantics",
    "performance",
}
BACKEND_FILTER_CARD_CLASSES = {
    "backend-fact",
    "graph-card",
    "graph-note",
    "hotspot-card",
    "kernel-card",
    "lane-card",
    "resource-card",
    "transfer-item",
}
APPROVED_GLOSSARY_TERMS = {
    "CSR",
    "SpMV",
    "primal",
    "dual",
    "reflection",
    "Halpern anchor",
    "resident",
    "H2D",
    "D2H",
    "pinned memory",
    "stream",
    "event",
    "warp",
    "occupancy",
    "coalescing",
    "divergence",
    "graph capture",
    "deterministic reduction",
    "host oracle",
    "postsolve",
    "KKT",
}
VOID_ELEMENTS = {
    "area",
    "base",
    "br",
    "col",
    "embed",
    "hr",
    "img",
    "input",
    "link",
    "meta",
    "param",
    "source",
    "track",
    "wbr",
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
        self.script_chunks = []
        self.in_script = False
        self.elements = {}
        self.lane_buttons = []
        self.evidence_filter_buttons = []
        self.evidence_cards = []
        self.kernel_inventory_rows = []
        self.backend_lane_cards = []
        self.backend_filter_cards = []
        self.shared_operation_cards = []
        self.style_chunks = []
        self.inline_styles = []
        self.in_style = False
        self.cuda_lessons = 0
        self.total_review_questions = 0
        self.review_questions = 0
        self.kernel_review_questions = 0
        self.answers = 0
        self.lesson_question_counts = {}
        self.lesson_answer_counts = {}
        self.lesson_part_counts = {}
        self.lesson_text = {}
        self.integrated_review_questions = 0
        self.integrated_answers = 0
        self.integrated_exercises = {}
        self.answer_texts = []
        self.answer_attributes = []
        self.review_pass_count = 0
        self.review_passes = {}
        self.review_content_capture = None
        self.glossary_entries = {}
        self.glossary_label_capture = None
        self.element_stack = []
        self.answer_chunks = None
        self.current_node = None
        self.detail_node = None
        self.node_text = {node_id: [] for node_id in SYSTEM_NODE_IDS}
        self.node_detail_count = {node_id: 0 for node_id in SYSTEM_NODE_IDS}
        self.node_detail_text = {node_id: [] for node_id in SYSTEM_NODE_IDS}
        self.node_detail_links = {node_id: [] for node_id in SYSTEM_NODE_IDS}
        self.node_detail_fields = {
            node_id: {} for node_id in SYSTEM_NODE_IDS
        }
        self.detail_field_capture = None

    def has_context(self, element_id=None, class_name=None):
        return any(
            (element_id is None or item["id"] == element_id)
            and (class_name is None or class_name in item["classes"])
            for item in self.element_stack
        )

    def lesson_context(self):
        for item in reversed(self.element_stack):
            if "cuda-lesson" in item["classes"]:
                return item["id"]
        return None

    def attribute_context(self, attribute):
        for item in reversed(self.element_stack):
            value = item["attributes"].get(attribute)
            if value is not None:
                return value
        return None

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        element_id = attributes.get("id")
        if element_id:
            self.ids.add(element_id)
            self.elements[element_id] = {"tag": tag, "attributes": attributes}
        if tag == "article" and element_id in SYSTEM_NODE_IDS:
            self.current_node = element_id
        classes = set(attributes.get("class", "").split())
        in_learning_view = (
            element_id == "learning-review"
            or self.has_context(element_id="learning-review")
        )
        in_kernel_card = (
            "kernel-card" in classes
            or self.has_context(class_name="kernel-card")
        )
        in_integrated_exercises = (
            "integrated-exercises" in classes
            or self.has_context(class_name="integrated-exercises")
        )
        in_cuda_course = (
            "cuda-course-section" in classes
            or self.has_context(class_name="cuda-course-section")
        )
        in_review_mode = (
            "review-mode" in classes
            or self.has_context(class_name="review-mode")
        )
        lesson_id = self.lesson_context()
        if (
            "cuda-lesson" in classes
            and in_learning_view
            and in_cuda_course
        ):
            self.cuda_lessons += 1
            lesson_id = element_id
            self.lesson_question_counts[lesson_id] = 0
            self.lesson_answer_counts[lesson_id] = 0
            self.lesson_text[lesson_id] = []
            self.lesson_part_counts[lesson_id] = {
                "visual": 0,
                "code_example": 0,
                "checklist": 0,
            }
        lesson_parts = {
            "concept-visual": "visual",
            "hipdlp-example": "code_example",
            "lesson-checklist": "checklist",
        }
        if lesson_id in self.lesson_part_counts:
            for class_name, part_name in lesson_parts.items():
                if class_name in classes:
                    self.lesson_part_counts[lesson_id][part_name] += 1
        exercise_topic = self.attribute_context("data-topic")
        if "integrated-exercise" in classes and in_integrated_exercises:
            exercise_topic = attributes.get("data-topic")
            if exercise_topic is not None:
                self.integrated_exercises.setdefault(
                    exercise_topic,
                    {"count": 0, "questions": 0, "answers": 0},
                )["count"] += 1
        if "review-question" in classes:
            self.total_review_questions += 1
            if in_learning_view:
                self.review_questions += 1
                if lesson_id in self.lesson_question_counts:
                    self.lesson_question_counts[lesson_id] += 1
                elif in_integrated_exercises:
                    self.integrated_review_questions += 1
                    if exercise_topic in self.integrated_exercises:
                        self.integrated_exercises[exercise_topic][
                            "questions"
                        ] += 1
            if in_kernel_card:
                self.kernel_review_questions += 1
        if tag == "details" and "answer" in classes and in_learning_view:
            self.answers += 1
            self.answer_chunks = []
            self.answer_attributes.append(attributes)
            if lesson_id in self.lesson_answer_counts:
                self.lesson_answer_counts[lesson_id] += 1
            elif in_integrated_exercises:
                self.integrated_answers += 1
                if exercise_topic in self.integrated_exercises:
                    self.integrated_exercises[exercise_topic]["answers"] += 1
        if (
            "review-pass" in classes
            and in_learning_view
            and in_review_mode
        ):
            self.review_pass_count += 1
            review_pass = attributes.get("data-review-pass")
            if review_pass is not None:
                self.review_passes.setdefault(
                    review_pass,
                    {
                        "count": 0,
                        "good_signals": [],
                        "danger_signals": [],
                        "human_prompts": [],
                        "ai_prompts": [],
                    },
                )["count"] += 1
        review_pass = self.attribute_context("data-review-pass")
        review_content_key = None
        if review_pass in self.review_passes and in_review_mode:
            if "good-signal" in classes:
                review_content_key = "good_signals"
            elif "danger-signal" in classes:
                review_content_key = "danger_signals"
            elif "review-prompt" in classes:
                audience = attributes.get("data-audience")
                if audience in {"human", "ai"}:
                    review_content_key = f"{audience}_prompts"
        if review_content_key is not None:
            self.review_content_capture = {
                "tag": tag,
                "review_pass": review_pass,
                "key": review_content_key,
                "chunks": [],
            }
        in_glossary = (
            "glossary" in classes
            or self.has_context(class_name="glossary")
        )
        if "glossary-entry" in classes and in_glossary:
            term = attributes.get("data-term")
            if term is not None:
                self.glossary_entries.setdefault(
                    term,
                    {"count": 0, "labels": []},
                )["count"] += 1
        glossary_term = self.attribute_context("data-term")
        if tag == "dt" and glossary_term is not None and in_glossary:
            self.glossary_label_capture = {
                "term": glossary_term,
                "chunks": [],
            }
        if tag == "button" and "data-lane" in attributes:
            self.lane_buttons.append(attributes)
        if tag == "button" and "data-evidence-filter" in attributes:
            self.evidence_filter_buttons.append(attributes)
        if "data-evidence" in attributes:
            self.evidence_cards.append({
                "tag": tag,
                "attributes": attributes,
            })
        if "data-kernel" in attributes:
            self.kernel_inventory_rows.append(attributes)
        if "data-backend" in attributes:
            self.backend_filter_cards.append({
                "attributes": attributes,
                "classes": classes,
                "tag": tag,
            })
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
        drawer_field = attributes.get("data-drawer-field")
        if self.detail_node is not None and drawer_field is not None:
            self.detail_field_capture = {
                "tag": tag,
                "node": self.detail_node,
                "field": drawer_field,
                "chunks": [],
            }
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
        srcset = attributes.get("srcset", "")
        for candidate in srcset.split(","):
            value = candidate.strip().split(" ", 1)[0]
            if value.startswith(("http://", "https://", "//")):
                self.remote_assets.append(value)
        inline_style = attributes.get("style")
        if inline_style is not None:
            self.inline_styles.append(inline_style)
        if tag == "script" and element_id == "architecture-data":
            self.in_data = True
        if tag == "script":
            self.in_script = True
        if tag == "style":
            self.in_style = True
        if tag not in VOID_ELEMENTS:
            self.element_stack.append({
                "tag": tag,
                "id": element_id,
                "classes": classes,
                "attributes": attributes,
            })

    def handle_endtag(self, tag):
        if (
            self.detail_field_capture is not None
            and tag == self.detail_field_capture["tag"]
        ):
            content = " ".join(
                "".join(self.detail_field_capture["chunks"]).split()
            )
            node = self.detail_field_capture["node"]
            field = self.detail_field_capture["field"]
            self.node_detail_fields[node].setdefault(field, []).append(content)
            self.detail_field_capture = None
        if (
            self.review_content_capture is not None
            and tag == self.review_content_capture["tag"]
        ):
            content = " ".join(
                "".join(
                    self.review_content_capture["chunks"]
                ).split()
            )
            review_pass = self.review_content_capture["review_pass"]
            key = self.review_content_capture["key"]
            self.review_passes[review_pass][key].append(content)
            self.review_content_capture = None
        if tag == "dt" and self.glossary_label_capture is not None:
            label = " ".join(
                "".join(self.glossary_label_capture["chunks"]).split()
            )
            term = self.glossary_label_capture["term"]
            self.glossary_entries[term]["labels"].append(label)
            self.glossary_label_capture = None
        if tag == "details" and self.answer_chunks is not None:
            self.answer_texts.append(
                " ".join("".join(self.answer_chunks).split())
            )
            self.answer_chunks = None
        if tag == "script" and self.in_data:
            self.in_data = False
        if tag == "script":
            self.in_script = False
        if tag == "style":
            self.in_style = False
        if tag == "details" and self.detail_node is not None:
            self.detail_node = None
        if tag == "article" and self.current_node is not None:
            self.current_node = None
        if self.element_stack:
            self.element_stack.pop()

    def handle_data(self, data):
        if self.in_data:
            self.data_chunks.append(data)
        if self.in_script:
            self.script_chunks.append(data)
        if self.in_style:
            self.style_chunks.append(data)
        lesson_id = self.lesson_context()
        if lesson_id in self.lesson_text:
            self.lesson_text[lesson_id].append(data)
        if self.current_node is not None:
            self.node_text[self.current_node].append(data)
        if self.detail_node is not None:
            self.node_detail_text[self.detail_node].append(data)
        if self.detail_field_capture is not None:
            self.detail_field_capture["chunks"].append(data)
        if self.answer_chunks is not None:
            self.answer_chunks.append(data)
        if self.review_content_capture is not None:
            self.review_content_capture["chunks"].append(data)
        if self.glossary_label_capture is not None:
            self.glossary_label_capture["chunks"].append(data)


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
        self.assertTrue(VIEW_IDS <= self.parser.ids)
        css = "\n".join(
            self.parser.style_chunks + self.parser.inline_styles
        )
        css_network_patterns = (
            r"@import\s+(?:url\(\s*)?[\"']?\s*(?:https?:)?//",
            r"url\(\s*[\"']?\s*(?:https?:)?//",
        )
        for pattern in css_network_patterns:
            with self.subTest(css_network_pattern=pattern):
                self.assertIsNone(re.search(pattern, css, re.IGNORECASE))

        javascript = "\n".join(self.parser.script_chunks)
        javascript_network_patterns = (
            r"\bfetch\s*\(\s*[\"'](?:https?:)?//",
            r"\.open\s*\(\s*[\"'][A-Z]+[\"']\s*,\s*"
            r"[\"'](?:https?:)?//",
            r"\b(?:new\s+)?WebSocket\s*\(\s*[\"'](?:wss?:)?//",
            r"\b(?:new\s+)?EventSource\s*\(\s*[\"'](?:https?:)?//",
            r"\bimport\s*\(\s*[\"'](?:https?:)?//",
            r"\bimport\s+(?:[^;\n]*?\s+from\s+)?"
            r"[\"'](?:https?:)?//",
        )
        for pattern in javascript_network_patterns:
            with self.subTest(javascript_network_pattern=pattern):
                self.assertIsNone(
                    re.search(pattern, javascript, re.IGNORECASE)
                )
        self.assertEqual(self.data["snapshot"], {
            "branch": "feature/hipdlp-rewrite",
            "commit": "49b85f4f3",
            "generated": "2026-09-23",
            "workingTree": "dirty",
        })

    def test_tabs_reference_existing_views_and_static_content_is_not_hidden(
            self):
        controls = {tab["aria-controls"] for tab in self.parser.tabs}
        self.assertEqual(controls, VIEW_IDS)
        self.assertTrue(controls <= self.parser.ids)
        for view_id in sorted(controls):
            with self.subTest(view_id=view_id):
                panel = self.parser.elements[view_id]
                self.assertEqual(panel["tag"], "section")
                self.assertEqual(
                    panel["attributes"].get("role"),
                    "tabpanel",
                )
                self.assertNotIn("hidden", panel["attributes"])
        self.assertIn("prefers-reduced-motion: reduce", self.source)
        self.assertIn("transition: none !important;", self.source)
        self.assertIn("@media print", self.source)
        self.assertIn('aria-live="polite"', self.source)

    def test_interaction_functions_are_present_without_external_frameworks(
            self):
        for token in (
            "function activateView(",
            "function activateLane(",
            "function estimateMemory(",
            "function updateEstimator(",
            "history.replaceState",
            "ArrowRight",
            "ArrowLeft",
            "Home",
            "End",
        ):
            with self.subTest(token=token):
                self.assertIn(token, self.source)
        for forbidden in ("React", "Vue", "d3.", "mermaid", "cdn."):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, self.source)

    def test_every_local_page_link_resolves(self):
        missing = []
        for link in self.parser.links:
            parts = urlsplit(link)
            if parts.scheme or link.startswith("#"):
                continue
            target = unquote(parts.path)
            if target and not (PAGE.parent / target).resolve().exists():
                missing.append(link)
        self.assertEqual(missing, [])

    def test_learning_view_has_complete_routes_lessons_and_answers(self):
        self.assertTrue({"route-30", "route-cpp-90", "route-cuda"}
                        <= self.parser.ids)
        self.assertEqual(self.parser.cuda_lessons, 8)
        self.assertEqual(
            set(self.parser.lesson_part_counts),
            CUDA_LESSON_IDS,
        )
        self.assertEqual(
            self.parser.lesson_part_counts,
            {
                lesson_id: {
                    "visual": 1,
                    "code_example": 1,
                    "checklist": 1,
                }
                for lesson_id in CUDA_LESSON_IDS
            },
        )
        self.assertEqual(self.parser.review_questions, 28)
        self.assertEqual(self.parser.answers, 28)
        self.assertEqual(
            self.parser.lesson_question_counts,
            {lesson_id: 3 for lesson_id in CUDA_LESSON_IDS},
        )
        self.assertEqual(
            self.parser.lesson_answer_counts,
            {lesson_id: 3 for lesson_id in CUDA_LESSON_IDS},
        )
        self.assertEqual(self.parser.integrated_review_questions, 4)
        self.assertEqual(self.parser.integrated_answers, 4)
        for heading in (
            "执行模型",
            "内存模型",
            "异步与 Stream",
            "性能基础",
            "cuSPARSE",
            "CUDA Graph",
            "数值与确定性",
            "调试工具地图",
        ):
            self.assertIn(heading, self.source)

    def test_lessons_two_and_three_cover_memory_and_timing_in_context(self):
        memory_lesson = self.normalized(
            self.parser.lesson_text["cuda-lesson-2"]
        )
        for topic in (
            "Host memory",
            "pinned Host memory",
            "global Device memory",
            "register",
            "static shared memory",
            "CheckPartial partial",
            "__shared__ CheckPartial shared[256]",
            "dynamic shared memory = 0",
        ):
            with self.subTest(lesson=2, topic=topic):
                self.assertIn(topic, memory_lesson)

        timing_lesson = self.normalized(
            self.parser.lesson_text["cuda-lesson-3"]
        )
        for topic in (
            "cudaEventRecord",
            "cudaEventSynchronize",
            "cudaEventElapsedTime",
            "cudaStreamSynchronize",
            "GPU elapsed time",
            "Host synchronization wall time",
        ):
            with self.subTest(lesson=3, topic=topic):
                self.assertIn(topic, timing_lesson)

    def test_kernel_prompts_remain_review_questions(self):
        self.assertEqual(self.parser.total_review_questions, 34)
        self.assertEqual(self.parser.kernel_review_questions, 6)
        self.assertNotIn("kernel-review-question", self.source)

    def test_learning_answers_separate_static_and_runtime_evidence(self):
        self.assertEqual(len(self.parser.answer_texts), 28)
        for index, answer_text in enumerate(self.parser.answer_texts, start=1):
            with self.subTest(answer=index):
                for label in ANSWER_LABELS:
                    self.assertIn(f"{label}：", answer_text)
        for index, attributes in enumerate(
            self.parser.answer_attributes,
            start=1,
        ):
            with self.subTest(answer=index):
                self.assertNotIn("open", attributes)

    def test_integrated_exercises_cover_each_required_topic_once(self):
        self.assertEqual(
            set(self.parser.integrated_exercises),
            INTEGRATED_EXERCISE_TOPICS,
        )
        self.assertEqual(
            self.parser.integrated_exercises,
            {
                topic: {"count": 1, "questions": 1, "answers": 1}
                for topic in INTEGRATED_EXERCISE_TOPICS
            },
        )

    def test_review_mode_has_four_complete_review_passes(self):
        self.assertEqual(self.parser.review_pass_count, 4)
        self.assertEqual(set(self.parser.review_passes), REVIEW_PASS_IDS)
        for review_pass, content in self.parser.review_passes.items():
            with self.subTest(review_pass=review_pass):
                self.assertEqual(content["count"], 1)
                for content_type in (
                    "good_signals",
                    "danger_signals",
                    "human_prompts",
                    "ai_prompts",
                ):
                    self.assertEqual(len(content[content_type]), 1)
                    self.assertTrue(content[content_type][0].strip())

    def test_glossary_marks_every_approved_term_contextually(self):
        self.assertEqual(
            set(self.parser.glossary_entries),
            APPROVED_GLOSSARY_TERMS,
        )
        for term, entry in self.parser.glossary_entries.items():
            with self.subTest(term=term):
                self.assertEqual(entry["count"], 1)
                self.assertEqual(entry["labels"], [term])

    def test_candidate_download_paths_distinguish_certificate_verification(self):
        self.assertNotIn(
            "threshold uncertainty / optional metric fallback / "
            "suspected certificate",
            self.source,
        )
        for evidence in (
            "普通 cadence 也可能触发完整 candidate D2H",
            "<code>requires_host_oracle</code>",
            "<code>certificate_possible</code>",
            "<code>detectCertificate()</code>",
            "raw convergence host oracle",
            "independent certificate verification",
        ):
            self.assertIn(evidence, self.source)

    def test_grid_stride_proof_is_scoped_to_grid_index_kernels(self):
        self.assertIn(
            "该证明只覆盖显式使用 <code>GridIndex</code> 的",
            self.source,
        )
        self.assertIn(
            "<code>differenceKernel</code> 和 "
            "<code>vectorPairBlendKernel</code>",
            self.source,
        )
        self.assertIn(
            "<code>projectedPrimalKernel</code> 等仍用 "
            "<code>HighsInt</code>",
            self.source,
        )
        self.assertNotIn(
            "宽索引 <code>GridIndex</code> 与 static_assert "
            "防止支持范围内的 stride 加法溢出",
            self.source,
        )

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
        for node_id in sorted(SYSTEM_NODE_IDS):
            with self.subTest(node_id=node_id):
                self.assertEqual(self.parser.node_detail_count[node_id], 1)
                fields = self.parser.node_detail_fields[node_id]
                self.assertEqual(set(fields), COMPONENT_DRAWER_FIELDS)
                for field_name, values in fields.items():
                    with self.subTest(
                        node_id=node_id,
                        field_name=field_name,
                    ):
                        self.assertEqual(len(values), 1)
                        self.assertTrue(values[0].strip())
                review_question = fields["review-questions"][0]
                self.assertGreaterEqual(len(review_question), 12)
                self.assertIn("？", review_question)
                source_order = fields["source-order"][0]
                self.assertIn("→", source_order)
                for signal_field in ("good-signals", "danger-signals"):
                    self.assertGreaterEqual(
                        len(fields[signal_field][0]),
                        12,
                    )
                links = self.parser.node_detail_links[node_id]
                self.assertEqual(len(links), 1)
                self.assertTrue(links[0].startswith("../highs/"))

    def test_evidence_filter_has_accessible_controls_and_immutable_cards(self):
        self.assertEqual(
            {
                button["data-evidence-filter"]
                for button in self.parser.evidence_filter_buttons
            },
            EVIDENCE_CATEGORIES,
        )
        self.assertEqual(len(self.parser.evidence_filter_buttons), 5)
        for button in self.parser.evidence_filter_buttons:
            category = button["data-evidence-filter"]
            with self.subTest(category=category):
                self.assertEqual(button.get("type"), "button")
                self.assertEqual(
                    button.get("aria-pressed"),
                    "true" if category == "all" else "false",
                )

        card_categories = {
            card["attributes"]["data-evidence"]
            for card in self.parser.evidence_cards
        }
        self.assertEqual(
            card_categories,
            EVIDENCE_CATEGORIES - {"all"},
        )
        self.assertTrue(self.parser.evidence_cards)
        self.assertTrue(all(
            "hidden" not in card["attributes"]
            for card in self.parser.evidence_cards
        ))
        status = self.parser.elements["evidence-filter-status"]["attributes"]
        self.assertEqual(status.get("role"), "status")
        self.assertEqual(status.get("aria-live"), "polite")
        self.assertEqual(status.get("aria-atomic"), "true")

        function_start = self.source.index(
            "function activateEvidenceFilter("
        )
        function_end = self.source.index(
            "function activateLane(", function_start
        )
        function_source = self.source[function_start:function_end]
        for behavior in (
            "card.hidden =",
            'button.setAttribute("aria-pressed",',
            "evidenceFilterStatus.textContent =",
        ):
            with self.subTest(behavior=behavior):
                self.assertIn(behavior, function_source)
        self.assertNotIn("card.textContent", function_source)
        self.assertNotIn("card.innerHTML", function_source)

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
        for name in SERIAL_1X1_KERNELS:
            with self.subTest(kernel=name):
                self.assertIn(name, kernels)
                self.assertEqual(kernels[name]["threads"], 1)

    def test_cuda_kernel_inventory_is_complete_and_source_accurate(self):
        kernels = {
            item["name"]: item
            for item in self.data["kernels"]
            if item.get("kind") == "custom"
        }
        self.assertEqual(set(kernels), PRODUCTION_CUSTOM_KERNELS)
        self.assertEqual(
            {
                row["data-kernel"]
                for row in self.parser.kernel_inventory_rows
            },
            PRODUCTION_CUSTOM_KERNELS | {"cusparseSpMV"},
        )
        for name, kernel in kernels.items():
            with self.subTest(kernel=name):
                self.assertEqual(kernel["dynamicSharedBytes"], 0)
                self.assertTrue(kernel["logicalElements"])
                self.assertTrue(kernel["buffers"])
                self.assertTrue(kernel["reviewRisk"])
                self.assertTrue(kernel["grid"])
                self.assertTrue(kernel["staticShared"])

        for name in GENERIC_256_KERNELS:
            with self.subTest(kernel=name, geometry="generic"):
                self.assertEqual(
                    kernels[name]["grid"],
                    GENERIC_KERNEL_GRIDS[name],
                )
                self.assertEqual(kernels[name]["threads"], 256)
                self.assertEqual(kernels[name]["maxBlocks"], 65535)
                self.assertEqual(kernels[name]["staticShared"], "none")
        for name in SERIAL_1X1_KERNELS:
            with self.subTest(kernel=name, geometry="serial"):
                self.assertEqual(kernels[name]["grid"], "1")
                self.assertEqual(kernels[name]["blocks"], 1)
                self.assertEqual(kernels[name]["threads"], 1)
                self.assertEqual(kernels[name]["staticShared"], "none")
        for name in CHECK_AND_CERTIFICATE_KERNELS:
            with self.subTest(kernel=name, geometry="reduction"):
                self.assertEqual(kernels[name]["threads"], 256)
                self.assertEqual(
                    kernels[name]["staticShared"],
                    STATIC_SHARED_BY_KERNEL[name],
                )
        for name in PARTIAL_1024_KERNELS:
            with self.subTest(kernel=name, geometry="partial"):
                self.assertEqual(
                    kernels[name]["grid"],
                    "checkPartialCount(rows, cols)",
                )
                self.assertEqual(kernels[name]["maxBlocks"], 1024)
        for name in FINAL_1X256_KERNELS:
            with self.subTest(kernel=name, geometry="final"):
                self.assertEqual(kernels[name]["grid"], "1")
                self.assertEqual(kernels[name]["blocks"], 1)

        vendor = {
            item["name"]: item
            for item in self.data["kernels"]
            if item.get("kind") == "vendor"
        }
        self.assertEqual(set(vendor), {"cusparseSpMV"})
        self.assertEqual(vendor["cusparseSpMV"]["grid"], "vendor-controlled")

        self.assertEqual(
            self.data["certificateKernelSequence"],
            CERTIFICATE_KERNEL_SEQUENCE,
        )
        self.assertEqual(
            self.data["testOnlyKernels"],
            [{
                "name": "injectedFailureKernel",
                "guard": "HIGHS_HIPDLP_TESTING",
                "launch": "<<<1, 1, 0, stream>>>",
            }],
        )
        self.assertNotIn(
            "所有 custom launches 使用 256 threads",
            self.source,
        )
        self.assertIn(
            "cuSPARSE launch geometry 由 vendor 控制",
            self.source,
        )

    def test_certificate_screen_has_explicit_eligibility_and_fallback(self):
        self.assertIn("certificateScreen", self.data)
        screen = self.data["certificateScreen"]
        self.assertEqual(screen["architecture"], "AArch64")
        self.assertEqual(screen["glibcCompile"], "2.39")
        self.assertEqual(screen["glibcRuntime"], "2.39")
        self.assertFalse(screen["fastMath"])
        self.assertEqual(screen["fltRadix"], 2)
        self.assertEqual(screen["dblMantDig"], 53)
        self.assertEqual(screen["dblMaxExp"], 1024)
        self.assertEqual(screen["fltEvalMethod"], 0)
        self.assertEqual(screen["rounding"], "FE_TONEAREST")
        self.assertEqual(screen["fpcrClearMask"], "((3 << 22) | (1 << 24) | 3)")
        self.assertEqual(
            screen["dimensionLimitExclusive"],
            4503599627370496,
        )
        self.assertEqual(
            screen["ineligibleDefaults"],
            {
                "certificateScreenValid": False,
                "primalRayPossible": True,
                "farkasPossible": True,
            },
        )
        self.assertIn("certificate-screen-gate", self.parser.elements)
        for limitation in (
            "__aarch64__",
            "__GLIBC__ == 2",
            "__GLIBC_MINOR__ == 39",
            "!defined(__FAST_MATH__)",
            "FLT_RADIX == 2",
            "DBL_MANT_DIG == 53",
            "DBL_MAX_EXP == 1024",
            "LDBL_MANT_DIG &gt;= 53",
            "LDBL_MAX_EXP &gt;= DBL_MAX_EXP",
            "FLT_EVAL_METHOD == 0",
            'gnu_get_libc_version() == "2.39"',
            "fegetround() == FE_TONEAREST",
            "num_rows + num_cols < 4503599627370496",
            "certificate_screen_valid = false",
            "primal_ray_possible = true",
            "farkas_possible = true",
            "certificate_possible = true",
            "requires_host_oracle = true",
            "full-candidate Host certificate verification",
        ):
            with self.subTest(limitation=limitation):
                self.assertIn(limitation, self.source)

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

    def test_lane_controls_filter_backend_cards_and_default_to_compare(self):
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
            'document.querySelectorAll("[data-backend]")',
            'button.setAttribute("aria-pressed", String(isActive))',
            'card.hidden = lane !== "compare" && card.dataset.backend !== lane',
            'button.addEventListener("click", () => activateLane(',
            'activateLane("compare")',
        ):
            with self.subTest(behavior=behavior):
                self.assertIn(behavior, self.source)

    def test_lane_filter_covers_relevant_views_and_announces_results(self):
        counts = Counter()
        for card in self.parser.backend_filter_cards:
            kinds = BACKEND_FILTER_CARD_CLASSES & card["classes"]
            self.assertEqual(len(kinds), 1, card)
            counts[(kinds.pop(), card["attributes"]["data-backend"])] += 1
            self.assertNotIn("hidden", card["attributes"])

        self.assertEqual(counts, Counter({
            ("backend-fact", "cuda"): 4,
            ("graph-card", "cuda"): 1,
            ("graph-note", "cuda"): 1,
            ("hotspot-card", "cpu"): 1,
            ("hotspot-card", "cuda"): 3,
            ("kernel-card", "cuda"): 6,
            ("lane-card", "cpu"): 5,
            ("lane-card", "cuda"): 5,
            ("resource-card", "cpu"): 1,
            ("resource-card", "cuda"): 1,
            ("transfer-item", "cuda"): 4,
        }))

        status = self.parser.elements["lane-filter-status"]["attributes"]
        self.assertEqual(status["role"], "status")
        self.assertEqual(status["aria-live"], "polite")
        self.assertEqual(status["aria-atomic"], "true")
        for behavior in (
            'const backendLaneCards = document.querySelectorAll(',
            '"[data-backend]"',
            'const laneFilterStatus = document.querySelector(',
            '"#lane-filter-status"',
            "laneFilterStatus.textContent =",
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
        for exact_output in (
            "332016 bytes",
            "96000 bytes",
            "24000 bytes",
            "112 bytes",
        ):
            with self.subTest(exact_output=exact_output):
                self.assertIn(exact_output, self.source)

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
            "`${bytes} bytes`",
            "function updateEstimator(root = document)",
            "估算结果超过 JavaScript 安全整数范围",
        ):
            with self.subTest(behavior=behavior):
                self.assertIn(behavior, self.source)

        update_start = self.source.index(
            "function updateEstimator(root = document)"
        )
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

        update_start = self.source.index(
            "function updateEstimator(root = document)"
        )
        update_end = self.source.index(
            'estimatorForm.addEventListener("input"', update_start
        )
        update_source = self.source[update_start:update_end]
        self.assertLess(
            update_source.index("estimatorResults.hidden = true;"),
            update_source.index("try {"),
        )

    def test_estimator_restricts_index_width_before_estimating(self):
        update_start = self.source.index(
            "function updateEstimator(root = document)"
        )
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

    def test_estimator_runs_detached_executable_contract_cases(self):
        update_signature = "function updateEstimator(root = document)"
        contract_signature = "function runEstimatorContractTests()"
        self.assertIn(update_signature, self.source)
        self.assertIn(contract_signature, self.source)

        update_start = self.source.index(update_signature)
        update_end = self.source.index(contract_signature, update_start)
        update_source = self.source[update_start:update_end]
        for element_id in (
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
            with self.subTest(root_selector=element_id):
                self.assertIn(
                    f'root.querySelector("#{element_id}")',
                    update_source,
                )

        contract_start = self.source.index(contract_signature)
        contract_end = self.source.index(
            'const estimatorForm = document.querySelector("#memory-estimator")',
            contract_start,
        )
        contract_source = self.source[contract_start:contract_end]
        for behavior in (
            'const fixture = document.createElement("div");',
            "fixture.innerHTML = `",
            "for (const testCase of cases)",
            "updateEstimator(fixture);",
            'name: "blank rows"',
            'rows: ""',
            'name: "4-byte indices"',
            'indexBytes: "4"',
            'matrix: "246008 bytes (240.24 KiB)"',
            'name: "8-byte indices"',
            'indexBytes: "8"',
            'matrix: "332016 bytes (324.23 KiB)"',
            'vector: "96000 bytes (93.75 KiB)"',
            'candidate: "24000 bytes (23.44 KiB)"',
            'packet: "112 bytes"',
            'name: "16-byte indices"',
            'indexBytes: "16"',
            "resultsHidden: true",
            "resultsHidden: false",
            "errorHidden: true",
            "errorHidden: false",
            'errorText: "rows 不能为空"',
            'errorText: "index width 必须是 4 或 8 bytes"',
            "testCase.expected.resultsHidden",
            "testCase.expected.errorHidden",
            'document.documentElement.dataset.estimatorSelftest = "pass";',
            'document.documentElement.dataset.estimatorSelftest = "fail";',
            'console.error("Estimator self-test failed", error);',
            "throw error;",
        ):
            with self.subTest(behavior=behavior):
                self.assertIn(behavior, contract_source)

        self.assertIn(
            'estimatorForm.addEventListener("input", () => updateEstimator());',
            self.source,
        )
        self.assertRegex(
            self.source,
            r"runEstimatorContractTests\(\);\s+updateEstimator\(\);",
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
            r'<article class="hotspot-card" id="hotspot-4"[^>]*>'
            r"(?P<body>.*?)</article>",
            self.source,
            re.DOTALL,
        )
        self.assertIsNotNone(card)
        self.assertIn("不能外推到全部模型", card.group("body"))
        self.assertIn("不是 release acceptance", card.group("body"))


if __name__ == "__main__":
    unittest.main()
