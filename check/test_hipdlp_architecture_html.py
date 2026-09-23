import json
from html.parser import HTMLParser
from pathlib import Path
import unittest
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs" / "hipdlp-architecture.html"
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


if __name__ == "__main__":
    unittest.main()
