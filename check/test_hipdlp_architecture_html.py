import json
from html.parser import HTMLParser
from pathlib import Path
import unittest
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs" / "hipdlp-architecture.html"


class ArchitecturePageParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()
        self.tabs = []
        self.remote_assets = []
        self.links = []
        self.data_chunks = []
        self.in_data = False

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        element_id = attributes.get("id")
        if element_id:
            self.ids.add(element_id)
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

    def handle_endtag(self, tag):
        if tag == "script" and self.in_data:
            self.in_data = False

    def handle_data(self, data):
        if self.in_data:
            self.data_chunks.append(data)


class ArchitecturePageTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(PAGE.is_file(), "Create docs/hipdlp-architecture.html")
        self.source = PAGE.read_text(encoding="utf-8")
        self.parser = ArchitecturePageParser()
        self.parser.feed(self.source)
        self.data = json.loads("".join(self.parser.data_chunks))

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
        expected = {
            "node-wrapper",
            "node-api",
            "node-model",
            "node-algorithm",
            "node-backend",
            "node-postsolve",
            "node-kkt",
        }
        self.assertTrue(expected <= self.parser.ids)
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
