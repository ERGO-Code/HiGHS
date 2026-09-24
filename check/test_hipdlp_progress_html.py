import json
from html.parser import HTMLParser
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs" / "hipdlp-rewrite-progress.html"


class ProgressPageParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()
        self.remote_assets = []
        self.tab_count = 0
        self.data_chunks = []
        self.in_data = False

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if attributes.get("id"):
            self.ids.add(attributes["id"])
        if attributes.get("role") == "tab":
            self.tab_count += 1
        for name in ("href", "src"):
            value = attributes.get(name, "")
            if value.startswith(("http://", "https://", "//")):
                self.remote_assets.append(value)
        if tag == "script" and attributes.get("id") == "progress-data":
            self.in_data = True

    def handle_endtag(self, tag):
        if tag == "script" and self.in_data:
            self.in_data = False

    def handle_data(self, data):
        if self.in_data:
            self.data_chunks.append(data)


class ProgressPageTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(PAGE.is_file(), "Create the standalone progress HTML")
        self.source = PAGE.read_text(encoding="utf-8")
        self.parser = ProgressPageParser()
        self.parser.feed(self.source)
        self.data = json.loads("".join(self.parser.data_chunks))

    def test_page_is_offline_and_has_all_views(self):
        self.assertEqual(self.parser.remote_assets, [])
        self.assertGreaterEqual(self.parser.tab_count, 4)
        self.assertTrue({"overview", "pairs", "engineering", "evidence"}
                        <= self.parser.ids)
        self.assertIn("Phase 8 BLOCKED", self.source)
        self.assertIn("为什么不是 COMPLETE", self.source)

    def test_embedded_results_match_strict_phase8_evidence(self):
        self.assertEqual(self.data["commit"], "618c34d589")
        self.assertEqual(self.data["formalJobs"], 92)
        self.assertEqual(self.data["formalRecords"], 552)
        self.assertEqual(self.data["qualification"]["candidateCount"], 30)
        self.assertEqual(self.data["qualification"]["selectedCount"], 23)
        self.assertEqual(self.data["qualification"]["excludedCount"], 7)
        self.assertTrue(self.data["qualification"]["frozenBeforeCandidateExecution"])
        self.assertEqual(self.data["protocol"]["solverTolerances"], {
            "kkt_tolerance": 1e-4,
            "pdlp_optimality_tolerance": 1e-4,
            "primal_feasibility_tolerance": 1e-4,
            "dual_feasibility_tolerance": 1e-4,
        })
        self.assertEqual(self.data["legacyGate"]["cpu"]["plannedPairs"], 23)
        self.assertEqual(self.data["legacyGate"]["cuda"]["plannedPairs"], 23)

    def test_pair_matrix_uses_the_frozen_common_qualified_set(self):
        pairs = {(row["backend"], row["instance"]): row
                 for row in self.data["pairs"]}
        selected = set(self.data["qualification"]["selectedNames"])
        self.assertEqual(len(selected), 23)
        self.assertEqual({instance for backend, instance in pairs}, selected)
        self.assertEqual({backend for backend, instance in pairs}, {"cpu", "cuda"})
        self.assertEqual(len(pairs), 46)
        self.assertNotIn("pypsa_medium", selected)
        self.assertNotIn("pypsa-eur-elec_50_168h", selected)


if __name__ == "__main__":
    unittest.main()