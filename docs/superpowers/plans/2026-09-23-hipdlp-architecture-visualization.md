# HiPDLP Architecture Explorer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a dependency-free Chinese HTML explorer that visualizes the current HiPDLP rewrite architecture, CPU/GPU data flow, solver-owned data structures, measured performance bottlenecks, and a guided C++/CUDA review course.

**Architecture:** One self-contained HTML file owns semantic content, embedded CSS, inline SVG diagrams, an embedded JSON fact model, and progressive-enhancement JavaScript. One Python standard-library test parses the file, validates its offline/static contracts and source facts, and supplies focused red-green coverage while the five views are added incrementally.

**Tech Stack:** HTML5, CSS, vanilla JavaScript, inline SVG, Python 3.12 standard library (`unittest`, `html.parser`, `json`, `pathlib`), VS Code integrated browser.

**Spec:** `docs/superpowers/specs/2026-09-23-hipdlp-architecture-visualization-design.md`

## Global Constraints

- Create exactly one deliverable page: `docs/hipdlp-architecture.html`.
- Create exactly one page contract test: `check/test_hipdlp_architecture_html.py`.
- Do not modify `docs/hipdlp-rewrite-progress.html` or `check/test_hipdlp_progress_html.py`.
- Do not modify solver source, build files, dependencies, or solver behavior.
- Keep the page fully usable through `file://`; make no remote request and use no CDN.
- Keep all core explanations in the HTML DOM. JavaScript may hide or filter content only after it loads.
- Label source facts, dated measurements, optimization hypotheses, and unavailable future backends distinctly.
- Identify the architecture snapshot as branch `feature/hipdlp-rewrite`, commit `49b85f4f3`, generated 2026-09-23, with a dirty-working-tree evidence caveat.
- Preserve the approved visual legend: host blue, device green, H2D orange, D2H red, control gray, correctness purple, measured hotspot amber.
- Preserve exact current facts: eight column vectors, eight row vectors, two SpMVs per production step, a 39-step replay, 117 doubles/936 bytes of replay parameters, a 112-byte check packet, 256 threads per custom-kernel block, 65,535 general blocks, and 1,024 check blocks.
- Treat cuSPARSE launch geometry as vendor-controlled rather than inventing a grid configuration.
- Do not describe backend convergence as sufficient for `OPTIMAL`; final original-model postsolve/KKT remains a separate CPU gate.
- Show explicit error exits for option/model validation, backend availability,
  allocation, dimension/alias checks, non-finite arithmetic, CUDA/library
  failures, graph failures, resident-check failures, postsolve, and KKT.
- Run Python checks with `/home/yanyu/github/hipdlp2/.venv/bin/python`.
- Preserve unrelated dirty work. Stage only the two files named by the current task.

## File Structure

| File | Responsibility |
| --- | --- |
| `docs/hipdlp-architecture.html` | All visual content, fact data, SVGs, styling, interactions, estimator, lessons, answers, accessibility behavior, and print behavior |
| `check/test_hipdlp_architecture_html.py` | Offline-asset checks, five-view contract, embedded-data assertions, source-link validation, content counts, fact checks, and accessibility/static-markup checks |

The single-file page is intentional. Splitting CSS or JavaScript would violate
the approved portable offline artifact. Internal page sections, JSON objects,
CSS component classes, and small JavaScript functions provide the boundaries
that separate responsibilities.

---

### Task 1: Establish The Offline Page Contract And Five-View Shell

**Files:**
- Create: `check/test_hipdlp_architecture_html.py`
- Create: `docs/hipdlp-architecture.html`

**Interfaces:**
- Consumes: no earlier task output.
- Produces: page IDs `system-map`, `data-residency`, `iteration-correctness`, `bottlenecks`, and `learning-review`; one `<script id="architecture-data" type="application/json">`; the reusable Python `ArchitecturePageParser`.

- [ ] **Step 1: Write the failing page-contract test**

Create `check/test_hipdlp_architecture_html.py` with this initial parser and
test:

```python
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


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run the test and observe the missing-page failure**

Run:

```bash
/home/yanyu/github/hipdlp2/.venv/bin/python \
  check/test_hipdlp_architecture_html.py
```

Expected: FAIL at `assertTrue(PAGE.is_file())` because the page does not exist.

- [ ] **Step 3: Create the semantic five-view shell**

Create `docs/hipdlp-architecture.html` with:

```html
<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="color-scheme" content="light">
  <link rel="icon" href="data:,">
  <title>HiPDLP Rewrite 架构与数据流</title>
  <style>
    :root {
      --paper: #f4f6f8;
      --surface: #ffffff;
      --ink: #17212b;
      --muted: #5d6875;
      --line: #cbd3dc;
      --host: #2563a6;
      --host-soft: #e5f0fb;
      --device: #138a63;
      --device-soft: #e3f5ee;
      --h2d: #c66a13;
      --d2h: #c43d4b;
      --control: #667085;
      --correctness: #7251b5;
      --hotspot: #b7791f;
      --radius: 8px;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      color: var(--ink);
      background: var(--paper);
      font: 15px/1.6 "Noto Sans CJK SC", "Microsoft YaHei", sans-serif;
    }
    .shell { width: min(1440px, calc(100% - 32px)); margin: 0 auto; }
    .tablist { display: flex; overflow-x: auto; gap: 4px; }
    [role="tab"] { border: 0; border-bottom: 3px solid transparent; padding: 14px 16px; }
    [role="tab"][aria-selected="true"] { border-bottom-color: var(--host); }
    .view { margin: 24px 0; padding: 24px; background: var(--surface); border: 1px solid var(--line); border-radius: var(--radius); }
  </style>
</head>
<body>
  <header class="masthead">
    <div class="shell">
      <p>CODE ATLAS · 2026-09-23</p>
      <h1>HiPDLP Rewrite 架构与数据流</h1>
      <p>源码快照 feature/hipdlp-rewrite @ 49b85f4f3；工作区包含未提交证据文档。</p>
    </div>
  </header>
  <nav class="shell tablist" role="tablist" aria-label="架构视图">
    <button role="tab" aria-selected="true" aria-controls="system-map" data-view="system-map">系统地图</button>
    <button role="tab" aria-selected="false" aria-controls="data-residency" data-view="data-residency">数据与驻留</button>
    <button role="tab" aria-selected="false" aria-controls="iteration-correctness" data-view="iteration-correctness">迭代与正确性</button>
    <button role="tab" aria-selected="false" aria-controls="bottlenecks" data-view="bottlenecks">卡点与优化</button>
    <button role="tab" aria-selected="false" aria-controls="learning-review" data-view="learning-review">学习与 Review</button>
  </nav>
  <main class="shell">
    <section class="view" id="system-map"><h2>系统地图</h2></section>
    <section class="view" id="data-residency"><h2>数据与驻留</h2></section>
    <section class="view" id="iteration-correctness"><h2>迭代与正确性</h2></section>
    <section class="view" id="bottlenecks"><h2>卡点与优化</h2></section>
    <section class="view" id="learning-review"><h2>学习与 Review</h2></section>
  </main>
  <script id="architecture-data" type="application/json">
  {
    "snapshot": {
      "branch": "feature/hipdlp-rewrite",
      "commit": "49b85f4f3",
      "generated": "2026-09-23",
      "workingTree": "dirty"
    }
  }
  </script>
</body>
</html>
```

Do not add `hidden` to the four inactive sections yet. With JavaScript absent,
all five sections must be readable.

- [ ] **Step 4: Run the contract test and verify green**

Run:

```bash
/home/yanyu/github/hipdlp2/.venv/bin/python \
  check/test_hipdlp_architecture_html.py
```

Expected: 1 test passes.

- [ ] **Step 5: Commit the shell and contract**

```bash
git add docs/hipdlp-architecture.html \
  check/test_hipdlp_architecture_html.py
git commit -m "Add HiPDLP architecture explorer shell" \
  -m "Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Add The End-To-End System Map

**Files:**
- Modify: `check/test_hipdlp_architecture_html.py`
- Modify: `docs/hipdlp-architecture.html`

**Interfaces:**
- Consumes: the five-view shell and embedded `architecture-data`.
- Produces: system node IDs `node-wrapper`, `node-api`, `node-model`, `node-algorithm`, `node-backend`, `node-postsolve`, and `node-kkt`; source links rooted at `../highs/`.

- [ ] **Step 1: Add failing system-map and source-link tests**

Add:

```python
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
```

- [ ] **Step 2: Run the focused tests and observe missing node IDs**

Run:

```bash
/home/yanyu/github/hipdlp2/.venv/bin/python \
  -m unittest \
  check.test_hipdlp_architecture_html.ArchitecturePageTests.test_system_map_has_the_real_solve_boundaries \
  check.test_hipdlp_architecture_html.ArchitecturePageTests.test_all_local_source_links_resolve
```

Expected: both tests fail because the map and links are absent.

- [ ] **Step 3: Add the control-flow SVG and explanatory cards**

In `#system-map`, add:

- a CPU control-plane lane;
- a backend fork with CPU and CUDA branches;
- a return path through postsolve and KKT;
- a visible legacy/rewrite decision at the API node;
- explicit `auto` behavior: try CUDA creation, otherwise select CPU;
- explicit requested-CUDA behavior: return an error, do not silently use CPU.

Use source cards following this exact shape:

```html
<article class="node host" id="node-api">
  <span class="fact-tag fact">已实现</span>
  <h3>API 路由与选项翻译</h3>
  <code>hipdlp::api::solve</code>
  <p>翻译 HighsOptions，选择 legacy/rewrite，并把剩余全局时间预算交给 rewrite。</p>
  <dl>
    <div><dt>输入</dt><dd>HighsOptions、HighsLp、MigrationOptions</dd></div>
    <div><dt>输出</dt><dd>SolveResult → HighsSolution/HighsInfo/status</dd></div>
    <div><dt>执行位置</dt><dd>CPU / Host</dd></div>
  </dl>
  <a href="../highs/pdlp/hipdlp/api/solver.cc">打开源码</a>
</article>
```

Create analogous cards for all seven IDs. Use inline SVG arrows with text
labels rather than color alone.

- [ ] **Step 4: Make every system node independently reviewable**

Inside each node, place ownership, input/output, residency, error paths, source,
and review questions in a native detail element:

```html
<details class="component-detail">
  <summary>展开接口、所有权与 Review 要点</summary>
  <dl>
    <div><dt>所有权</dt><dd>API 借用 HighsLp；SolveResult 拥有返回解。</dd></div>
    <div><dt>失败</dt><dd>选项翻译错误直接形成结构化失败结果。</dd></div>
  </dl>
  <p><strong>Review：</strong>剩余时间预算是否只扣除一次？</p>
</details>
```

Keep the one-sentence responsibility and CPU/GPU label visible while details
are collapsed.

- [ ] **Step 5: Add evidence-category legend and source-reading hint**

Add visible badges for:

```html
<span class="fact-tag fact">已实现</span>
<span class="fact-tag measured">已测量</span>
<span class="fact-tag hypothesis">待测假设</span>
<span class="fact-tag unavailable">未实现后端</span>
```

Add a first-reading hint that points from Wrapper to API, Model, Algorithm,
Backend, Postsolve, and KKT.

- [ ] **Step 6: Run the complete static test**

Run:

```bash
/home/yanyu/github/hipdlp2/.venv/bin/python \
  check/test_hipdlp_architecture_html.py
```

Expected: all tests pass and every local source link resolves.

- [ ] **Step 7: Commit the system map**

```bash
git add docs/hipdlp-architecture.html \
  check/test_hipdlp_architecture_html.py
git commit -m "Visualize HiPDLP solve boundaries" \
  -m "Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: Add Data Structures, Ownership, And Residency

**Files:**
- Modify: `check/test_hipdlp_architecture_html.py`
- Modify: `docs/hipdlp-architecture.html`

**Interfaces:**
- Consumes: `architecture-data` and the host/device visual legend.
- Produces: `data.structures`, `data.vectorInventory`, residency cards, and ownership/lifetime explanations used by later estimator and learning views.

- [ ] **Step 1: Add failing data-model assertions**

Add:

```python
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
        self.assertEqual(inventory["candidateDownloadDoubles"], "2*(rows+cols)")
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
```

- [ ] **Step 2: Run the new test and observe the missing JSON keys**

Run:

```bash
/home/yanyu/github/hipdlp2/.venv/bin/python \
  -m unittest \
  check.test_hipdlp_architecture_html.ArchitecturePageTests.test_data_model_matches_current_problem_and_iteration_storage
```

Expected: ERROR or FAIL because `structures` and `vectorInventory` are absent.

- [ ] **Step 3: Extend the embedded fact model**

Add:

```json
"structures": [
  {"id": "highs-lp", "residency": "host", "owner": "HiGHS caller"},
  {"id": "original-model", "residency": "host", "owner": "non-owning view"},
  {"id": "problem-data", "residency": "host", "owner": "model::Problem"},
  {"id": "handles", "residency": "neutral", "owner": "move-only opaque owners"},
  {"id": "cpu-resources", "residency": "host", "owner": "CpuBackend resources"},
  {"id": "cuda-resources", "residency": "device", "owner": "CudaBackend resources"},
  {"id": "host-iterate", "residency": "host", "owner": "algorithm::solve"},
  {"id": "highs-solution", "residency": "host", "owner": "HiGHS result"}
],
"vectorInventory": {
  "columnSized": 8,
  "rowSized": 8,
  "workspaceColumnSized": 1,
  "workspaceRowSized": 1,
  "candidateDownloadDoubles": "2*(rows+cols)"
}
```

- [ ] **Step 4: Build the residency and lifetime view**

Add:

- a transformation chain from `HighsLp` to `HighsSolution`;
- an expandable `Problem::Data` field map;
- separate CPU and CUDA resource columns;
- ownership labels (`owns`, `borrows`, `shared lifetime token`);
- a 16-vector inventory table using the exact names from `IterationVectors`;
- a graph-dependency card explaining weak lifetime dependencies;
- a transfer legend showing setup H2D and candidate/check D2H paths.

Use these source links:

```text
../highs/pdlp/hipdlp/model/problem.hpp
../highs/pdlp/hipdlp/model/problem.cc
../highs/pdlp/hipdlp/backend/types.hpp
../highs/pdlp/hipdlp/backend/cpu/cpu_backend.cc
../highs/pdlp/hipdlp/backend/cuda/cuda_backend.cu
../highs/pdlp/hipdlp/algorithm/solver.cc
```

- [ ] **Step 5: Run all static tests**

Run:

```bash
/home/yanyu/github/hipdlp2/.venv/bin/python \
  check/test_hipdlp_architecture_html.py
```

Expected: all tests pass.

- [ ] **Step 6: Commit data-residency documentation**

```bash
git add docs/hipdlp-architecture.html \
  check/test_hipdlp_architecture_html.py
git commit -m "Document HiPDLP data residency" \
  -m "Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: Visualize The Production Step, Checks, And CUDA Kernels

**Files:**
- Modify: `check/test_hipdlp_architecture_html.py`
- Modify: `docs/hipdlp-architecture.html`

**Interfaces:**
- Consumes: host/device lanes and structure IDs.
- Produces: `iteration` constants, `kernels` records, lane buttons with `data-lane="cpu|cuda|compare"`, and detailed operation/kernel cards.

- [ ] **Step 1: Add failing iteration and launch-configuration tests**

Add:

```python
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
```

- [ ] **Step 2: Run the new tests and observe missing iteration data**

Run:

```bash
/home/yanyu/github/hipdlp2/.venv/bin/python \
  -m unittest \
  check.test_hipdlp_architecture_html.ArchitecturePageTests.test_iteration_and_cuda_constants_match_source \
  check.test_hipdlp_architecture_html.ArchitecturePageTests.test_iteration_view_names_real_backend_and_kernel_symbols \
  check.test_hipdlp_architecture_html.ArchitecturePageTests.test_error_rail_covers_explicit_failure_surfaces
```

Expected: FAIL because `iteration`, `kernels`, and symbol cards do not exist.

- [ ] **Step 3: Add exact iteration and kernel facts**

Extend `architecture-data` with:

```json
"iteration": {
  "spmvPerStep": 2,
  "convergenceCadence": 40,
  "recordedSteps": 39,
  "parameterDoubles": 117,
  "parameterBytes": 936,
  "checkPacketBytes": 112
},
"kernels": [
  {"name": "projectedPrimalKernel", "threads": 256, "maxBlocks": 65535, "dynamicSharedBytes": 0},
  {"name": "projectedDualKernel", "threads": 256, "maxBlocks": 65535, "dynamicSharedBytes": 0},
  {"name": "vectorPairBlendKernel", "threads": 256, "maxBlocks": 65535, "dynamicSharedBytes": 0},
  {"name": "checkPartialsKernel", "threads": 256, "maxBlocks": 1024, "dynamicSharedBytes": 0},
  {"name": "checkPacketKernel", "threads": 256, "blocks": 1, "dynamicSharedBytes": 0},
  {"name": "cusparseSpMV", "grid": "vendor-controlled"}
]
```

- [ ] **Step 4: Add CPU/CUDA/comparison lane controls and the two-SpMV step**

Add three buttons:

```html
<div class="segment" role="group" aria-label="执行路径">
  <button type="button" data-lane="cpu" aria-pressed="false">CPU</button>
  <button type="button" data-lane="cuda" aria-pressed="false">CUDA</button>
  <button type="button" data-lane="compare" aria-pressed="true">对比</button>
</div>
```

Render the five algorithm operations in order and connect them to CPU loops or
CUDA implementations. Make `A * reflected_primal` part of the combined dual
backend operation, and make `A^T * next_dual` the second SpMV. Do not draw
`A * next_primal` inside every production step.

- [ ] **Step 5: Add boundary-check and correctness paths**

Draw:

```text
candidate products
  -> fixed-point preparation
  -> resident check
  -> 112-byte packet
  -> algorithm stopping/PID/restart
```

Add a red fallback branch:

```text
threshold uncertainty / optional metric fallback / suspected certificate
  -> full primal + dual + Ax + A^T y D2H
  -> host oracle
```

Add a distinct final purple path:

```text
postsolve -> original HighsLp -> independent verifyKkt -> OPTIMAL gate
```

- [ ] **Step 6: Add the CUDA kernel drill-down cards**

For each custom kernel, show:

- `<<<blocks, 256, 0, stream>>>`;
- the exact block formula;
- grid-stride behavior;
- buffers read/written;
- branches and error-flag atomics;
- graph-capture participation;
- performance dimension;
- review questions.

For `checkPartialsKernel`, document static shared reduction storage rather than
claiming dynamic shared memory. For `cusparseSpMV`, state that internal launch
geometry requires profiler observation.

- [ ] **Step 7: Add an explicit error rail beside the happy path**

Show all ten tested failure categories. Connect each category to the layer that
creates the structured `api::Error`. Explain that explicit CUDA selection does
not fall back, while `auto` may choose CPU only when CUDA backend creation is
unavailable. Do not turn any failure into a success-colored node.

- [ ] **Step 8: Run all static tests**

Run:

```bash
/home/yanyu/github/hipdlp2/.venv/bin/python \
  check/test_hipdlp_architecture_html.py
```

Expected: all tests pass.

- [ ] **Step 9: Commit the iteration and CUDA view**

```bash
git add docs/hipdlp-architecture.html \
  check/test_hipdlp_architecture_html.py
git commit -m "Visualize HiPDLP CPU and CUDA steps" \
  -m "Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: Add The Performance Workbench And Memory Estimator

**Files:**
- Modify: `check/test_hipdlp_architecture_html.py`
- Modify: `docs/hipdlp-architecture.html`

**Interfaces:**
- Consumes: iteration constants and residency counts.
- Produces: `evidence`, `estimator`, `estimateMemory(input)`, `parseCount(input, label)`, `formatBytes(bytes)`, and estimator output IDs `matrix-bytes`, `vector-bytes`, `candidate-bytes`, `packet-bytes`.

- [ ] **Step 1: Add failing evidence and estimator-data tests**

Add:

```python
    def test_performance_evidence_retains_scope_and_limits(self):
        evidence = {item["id"]: item for item in self.data["evidence"]}
        self.assertEqual(evidence["phase8-gate"]["date"], "2026-09-08")
        self.assertEqual(evidence["cuda-fusion-replay"]["date"], "2026-09-22")
        self.assertEqual(evidence["cuda-fusion-replay"]["reductionPercent"], [41.5, 52.4])
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
```

- [ ] **Step 2: Run the tests and observe missing evidence/estimator data**

Run:

```bash
/home/yanyu/github/hipdlp2/.venv/bin/python \
  -m unittest \
  check.test_hipdlp_architecture_html.ArchitecturePageTests.test_performance_evidence_retains_scope_and_limits \
  check.test_hipdlp_architecture_html.ArchitecturePageTests.test_estimator_constants_and_reference_case
```

Expected: FAIL because the new JSON sections and workbench text are absent.

- [ ] **Step 3: Add scoped evidence records**

Add:

```json
"evidence": [
  {"id": "phase8-gate", "date": "2026-09-08", "releaseGate": true, "status": "blocked"},
  {"id": "two-spmv-cpu", "date": "2026-09-14", "releaseGate": false},
  {"id": "cuda-fusion-replay", "date": "2026-09-22", "releaseGate": false, "reductionPercent": [41.5, 52.4]},
  {"id": "large20", "date": "2026-09-22", "releaseGate": false, "acceptedPairs": 0}
],
"estimator": {
  "doubleBytes": 8,
  "columnVectors": 8,
  "rowVectors": 8,
  "matrixCopies": 2,
  "checkPacketBytes": 112
}
```

- [ ] **Step 4: Build six hotspot cards**

Create cards for:

1. host preprocessing and small-model fixed cost;
2. serial CPU SpMV/fused loops;
3. CUDA setup/uploads/descriptors/graph capture;
4. steady-state CUDA SpMV/fused kernels/replay;
5. cadence boundary, fallback D2H, and synchronization;
6. algorithm trajectory from scaling, norm guard, PID, restart, and KKT.

Each card must have visible headings:

```html
<h4>症状</h4>
<h4>当前证据</h4>
<h4>看什么指标</h4>
<h4>下一步实验</h4>
<h4>不能破坏</h4>
```

Link each measurement card to its dated report under `docs/`.

- [ ] **Step 5: Implement exact estimator validation and formulas**

Add:

```javascript
function parseCount(input, label) {
  const value = Number(input.value);
  if (!Number.isSafeInteger(value) || value < 0) {
    throw new Error(`${label} 必须是非负安全整数`);
  }
  return value;
}

function estimateMemory({ rows, cols, nnz, indexBytes }) {
  const matrixBytes =
    2 * nnz * (8 + indexBytes) +
    ((rows + 1) + (cols + 1)) * indexBytes;
  const vectorBytes = 8 * (rows + cols) * 8;
  const candidateBytes = 2 * (rows + cols) * 8;
  const packetBytes = 112;
  for (const value of [matrixBytes, vectorBytes, candidateBytes]) {
    if (!Number.isSafeInteger(value)) {
      throw new Error("估算结果超过 JavaScript 安全整数范围");
    }
  }
  return { matrixBytes, vectorBytes, candidateBytes, packetBytes };
}
```

Inputs are `rows`, `cols`, `nnz`, and an index-width select containing 4 and 8.
Render errors in a `role="alert"` element. Label all values as lower bounds and
list the exclusions from the approved spec.

- [ ] **Step 6: Map SolveDiagnostics to profiling tools**

Add a table mapping setup/iteration/boundary/download/verification timings,
candidate downloads, ordinary D2H bytes/transfers/synchronizations, host
oracles, restarts, and command-buffer replays to:

- solver trace;
- CPU profiler;
- Nsight Systems;
- Nsight Compute;
- Compute Sanitizer.

- [ ] **Step 7: Run all static tests**

Run:

```bash
/home/yanyu/github/hipdlp2/.venv/bin/python \
  check/test_hipdlp_architecture_html.py
```

Expected: all tests pass.

- [ ] **Step 8: Commit the performance workbench**

```bash
git add docs/hipdlp-architecture.html \
  check/test_hipdlp_architecture_html.py
git commit -m "Add HiPDLP performance workbench" \
  -m "Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 6: Add Guided C++ Review And The Eight-Lesson CUDA Course

**Files:**
- Modify: `check/test_hipdlp_architecture_html.py`
- Modify: `docs/hipdlp-architecture.html`

**Interfaces:**
- Consumes: source links, structures, kernels, and hotspot cards.
- Produces: route IDs `route-30`, `route-cpp-90`, `route-cuda`; lesson IDs `cuda-lesson-1` through `cuda-lesson-8`; 24 lesson question/answer pairs plus four integrated exercise answers; glossary and four-pass review checklist.

- [ ] **Step 1: Teach the parser to count lessons, questions, and answers**

Add fields:

```python
        self.cuda_lessons = 0
        self.review_questions = 0
        self.answers = 0
```

In `handle_starttag`, add:

```python
        classes = set(attributes.get("class", "").split())
        if "cuda-lesson" in classes:
            self.cuda_lessons += 1
        if "review-question" in classes:
            self.review_questions += 1
        if tag == "details" and "answer" in classes:
            self.answers += 1
```

Add the failing test:

```python
    def test_learning_view_has_complete_routes_lessons_and_answers(self):
        self.assertTrue({"route-30", "route-cpp-90", "route-cuda"}
                        <= self.parser.ids)
        self.assertEqual(self.parser.cuda_lessons, 8)
        self.assertEqual(self.parser.review_questions, 28)
        self.assertEqual(self.parser.answers, 28)
        self.assertNotIn('class="answer" open', self.source)
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
```

- [ ] **Step 2: Run the learning-view test and observe zero lessons**

Run:

```bash
/home/yanyu/github/hipdlp2/.venv/bin/python \
  -m unittest \
  check.test_hipdlp_architecture_html.ArchitecturePageTests.test_learning_view_has_complete_routes_lessons_and_answers
```

Expected: FAIL with missing route IDs and counts of zero.

- [ ] **Step 3: Add the 30-minute and 90-minute reading routes**

The 30-minute route follows:

```text
HiPdlpWrapper.cpp
-> api/options.cc + api/solver.cc
-> model/problem.cc
-> algorithm/solver.cc
-> halpern_pdhg_step.cc
-> backend implementation
-> postsolve + kkt_verifier.cc
```

The 90-minute C++ route covers:

- move-only opaque handles;
- virtual backend contract;
- `unique_ptr` and resource ownership;
- phase-guard RAII;
- structured `Error`/`SolveResult`;
- original/transformed model boundary;
- CPU/CUDA contract parity;
- independent verification.

Each stop gets “先猜职责”, “打开源码”, three questions, and a short expected
outcome.

- [ ] **Step 4: Add all eight CUDA lessons**

Use IDs `cuda-lesson-1` through `cuda-lesson-8`. Each lesson must contain:

- one concept diagram or compact visual;
- one direct HiPDLP example;
- exactly three elements with class `review-question`;
- exactly three collapsed `<details class="answer">`;
- one checklist.

Answer shape:

```html
<details class="answer">
  <summary>展开参考答案</summary>
  <p><strong>结论：</strong>每个线程处理一个起始索引，并按整个 grid 的 stride 继续前进。</p>
  <p><strong>代码证据：</strong><code>index += blockDim.x * gridDim.x</code>。</p>
  <p><strong>Review 推理：</strong>静态代码能证明覆盖方式；实际吞吐和 occupancy 仍需 profiler。</p>
</details>
```

Cover the eight approved topics without collapsing two lessons into one.

- [ ] **Step 5: Add integrated exercises and four-pass review mode**

Integrated exercises:

1. trace `primal_candidate` from allocation through kernel, resident check,
   download, and postsolve;
2. decide whether a full D2H is needed for a normal check;
3. prove grid-stride coverage when `N > 65535 * 256`;
4. check whether graph dependencies can expire.

Mark each integrated exercise with `class="review-question"` and follow it
with one collapsed `<details class="answer">` using the same conclusion,
code-evidence, and review-reasoning structure as lesson answers.

Review passes:

1. interface and ownership;
2. numerical correctness;
3. concurrency and CUDA semantics;
4. performance.

Add “good signal”, “danger signal”, and reusable human/AI review prompts to
each pass.

- [ ] **Step 6: Add the glossary**

Include at least:

```text
CSR, SpMV, primal, dual, reflection, Halpern anchor, resident,
H2D, D2H, pinned memory, stream, event, warp, occupancy,
coalescing, divergence, graph capture, deterministic reduction,
host oracle, postsolve, KKT
```

- [ ] **Step 7: Run all static tests**

Run:

```bash
/home/yanyu/github/hipdlp2/.venv/bin/python \
  check/test_hipdlp_architecture_html.py
```

Expected: all tests pass with eight lessons, 28 questions, and 28 collapsed
answers.

- [ ] **Step 8: Commit the learning and review material**

```bash
git add docs/hipdlp-architecture.html \
  check/test_hipdlp_architecture_html.py
git commit -m "Add HiPDLP C++ and CUDA review guide" \
  -m "Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 7: Finish Progressive Enhancement, Accessibility, And Browser Validation

**Files:**
- Modify: `check/test_hipdlp_architecture_html.py`
- Modify: `docs/hipdlp-architecture.html`

**Interfaces:**
- Consumes: all five populated views, lane controls, estimator, lessons, and source links.
- Produces: `activateView(viewId, updateHash)`, `activateLane(lane)`, hash navigation, keyboard tab navigation, evidence filtering, source-detail interactions, reduced-motion/print behavior, and final browser-verified artifact.

- [ ] **Step 1: Add failing accessibility and progressive-enhancement tests**

Add:

```python
    def test_tabs_reference_existing_views_and_static_content_is_not_hidden(self):
        controls = {tab["aria-controls"] for tab in self.parser.tabs}
        self.assertEqual(controls, {
            "system-map",
            "data-residency",
            "iteration-correctness",
            "bottlenecks",
            "learning-review",
        })
        self.assertTrue(controls <= self.parser.ids)
        self.assertNotIn('<section class="view" id="system-map" hidden', self.source)
        self.assertIn("prefers-reduced-motion: reduce", self.source)
        self.assertIn("@media print", self.source)
        self.assertIn("aria-live=\"polite\"", self.source)

    def test_interaction_functions_are_present_without_external_frameworks(self):
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
            self.assertIn(token, self.source)
        for forbidden in ("React", "Vue", "d3.", "mermaid", "cdn."):
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
```

- [ ] **Step 2: Run the tests and observe missing interaction/accessibility tokens**

Run:

```bash
/home/yanyu/github/hipdlp2/.venv/bin/python \
  -m unittest \
  check.test_hipdlp_architecture_html.ArchitecturePageTests.test_tabs_reference_existing_views_and_static_content_is_not_hidden \
  check.test_hipdlp_architecture_html.ArchitecturePageTests.test_interaction_functions_are_present_without_external_frameworks
```

Expected: FAIL until the JavaScript and media styles are added.

- [ ] **Step 3: Implement progressive tab navigation**

At startup, add an `enhanced` class to `<html>`, then:

```javascript
function activateView(viewId, updateHash = true) {
  const target = document.getElementById(viewId);
  if (!target || !target.classList.contains("view")) return;
  views.forEach(view => {
    view.hidden = view !== target;
  });
  tabs.forEach(tab => {
    const selected = tab.getAttribute("aria-controls") === viewId;
    tab.setAttribute("aria-selected", String(selected));
    tab.tabIndex = selected ? 0 : -1;
  });
  if (updateHash) history.replaceState(null, "", `#${viewId}`);
}
```

Implement ArrowLeft, ArrowRight, Home, and End for the tab list. Use the hash
only when it matches one of the five view IDs. Without JavaScript, all sections
remain visible.

- [ ] **Step 4: Implement lane and evidence filters**

`activateLane("cpu"|"cuda"|"compare")` must:

- set `aria-pressed`;
- hide only cards explicitly tagged for the other backend;
- retain shared algorithm/control cards;
- update a polite live-region summary.

Evidence filters must preserve each card's date, classification, and limitation
text rather than rewriting factual content.

- [ ] **Step 5: Complete responsive, reduced-motion, focus, and print CSS**

Add:

- visible `:focus-visible`;
- a stacked layout below 760px;
- intentionally scrollable SVG wrappers rather than page-level overflow;
- `prefers-reduced-motion: reduce` disabling transitions and animated flow;
- print CSS that shows all five views, all source facts, and all answer content;
- sufficient non-color labels for every lane and transfer.

- [ ] **Step 6: Run the complete Python test suite for the page**

Run:

```bash
/home/yanyu/github/hipdlp2/.venv/bin/python \
  check/test_hipdlp_architecture_html.py
```

Expected: all tests pass with no failure or error.

- [ ] **Step 7: Open the page directly and validate navigation**

Open:

```text
file:///home/yanyu/github/hipdlp2/docs/hipdlp-architecture.html
```

In the integrated browser:

1. click each of the five tabs;
2. use ArrowLeft/ArrowRight/Home/End while a tab has focus;
3. reload each of the five hash URLs and confirm the expected view opens;
4. confirm CPU, CUDA, and comparison lane buttons change only backend-specific
   cards;
5. confirm source links resolve to local files.

Expected: no console error and no network request.

- [ ] **Step 8: Validate the estimator with the reference case and invalid input**

Enter:

```text
rows = 1000
cols = 500
nnz = 10000
index width = 8 bytes
```

Expected outputs:

```text
A + A^T core storage = 332016 bytes
16 iteration vectors = 96000 bytes
full candidate D2H = 24000 bytes
normal check packet = 112 bytes
```

Then enter `rows = -1`.

Expected: a visible alert saying the value must be a nonnegative safe integer;
no success-shaped numeric output remains visible.

- [ ] **Step 9: Validate lessons, narrow layout, reduced motion, and print**

Check:

- all eight lessons;
- answers collapsed initially and keyboard-toggleable;
- 390x844 viewport without page-level horizontal scrolling;
- reduced-motion emulation without animated flow;
- print preview showing all views and answer content;
- source/evidence labels still understandable in grayscale.

- [ ] **Step 10: Perform the final factual and repository-scope audit**

Run:

```bash
git --no-pager diff --check
git --no-pager status --short
/home/yanyu/github/hipdlp2/.venv/bin/python \
  check/test_hipdlp_architecture_html.py
```

Compare the rendered facts against:

```text
highs/pdlp/hipdlp/algorithm/solver.cc
highs/pdlp/hipdlp/algorithm/halpern_pdhg_step.cc
highs/pdlp/hipdlp/backend/cuda/cuda_kernels.cu
highs/pdlp/hipdlp/backend/cuda/cuda_backend.cu
highs/pdlp/hipdlp/model/problem.cc
highs/pdlp/hipdlp/verification/kkt_verifier.cc
docs/hipdlp-phase8-gpu-performance-report.md
docs/hipdlp-phase8-step-fusion.md
docs/hipdlp-phase8-large20-screen.md
```

Expected: only the two planned implementation files differ from the
pre-existing dirty state; tests pass; no whitespace errors; browser checks
pass; every dated measurement keeps its scope limitation.

- [ ] **Step 11: Commit the completed explorer**

```bash
git add docs/hipdlp-architecture.html \
  check/test_hipdlp_architecture_html.py
git commit -m "Complete HiPDLP architecture explorer" \
  -m "Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```
