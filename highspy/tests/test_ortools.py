import subprocess
import sys
import tempfile
import unittest
from importlib.util import find_spec

IMPORT_HIGHSPY = "import highspy\n"
IMPORT_ORTOOLS = "from ortools.linear_solver import pywraplp\n"
SOLVE_WITH_BOTH = """
h = highspy.Highs()
h.setOptionValue("output_flag", False)
x = h.addVariable(lb=0, ub=4)
h.minimize(-x)
assert abs(h.getObjectiveValue() + 4) < 1e-6

solver = pywraplp.Solver.CreateSolver("HIGHS")
assert solver is not None
v = solver.NumVar(0, 4, "x")
solver.Minimize(-v)
assert solver.Solve() == solver.OPTIMAL
assert abs(solver.Objective().Value() + 4) < 1e-6
"""


@unittest.skipIf(find_spec("ortools") is None, "ortools is not installed")
class TestHighsPyWithOrTools(unittest.TestCase):
    """OR-Tools bundles its own HiGHS, which must not be mixed up with highspy's"""

    def run_in_new_process(self, imports):
        # Libraries stay loaded for the life of a process, so each import order
        # needs its own. The temporary directory stops Python from finding the
        # highspy source directory when the tests are run from the repository.
        with tempfile.TemporaryDirectory() as cwd:
            result = subprocess.run(
                [sys.executable, "-c", imports + SOLVE_WITH_BOTH],
                cwd=cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                timeout=300,
            )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_highspy_imported_first(self):
        self.run_in_new_process(IMPORT_HIGHSPY + IMPORT_ORTOOLS)

    def test_ortools_imported_first(self):
        self.run_in_new_process(IMPORT_ORTOOLS + IMPORT_HIGHSPY)
