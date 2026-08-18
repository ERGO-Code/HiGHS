"""Functional smoke test for the pyodide/wasm32 highspy wheel.

The excluded pytest cases (test_user_interrupts, test_with_context_manager)
only cover Python-level threading/signal behaviour that Pyodide doesn't
support. This script instead checks that the solver itself -- LP and, more
importantly, MIP branch-and-bound (the part of HiGHS that drives the
parallel task scheduler) -- still produces correct results when compiled
single-threaded for wasm32.
"""

import sys

import highspy


def check(label, condition):
    print(f"{'PASS' if condition else 'FAIL'}: {label}")
    if not condition:
        sys.exit(1)


# LP: minimize x2 s.t. x2 - x1 >= 2, x1 + x2 >= 0 (unbounded x1)
# optimum at x1=-1, x2=1
h = highspy.Highs()
h.silent()

x1 = h.addVariable(lb=-h.inf)
x2 = h.addVariable(lb=-h.inf)
h.addConstrs(x2 - x1 >= 2, x1 + x2 >= 0)
h.minimize(x2)

check("LP solved to optimality", h.getModelStatus() == highspy.HighsModelStatus.kOptimal)
check("LP objective correct", abs(h.getObjectiveValue() - 1.0) < 1e-6)

# MIP (binary knapsack): max 8x1 + 5x2 + 3x3 + 11x4 + 7x5
#                         s.t. 4x1 + 3x2 + x3 + 5x4 + 4x5 <= 11
# known optimum: x = [1, 0, 1, 1, 0], objective = 22
h2 = highspy.Highs()
h2.silent()

x = h2.addBinaries(5)
h2.addConstr((x * [4, 3, 1, 5, 4]).sum() <= 11)
h2.maximize((x * [8, 5, 3, 11, 7]).sum())

check("MIP solved to optimality", h2.getModelStatus() == highspy.HighsModelStatus.kOptimal)
check("MIP objective correct", abs(h2.getObjectiveValue() - 22.0) < 1e-6)
solution = [round(abs(v)) for v in h2.val(x)]
check(f"MIP solution correct ({solution})", solution == [1, 0, 1, 1, 0])

print("All wasm smoke-test checks passed.")
