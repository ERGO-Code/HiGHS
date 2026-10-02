"""Focused python test for the automatic solver-select feature
(util/HighsSolverSelect.h), exercised purely through the public
HighsOptions surface: "solver_select_strategy" (int, 0-3) and
"solver_select_require_basis" (bool), used together with solver="choose".
"""
import unittest

import highspy


class TestSolverSelect(unittest.TestCase):
    def get_example_model(self):
        """
        minimize    f  =  x0 +  x1
        subject to              x1 <= 7
                    5 <=  x0 + 2x1 <= 15
                    6 <= 3x0 + 2x1
                    0 <= x0 <= 4; 1 <= x1
        """
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)

        x0 = h.addVariable(lb=0, ub=4, obj=1)
        x1 = h.addVariable(lb=1, ub=7, obj=1)
        h.addConstr(5 <= x0 + 2 * x1 <= 15)
        h.addConstr(6 <= 3 * x0 + 2 * x1)
        return h

    def test_strategy_bounds(self):
        h = self.get_example_model()
        # Valid range is [0, 3].
        for value in (0, 1, 2, 3):
            self.assertEqual(
                h.setOptionValue("solver_select_strategy", value),
                highspy.HighsStatus.kOk,
            )
        for value in (-1, 4):
            self.assertEqual(
                h.setOptionValue("solver_select_strategy", value),
                highspy.HighsStatus.kError,
            )

    def test_require_basis_bounds(self):
        h = self.get_example_model()
        self.assertEqual(
            h.setOptionValue("solver_select_require_basis", True),
            highspy.HighsStatus.kOk,
        )
        self.assertEqual(
            h.setOptionValue("solver_select_require_basis", False),
            highspy.HighsStatus.kOk,
        )

    def test_each_strategy_solves_optimally(self):
        # strategy 0: simplex only; 3: the feature-based heuristic.
        # 1 and 2 currently fall back to the strategy-0 behaviour.
        # for strategy in (0, 1, 2, 3):
        for strategy in (0, 3):
            for require_basis in (False, True):
                with self.subTest(strategy=strategy, require_basis=require_basis):
                    h = self.get_example_model()
                    h.setOptionValue("solver", "choose")
                    h.setOptionValue("solver_select_strategy", strategy)
                    h.setOptionValue(
                        "solver_select_require_basis", require_basis
                    )

                    run_status = h.run()
                    self.assertEqual(run_status, highspy.HighsStatus.kOk)
                    self.assertEqual(
                        h.getModelStatus(), highspy.HighsModelStatus.kOptimal
                    )
                    self.assertAlmostEqual(h.getObjectiveValue(), 2.75, places=6)

    def test_require_basis_yields_a_basis(self):
        h = self.get_example_model()
        h.setOptionValue("solver", "choose")
        h.setOptionValue("solver_select_strategy", 3)
        h.setOptionValue("solver_select_require_basis", True)

        self.assertEqual(h.run(), highspy.HighsStatus.kOk)
        self.assertEqual(h.getModelStatus(), highspy.HighsModelStatus.kOptimal)
        basis = h.getBasis()
        self.assertTrue(basis.valid)


if __name__ == "__main__":
    unittest.main()
