import tempfile
import unittest
import highspy
from highspy.highs import highs_linear_expression, qsum
import numpy as np
from sys import platform
import signal


class TestHighsPy(unittest.TestCase):
    def get_example_model(self):
        """
        minimize    f  =  x0 +  x1
        subject to              x1 <= 7
                    5 <=  x0 + 2x1 <= 15
                    6 <= 3x0 + 2x1
                    0 <= x0 <= 4; 1 <= x1
        """
        inf = highspy.kHighsInf
        h = highspy.Highs()
        # Define a HighsLp instance
        lp = highspy.HighsLp()
        lp.num_col_ = 2
        lp.num_row_ = 3
        lp.col_cost_ = np.array([1, 1], dtype=np.double)
        lp.col_lower_ = np.array([0, 1], dtype=np.double)
        lp.col_upper_ = np.array([4, inf], dtype=np.double)
        lp.row_lower_ = np.array([-inf, 5, 6], dtype=np.double)
        lp.row_upper_ = np.array([7, 15, inf], dtype=np.double)
        lp.a_matrix_.start_ = np.array([0, 2, 5])
        lp.a_matrix_.index_ = np.array([1, 2, 0, 1, 2])
        lp.a_matrix_.value_ = np.array([1, 3, 1, 2, 2], dtype=np.double)
        h.passModel(lp)
        return h

    def test_example_model_builder(self):
        """
        minimize    f  =  x0 +  x1
        subject to              x1 <= 7
                    5 <=  x0 + 2x1 <= 15
                    6 <= 3x0 + 2x1
                    0 <= x0 <= 4; 1 <= x1
        """
        h = highspy.Highs()

        x0 = h.addVariable(lb=0, ub=4, obj=1)
        x1 = h.addVariable(lb=1, ub=7, obj=1)

        h.addConstr(5 <= x0 + 2 * x1 <= 15)
        h.addConstr(6 <= 3 * x0 + 2 * x1)

        lp = h.getLp()

        self.assertEqual(lp.num_col_, 2)
        self.assertEqual(lp.num_row_, 2)
        self.assertAlmostEqual(lp.col_cost_[0], 1)
        self.assertAlmostEqual(lp.col_lower_[0], 0)
        self.assertAlmostEqual(lp.col_upper_[0], 4)
        self.assertAlmostEqual(lp.row_lower_[0], 5)
        self.assertAlmostEqual(lp.row_upper_[0], 15)
        self.assertAlmostEqual(lp.row_lower_[1], 6)
        self.assertAlmostEqual(lp.row_upper_[1], highspy.kHighsInf)

    def test_hipo(self):
        print("running hipo test")
        h = self.get_example_model()
        h.setOptionValue("solver", "hipo")
        h.setOptionValue("output_flag", True)

        [status, output_flag] = h.getOptionValue("solver")
        self.assertEqual(output_flag, "hipo")

        h.run()
        self.assertEqual(status, highspy.HighsStatus.kOk)

        status = h.getModelStatus()
        self.assertEqual(status, highspy.HighsModelStatus.kOptimal)