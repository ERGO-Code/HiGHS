import tempfile
import unittest
import highspy
from highspy.highs import highs_linear_expression, qsum
import numpy as np
from sys import platform
import signal


class TestHighsPy(unittest.TestCase):
    def assertEqualExpr(self, expr, idxs, vals, constant=None, bounds=None):
        self.assertEqual(list(map(int, expr.idxs)), list(map(int, idxs)), "variable index")
        self.assertEqual(expr.vals, vals, "variable values")
        self.assertEqual(expr.constant, constant, "constant")
        self.assertEqual(expr.bounds, (bounds[0], bounds[1]) if bounds is not None else None, "bounds")

    def get_basic_model(self):
        """
        min y
        s.t.
        -x + y >= 2
         x + y >= 0
        """
        inf = highspy.kHighsInf
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)
        h.addVars(2, np.array([-inf, -inf]), np.array([inf, inf]))
        h.changeColsCost(2, np.array([0, 1]), np.array([0, 1], dtype=np.double))
        num_cons = 2
        lower = np.array([2, 0], dtype=np.double)
        upper = np.array([inf, inf], dtype=np.double)
        num_new_nz = 4
        starts = np.array([0, 2])
        indices = np.array([0, 1, 0, 1])
        values = np.array([-1, 1, 1, 1], dtype=np.double)
        h.addRows(num_cons, lower, upper, num_new_nz, starts, indices, values)
        return h

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

    def get_infeasible_model(self):
        inf = highspy.kHighsInf
        lp = highspy.HighsLp()
        lp.num_col_ = 2
        lp.num_row_ = 2
        lp.col_cost_ = np.array([10, 15], dtype=np.double)
        lp.col_lower_ = np.array([0, 0], dtype=np.double)
        lp.col_upper_ = np.array([inf, inf], dtype=np.double)
        lp.row_lower_ = np.array([3, 1], dtype=np.double)
        lp.row_upper_ = np.array([3, 1], dtype=np.double)
        lp.a_matrix_.start_ = np.array([0, 2, 4])
        lp.a_matrix_.index_ = np.array([0, 1, 0, 1])
        lp.a_matrix_.value_ = np.array([2, 1, 1, 3], dtype=np.double)
        lp.offset_ = 0
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)
        status = h.passModel(lp)
        self.assertEqual(status, highspy.HighsStatus.kOk)
        h.setOptionValue("presolve", "off")
        return h

    def test_version(self):
        h = self.get_basic_model()
        self.assertEqual(h.version(), "1.7.2")
        self.assertEqual(h.versionMajor(), 1)
        self.assertEqual(h.versionMinor(), 7)
        self.assertEqual(h.versionPatch(), 2)

    def test_basics(self):
        h = self.get_basic_model()
        h.passColName(0, "Col0")
        h.passColName(1, "Col1")
        h.passRowName(0, "Row0")
        h.passRowName(1, "Row1")
        #        h.setOptionValue('output_flag', True)
        h.writeModel("")
        h.setOptionValue("output_flag", False)
        self.assertEqual(h.setOptionValue("presolve", "off"), highspy.HighsStatus.kOk)
        #        h.setOptionValue('output_flag', True)
        h.run()

        # Info can be obtained from the class instance, specific call
        # and, in the case of objective_function_value,
        # h.getObjectiveValue()
        info = h.getInfo()
        objective_function_value0 = info.objective_function_value
        self.assertAlmostEqual(objective_function_value0, 1)
        [status, objective_function_value1] = h.getInfoValue("objective_function_value")
        self.assertAlmostEqual(objective_function_value0, objective_function_value1)
        self.assertAlmostEqual(h.getObjectiveValue(), objective_function_value0)

        simplex_iteration_count0 = info.simplex_iteration_count
        self.assertAlmostEqual(simplex_iteration_count0, 2)
        [status, simplex_iteration_count1] = h.getInfoValue("simplex_iteration_count")
        self.assertAlmostEqual(simplex_iteration_count0, simplex_iteration_count1)

        sol = h.getSolution()
        self.assertAlmostEqual(sol.col_value[0], -1)
        self.assertAlmostEqual(sol.col_value[1], 1)

        #        h.setOptionValue('output_flag', False)
        """
        min y
        s.t.
        -x + y >= 3
        x + y >= 0
        """
        inf = highspy.kHighsInf
        h.changeRowBounds(0, 3, inf)
        h.run()
        sol = h.getSolution()
        self.assertAlmostEqual(sol.col_value[0], -1.5)
        self.assertAlmostEqual(sol.col_value[1], 1.5)

        # now make y integer
        h.changeColsIntegrality(1, np.array([1]), np.array([highspy.HighsVarType.kInteger]))
        h.run()
        sol = h.getSolution()
        self.assertAlmostEqual(sol.col_value[0], -1.5)
        self.assertAlmostEqual(sol.col_value[1], 2)

        """
        now delete the first constraint and add a new one
        
        min y
        s.t.
        x + y >= 0
        -x + y >= 0
        """
        h.deleteRows(1, np.array([0]))
        h.addRows(1, np.array([0], dtype=np.double), np.array([inf]), 2, np.array([0]), np.array([0, 1]), np.array([-1, 1], dtype=np.double))
        h.run()
        sol = h.getSolution()
        self.assertAlmostEqual(sol.col_value[0], 0)
        self.assertAlmostEqual(sol.col_value[1], 0)

        # change the upper bound of x to -5
        h.changeColsBounds(1, np.array([0]), np.array([-inf], dtype=np.double), np.array([-5], dtype=np.double))
        h.run()
        sol = h.getSolution()
        self.assertAlmostEqual(sol.col_value[0], -5)
        self.assertAlmostEqual(sol.col_value[1], 5)

        # now maximize
        h.changeColCost(1, -1)
        h.changeRowBounds(0, -inf, 0)
        h.changeRowBounds(1, -inf, 0)
        h.run()
        sol = h.getSolution()
        self.assertAlmostEqual(sol.col_value[0], -5)
        self.assertAlmostEqual(sol.col_value[1], -5)

        h.changeColCost(1, 1)
        [status, sense] = h.getObjectiveSense()
        self.assertEqual(sense, highspy.ObjSense.kMinimize)
        h.changeObjectiveSense(highspy.ObjSense.kMaximize)
        [status, sense] = h.getObjectiveSense()
        self.assertEqual(sense, highspy.ObjSense.kMaximize)
        h.run()
        sol = h.getSolution()
        self.assertAlmostEqual(sol.col_value[0], -5)
        self.assertAlmostEqual(sol.col_value[1], -5)

        self.assertAlmostEqual(h.getObjectiveValue(), -5)

        h.changeObjectiveOffset(1)
        [status, offset] = h.getObjectiveOffset()
        self.assertAlmostEqual(offset, 1)
        h.run()
        self.assertAlmostEqual(h.getObjectiveValue(), -4)

        info = h.getInfo()
        mip_node_count0 = info.mip_node_count
        self.assertAlmostEqual(mip_node_count0, 0)
        [status, mip_node_count1] = h.getInfoValue("mip_node_count")
        self.assertEqual(status, highspy.HighsStatus.kOk)
        self.assertAlmostEqual(mip_node_count0, mip_node_count1)

    def test_example(self):
        h = self.get_example_model()
        lp = h.getLp()
        #
        # Extract column 0
        iCol = 0
        [status, cost, lower, upper, get_num_nz] = h.getCol(iCol)
        self.assertEqual(cost, lp.col_cost_[iCol])
        self.assertEqual(lower, lp.col_lower_[iCol])
        self.assertEqual(upper, lp.col_upper_[iCol])
        index = np.empty(get_num_nz)
        value = np.empty(get_num_nz, dtype=np.double)
        [status, index, value] = h.getColEntries(iCol)
        for iEl in range(get_num_nz):
            self.assertEqual(index[iEl], lp.a_matrix_.index_[iEl])
            self.assertEqual(value[iEl], lp.a_matrix_.value_[iEl])
        #
        # Extract columns 0 and 1
        indices = np.array([0, 1])
        [status, get_num_col, cost, lower, upper, get_num_nz] = h.getCols(2, indices)
        for get_col in range(get_num_col):
            iCol = indices[get_col]
            self.assertEqual(cost[get_col], lp.col_cost_[iCol])
            self.assertEqual(lower[get_col], lp.col_lower_[iCol])
            self.assertEqual(upper[get_col], lp.col_upper_[iCol])
        start = np.empty(get_num_col)
        index = np.empty(get_num_nz)
        value = np.empty(get_num_nz, dtype=np.double)
        [status, start, index, value] = h.getColsEntries(2, indices)
        for iCol in range(lp.num_col_):
            self.assertEqual(start[iCol], lp.a_matrix_.start_[iCol])
        for iEl in range(get_num_nz):
            self.assertEqual(index[iEl], lp.a_matrix_.index_[iEl])
            self.assertEqual(value[iEl], lp.a_matrix_.value_[iEl])
        #
        # Extract row 1
        iRow = 1
        [status, lower, upper, get_num_nz] = h.getRow(iRow)
        self.assertEqual(lower, lp.row_lower_[iRow])
        self.assertEqual(upper, lp.row_upper_[iRow])
        index = np.empty(get_num_nz)
        value = np.empty(get_num_nz, dtype=np.double)
        [status, index, value] = h.getRowEntries(iRow)
        #
        # Extract rows 0 and 2
        indices = np.array([0, 2])
        [status, get_num_row, lower, upper, get_num_nz] = h.getRows(2, indices)
        for get_row in range(get_num_row):
            iRow = indices[get_row]
            self.assertEqual(lower[get_row], lp.row_lower_[iRow])
            self.assertEqual(upper[get_row], lp.row_upper_[iRow])
        start = np.empty(get_num_row)
        index = np.empty(get_num_nz)
        value = np.empty(get_num_nz, dtype=np.double)
        [status, start, index, value] = h.getRowsEntries(2, indices)

    def test_options(self):
        h = highspy.Highs()

        # test vanilla get option value method

        [status, output_flag] = h.getOptionValue("output_flag")
        [status, solver] = h.getOptionValue("solver")
        [status, primal_feasibility_tolerance] = h.getOptionValue("primal_feasibility_tolerance")
        [status, simplex_update_limit] = h.getOptionValue("simplex_update_limit")
        self.assertEqual(output_flag, True)
        self.assertEqual(solver, "choose")
        self.assertEqual(primal_feasibility_tolerance, 1e-7)
        self.assertEqual(simplex_update_limit, 5000)
        # Illegal name
        option_value = h.getOptionValue("simplex_limit")
        self.assertEqual(option_value[0], highspy.HighsStatus.kError)

        # test bool option
        [status, type] = h.getOptionType("output_flag")
        self.assertEqual(type, highspy.HighsOptionType.kBool)

        h.setOptionValue("output_flag", True)
        [status, value] = h.getOptionValue("output_flag")
        self.assertTrue(value)
        h.setOptionValue("output_flag", False)
        [status, value] = h.getOptionValue("output_flag")
        self.assertFalse(value)

        # test string option
        [status, type] = h.getOptionType("presolve")
        self.assertEqual(type, highspy.HighsOptionType.kString)
        h.setOptionValue("presolve", "off")
        [status, value] = h.getOptionValue("presolve")
        self.assertEqual(value, "off")
        h.setOptionValue("presolve", "on")
        [status, value] = h.getOptionValue("presolve")
        self.assertEqual(value, "on")

        # test int option
        [status, type] = h.getOptionType("threads")
        self.assertEqual(type, highspy.HighsOptionType.kInt)
        h.setOptionValue("threads", 1)
        [status, value] = h.getOptionValue("threads")
        self.assertEqual(value, 1)
        h.setOptionValue("threads", 2)
        [status, value] = h.getOptionValue("threads")
        self.assertEqual(value, 2)

        # test double option
        [status, type] = h.getOptionType("time_limit")
        self.assertEqual(type, highspy.HighsOptionType.kDouble)
        h.setOptionValue("time_limit", 1.7)
        [status, value] = h.getOptionValue("time_limit")
        self.assertAlmostEqual(value, 1.7)
        h.setOptionValue("time_limit", 2.7)
        [status, value] = h.getOptionValue("time_limit")
        self.assertAlmostEqual(value, 2.7)

    def test_clear(self):
        h = self.get_basic_model()
        self.assertEqual(h.getNumCol(), 2)
        self.assertEqual(h.getNumRow(), 2)
        self.assertEqual(h.getNumNz(), 4)

        [status, orig_feas_tol] = h.getOptionValue("primal_feasibility_tolerance")
        new_feas_tol = orig_feas_tol + 1
        h.setOptionValue("primal_feasibility_tolerance", new_feas_tol)
        [status, value] = h.getOptionValue("primal_feasibility_tolerance")
        self.assertAlmostEqual(value, new_feas_tol)
        h.clear()
        self.assertEqual(h.getNumCol(), 0)
        self.assertEqual(h.getNumRow(), 0)
        self.assertEqual(h.getNumNz(), 0)
        [status, value] = h.getOptionValue("primal_feasibility_tolerance")
        self.assertAlmostEqual(value, orig_feas_tol)

        h = self.get_basic_model()
        h.setOptionValue("primal_feasibility_tolerance", new_feas_tol)
        [status, value] = h.getOptionValue("primal_feasibility_tolerance")
        self.assertAlmostEqual(value, new_feas_tol)
        h.clearModel()
        self.assertEqual(h.getNumCol(), 0)
        self.assertEqual(h.getNumRow(), 0)
        self.assertEqual(h.getNumNz(), 0)
        [status, value] = h.getOptionValue("primal_feasibility_tolerance")
        self.assertAlmostEqual(value, new_feas_tol)

        h = self.get_basic_model()
        h.run()
        sol = h.getSolution()
        self.assertAlmostEqual(sol.col_value[0], -1)
        self.assertAlmostEqual(sol.col_value[1], 1)
        h.clearSolver()
        self.assertEqual(h.getNumCol(), 2)
        self.assertEqual(h.getNumRow(), 2)
        self.assertEqual(h.getNumNz(), 4)
        sol = h.getSolution()
        self.assertFalse(sol.value_valid)
        self.assertFalse(sol.dual_valid)

        h = self.get_basic_model()
        [status, orig_feas_tol] = h.getOptionValue("primal_feasibility_tolerance")
        new_feas_tol = orig_feas_tol + 1
        h.setOptionValue("primal_feasibility_tolerance", new_feas_tol)
        [status, value] = h.getOptionValue("primal_feasibility_tolerance")
        self.assertAlmostEqual(value, new_feas_tol)
        h.resetOptions()
        [status, value] = h.getOptionValue("primal_feasibility_tolerance")
        self.assertAlmostEqual(value, orig_feas_tol)

    def test_ranging(self):
        inf = highspy.kHighsInf
        h = self.get_basic_model()
        # Cost ranging
        # c0 2 -1 1 0
        # c1 0 0 inf inf
        #
        ## Bound ranging
        ## Columns
        # c0 1 -inf inf 1
        # c1 1 1 inf 1
        ## Rows
        # r0 -inf -inf inf inf
        # r1 -inf -inf inf inf
        h.run()
        [status, ranging] = h.getRanging()
        self.assertEqual(ranging.col_cost_dn.objective_[0], 2)
        self.assertEqual(ranging.col_cost_dn.value_[0], -1)
        self.assertEqual(ranging.col_cost_up.value_[0], 1)
        self.assertEqual(ranging.col_cost_up.objective_[0], 0)
        self.assertEqual(ranging.col_cost_dn.objective_[1], 0)
        self.assertEqual(ranging.col_cost_dn.value_[1], 0)
        self.assertEqual(ranging.col_cost_up.value_[1], inf)
        self.assertEqual(ranging.col_cost_up.objective_[1], inf)
        #
        self.assertEqual(ranging.col_bound_dn.objective_[0], 1)
        self.assertEqual(ranging.col_bound_dn.value_[0], -inf)
        self.assertEqual(ranging.col_bound_up.value_[0], inf)
        self.assertEqual(ranging.col_bound_up.objective_[0], 1)
        self.assertEqual(ranging.col_bound_dn.objective_[1], 1)
        self.assertEqual(ranging.col_bound_dn.value_[1], 1)
        self.assertEqual(ranging.col_bound_up.value_[1], inf)
        self.assertEqual(ranging.col_bound_up.objective_[1], 1)
        #
        self.assertEqual(ranging.row_bound_dn.objective_[0], -inf)
        self.assertEqual(ranging.row_bound_dn.value_[0], -inf)
        self.assertEqual(ranging.row_bound_up.value_[0], inf)
        self.assertEqual(ranging.row_bound_up.objective_[0], inf)
        self.assertEqual(ranging.row_bound_dn.objective_[1], -inf)
        self.assertEqual(ranging.row_bound_dn.value_[1], -inf)
        self.assertEqual(ranging.row_bound_up.value_[1], inf)
        self.assertEqual(ranging.row_bound_up.objective_[1], inf)

    # this catches error in our highs_bindings (pybind11) when passed arrays are not contiguous
    def test_numpy_slice(self):
        h = highspy.Highs()

        N = 10
        zero, ones = np.zeros(N), np.ones(N)
        h.addCols(N, zero, zero, ones, 0, [], [], [])

        x = np.arange(N, dtype=np.int32)  # [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        tmp = x[1::2]  # [1, 3, 5, 7, 9]

        h.addRow(1, 1, len(tmp), tmp, ones)

        # needs to be rowwise for test to work correctly
        self.assertEqual(h.getLp().a_matrix_.format_, highspy.MatrixFormat.kRowwise)
        self.assertEqual(h.getLp().a_matrix_.index_, list(tmp))

    def test_constraint_removal(self):
        h = highspy.Highs()
        x = h.addVariable(lb=-highspy.kHighsInf)
        y = h.addVariable(lb=-highspy.kHighsInf)
        c1 = h.addConstr(-x + y >= 2)
        c2 = h.addConstr(x + y >= 0)
        self.assertEqual(h.numConstrs, 2)
        h.removeConstr(c1)
        self.assertEqual(h.numConstrs, 1)

    def test_infeasible_model(self):
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)
        h.setOptionValue("presolve", "off")

        x = h.addVariable()
        y = h.addVariable()

        c1 = h.addConstr(x + y == 3)
        c2 = h.addConstr(x + y == 1)

        status = h.minimize(10 * x + 15 * y)
        self.assertEqual(status, highspy.HighsStatus.kOk)

        status = h.getModelStatus()
        self.assertEqual(status, highspy.HighsModelStatus.kInfeasible)

        # change the model to be feasible
        h.chgCoeff(c1, y, 3)

        status = h.run()
        self.assertEqual(status, highspy.HighsStatus.kOk)

        status = h.getModelStatus()
        self.assertEqual(status, highspy.HighsModelStatus.kOptimal)

    def test_simple_basics_builder(self):
        h = highspy.Highs()
        h.silent()

        x, y = h.addVariables(2, lb=-highspy.kHighsInf)
        c = h.addConstrs(-x + y >= 2, x + y >= 0)
        h.minimize(y)
        self.assertAlmostEqual(list(h.val([x, y])), [-1, 1])
        self.assertAlmostEqual(list(h.variableValue([x, y])), [-1, 1])
        self.assertAlmostEqual(list(h.variableValues([x, y])), [-1, 1])
        self.assertAlmostEqual(list(h.allVariableValues()), [-1, 1])

        # -x + y >= 3
        h.changeRowBounds(0, 3, highspy.kHighsInf)
        h.run()
        self.assertAlmostEqual(list(h.val([x, y])), [-1.5, 1.5])

        # make y integer
        h.setInteger(y)
        h.run()
        self.assertAlmostEqual(list(h.val([x, y])), [-1, 2])

        # delete the first constraint and add a new one
        h.removeConstr(c[0], c)
        self.assertEqual(list(map(int, c)), [-1, 0])

        h.addConstr(-x + y >= 0)
        h.run()
        self.assertAlmostEqual(list(h.val([x, y])), [0, 0])

    def test_basics_builder(self):
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)

        x = h.addVariable(lb=-highspy.kHighsInf)
        y = h.addVariable(lb=-highspy.kHighsInf)

        c1 = h.addConstr(-x + y >= 2)
        c2 = h.addConstr(x + y >= 0)

        h.minimize(y)

        self.assertAlmostEqual(h.val(x), -1)
        self.assertAlmostEqual(h.val(y), 1)

        """
        min y
        s.t.
        -x + y >= 3
        x + y >= 0
        """
        h.changeRowBounds(0, 3, highspy.kHighsInf)
        h.run()

        self.assertAlmostEqual(h.val(x), -1.5)
        self.assertAlmostEqual(h.val(y), 1.5)

        sol = h.getSolution()
        self.assertAlmostEqual(h.variableDual(x), sol.col_dual[0])

        self.assertAlmostEqual(h.variableDuals(x), sol.col_dual[0])
        self.assertAlmostEqual(list(h.variableDuals([x, y])), [sol.col_dual[0], sol.col_dual[1]])
        self.assertAlmostEqual(h.variableDuals({"x": x, "y": y}), {"x": sol.col_dual[0], "y": sol.col_dual[1]})
        self.assertAlmostEqual(list(h.variableDuals(np.asarray([x, y]))), [sol.col_dual[0], sol.col_dual[1]])

        self.assertAlmostEqual(h.allVariableDuals(), sol.col_dual)

        c1, c2 = h.getConstrs()
        self.assertAlmostEqual(h.constrValue(c1), sol.row_value[0])

        self.assertAlmostEqual(h.constrValues(c1), sol.row_value[0])
        self.assertAlmostEqual(list(h.constrValues([c1, c2])), sol.row_value)
        self.assertAlmostEqual(list(h.constrValues([c2, c1])), sol.row_value[::-1])  # order matters
        self.assertAlmostEqual(h.constrValues({"c1": c1, "c2": c2}), {"c1": sol.row_value[0], "c2": sol.row_value[1]})
        self.assertAlmostEqual(list(h.constrValues(np.asarray([c1, c2]))), sol.row_value)

        self.assertAlmostEqual(h.allConstrValues(), sol.row_value)

        self.assertAlmostEqual(h.constrDual(c1), sol.row_dual[0])

        self.assertAlmostEqual(h.constrDuals(c1), sol.row_dual[0])
        self.assertAlmostEqual(list(h.constrDuals([c1, c2])), sol.row_dual)
        self.assertAlmostEqual(list(h.constrDuals([c2, c1])), sol.row_dual[::-1])  # order matters
        self.assertAlmostEqual(h.constrDuals({"c1": c1, "c2": c2}), {"c1": sol.row_dual[0], "c2": sol.row_dual[1]})
        self.assertAlmostEqual(list(h.constrDuals(np.asarray([c1, c2]))), sol.row_dual)

        self.assertAlmostEqual(h.allConstrDuals(), sol.row_dual)

        # now make y integer
        h.changeColsIntegrality(1, np.array([1]), np.array([highspy.HighsVarType.kInteger]))
        h.run()
        sol = h.getSolution()
        self.assertAlmostEqual(sol.col_value[0], -1)
        self.assertAlmostEqual(sol.col_value[1], 2)

        """
        now delete the first constraint and add a new one
        
        min y
        s.t.
        x + y >= 0
        -x + y >= 0
        """
        h.removeConstr(c1)

        c1 = h.addConstr(-x + y >= 0)

        h.run()

        self.assertAlmostEqual(h.val(x), 0)
        self.assertAlmostEqual(h.val(y), 0)

        # change the upper bound of x to -5
        h.changeColsBounds(1, np.array([0]), np.array([-highspy.kHighsInf], dtype=np.double), np.array([-5], dtype=np.double))
        h.run()
        self.assertAlmostEqual(h.val(x), -5)
        self.assertAlmostEqual(h.val(y), 5)

        # now maximize
        h.changeRowBounds(1, -highspy.kHighsInf, 0)
        h.changeRowBounds(0, -highspy.kHighsInf, 0)
        h.minimize(-y)

        self.assertAlmostEqual(h.val(x), -5)
        self.assertAlmostEqual(h.val(y), -5)

        self.assertEqual(h.getObjectiveSense()[1], highspy.ObjSense.kMinimize)
        h.maximize(y)
        self.assertEqual(h.getObjectiveSense()[1], highspy.ObjSense.kMaximize)

        self.assertAlmostEqual(h.val(x), -5)
        self.assertAlmostEqual(h.val(y), -5)

        self.assertAlmostEqual(h.getObjectiveValue(), -5)

        h.maximize(y + 1)
        self.assertAlmostEqual(h.getObjectiveOffset()[1], 1)
        self.assertAlmostEqual(h.getObjectiveValue(), -4)

    def test_addVariable(self):
        h = highspy.Highs()
        h.addVariable()
        self.assertEqual(h.numVariables, 1)

        # exception
        self.assertRaises(Exception, lambda: h.addVariable(lb=h.inf))
        self.assertRaises(Exception, lambda: h.addVariables(2, lb=h.inf))

    def test_addConstr(self):
        h = highspy.Highs()
        x = h.addVariable()
        y = h.addVariable()

        c1 = h.addConstr(2 * x + 3 * y <= 10, name="c1")
        self.assertEqual(h.numVariables, 2)
        self.assertEqual(h.numConstrs, 1)
        self.assertEqual(h.getNumNz(), 2)
        self.assertEqual(c1.name, "c1")
        c1.name = "c2"
        self.assertEqual(c1.name, "c2")

        lp = h.getLp()
        self.assertAlmostEqual(lp.row_lower_[0], -highspy.kHighsInf)
        self.assertAlmostEqual(lp.row_upper_[0], 10)

        self.assertEqual(lp.a_matrix_.index_[0], 0)
        self.assertEqual(lp.a_matrix_.index_[1], 1)

        self.assertAlmostEqual(lp.a_matrix_.value_[0], 2)
        self.assertAlmostEqual(lp.a_matrix_.value_[1], 3)

    def test_removeConstr(self):
        h = highspy.Highs()
        x = h.addVariable()
        y = h.addVariable()
        c = h.addConstr(2 * x + 3 * y <= 10)
        self.assertEqual(h.numConstrs, 1)

        h.removeConstr(c)
        self.assertEqual(h.numVariables, 2)
        self.assertEqual(h.numConstrs, 0)
        self.assertRaises(Exception, lambda: c.name("c"))

    def test_val(self):
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)

        x = [h.addVariable(), h.addVariable()]
        h.addConstr(2 * x[0] + 3 * x[1] <= 10)
        h.maximize(x[0])

        self.assertAlmostEqual(h.val(x[0]), 5)

        vals = h.vals(x)
        self.assertAlmostEqual(vals[0], 5)
        self.assertAlmostEqual(vals[1], 0)

        # test linear expr
        self.assertAlmostEqual(h.val(2 * x[0] + 3 * x[1]), 10)
        self.assertAlmostEqual(h.val(1 * x[0] + 3 * x[1]), 5)

        self.assertAlmostEqual(h.val(1 * x[0] + 3 * x[1] <= 4), False)
        self.assertAlmostEqual(h.val(1 * x[0] + 3 * x[1] == 5), True)
        self.assertAlmostEqual(h.val(1 * x[0] + 3 * x[1] >= 6), False)

    def test_var_name(self):
        h = highspy.Highs()

        # name set, but not in the model
        x = h.addVariable(name="x")
        self.assertEqual(x.name, "x")

        # change name before adding to the model
        x.name = "y"
        self.assertEqual(x.name, "y")

        # add to the model
        self.assertEqual(h.numVariables, 1)
        self.assertEqual(h.getLp().col_names_[0], "y")

        # change name after adding to the model
        x.name = "z"
        self.assertEqual(h.getLp().col_names_[0], "z")

        # change name via the model
        h.passColName(0, "a")
        self.assertEqual(h.getLp().col_names_[0], "a")
        self.assertEqual(x.name, "a")
        self.assertEqual(h.variableName(x), "a")
        self.assertRaises(Exception, lambda: h.variableName(1))
        self.assertEqual(h.variableNames([x]), ["a"])
        self.assertEqual(h.variableNames({"key": x}), {"key": "a"})
        self.assertEqual(h.allVariableNames(), ["a"])

        y = h.addVariable()
        h.deleteVariable(y)
        self.assertRaises(Exception, lambda: y.name)

    def test_binary(self):
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)

        x = [h.addBinary(), h.addBinary()]
        h.addConstr(2 * x[0] + 3 * x[1] <= 10)
        h.maximize(x[0])

        lp = h.getLp()
        self.assertAlmostEqual(lp.col_lower_[0], 0)
        self.assertAlmostEqual(lp.col_upper_[0], 1)
        self.assertEqual(lp.integrality_[0], highspy.HighsVarType.kInteger)

        self.assertAlmostEqual(h.val(x[0]), 1)

        vals = h.vals(x)
        self.assertAlmostEqual(vals[0], 1)
        self.assertAlmostEqual(vals[1], 0)

    def test_integer(self):
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)

        x = [h.addIntegral(), h.addVariable()]
        h.addConstr(2 * x[0] + 3 * x[1] <= 10.6)
        h.maximize(x[0] + x[1])

        lp = h.getLp()
        self.assertEqual(lp.integrality_[0], highspy.HighsVarType.kInteger)
        self.assertEqual(lp.integrality_[1], highspy.HighsVarType.kContinuous)

        self.assertAlmostEqual(h.val(x[0]), 5)

        vals = h.vals(x)
        self.assertAlmostEqual(vals[0], 5)
        self.assertAlmostEqual(vals[1], 0.2)

        y = h.addIntegrals(2)
        self.assertEqual(h.getLp().integrality_[2], highspy.HighsVarType.kInteger)
        self.assertEqual(h.getLp().integrality_[3], highspy.HighsVarType.kInteger)

        h.setContinuous(y[0])
        self.assertEqual(h.getLp().integrality_[2], highspy.HighsVarType.kContinuous)
        self.assertEqual(h.getLp().integrality_[3], highspy.HighsVarType.kInteger)

        h.setInteger(y)
        self.assertEqual(h.getLp().integrality_[2], highspy.HighsVarType.kInteger)
        self.assertEqual(h.getLp().integrality_[3], highspy.HighsVarType.kInteger)

        h.setContinuous(y)
        self.assertEqual(h.getLp().integrality_[2], highspy.HighsVarType.kContinuous)
        self.assertEqual(h.getLp().integrality_[3], highspy.HighsVarType.kContinuous)

    def test_objective(self):
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)

        x = [h.addVariable(), h.addVariable()]
        h.addConstr(2 * x[0] + 3 * x[1] <= 10)

        self.assertRaises(Exception, h.maximize, x[0] + x[1] <= 3)
        self.assertRaises(Exception, h.minimize, x[0] + x[1] <= 3)

    def test_constraint_builder(self):
        h = highspy.Highs()
        (x, y) = [h.addVariable(), h.addVariable()]

        # -inf <= 2x + 3y <= inf
        c1 = 2 * x + 3 * y
        self.assertEqualExpr(c1, [x, y], [2, 3])

        # -inf <= 2x + 3y <= 2x
        c1 = 2 * x + 3 * y <= 2 * x
        self.assertEqualExpr(c1, [x, y, x], [2, 3, -2], None, [-highspy.kHighsInf, 0])

        # -inf <= 2x + 3y <= 2x
        c1 = 2 * x >= 2 * x + 3 * y
        self.assertEqualExpr(c1, [x, y, x], [2, 3, -2], None, [-highspy.kHighsInf, 0])

        # failure add constraint without inequality
        self.assertRaises(Exception, lambda: h.addConstr(x + 3 * y))
        self.assertRaises(Exception, lambda: h.addConstrs(x == 1, x + 3 * y))
        self.assertRaises(Exception, lambda: h.addConstrs({"a": x == 1, "b": x + 3 * y}))

        # ensure model is rolled back on error
        self.assertEqual(h.numConstrs, 0)
        self.assertRaises(Exception, lambda: h.addConstrs([x == 1] * 100 + [x + 3 * y]))
        self.assertEqual(h.numConstrs, 0)

        # failure - bounds already set
        self.assertRaises(Exception, lambda: 1 <= (4 <= 2 * x + 3 * y))
        self.assertRaises(Exception, lambda: 2 >= (4 >= 2 * x + 3 * y))
        self.assertRaises(Exception, lambda: (2 * x + 3 * y <= 2) <= 4)
        self.assertRaises(Exception, lambda: (1 <= 2 * x + 3 * y) <= 5)
        self.assertRaises(Exception, lambda: 1 <= (2 * x + 3 * y <= 5))
        self.assertRaises(Exception, lambda: 1 <= (5 >= 2 * x + 3 * y))
        self.assertRaises(Exception, lambda: (5 >= 2 * x + 3 * y) >= 1)
        self.assertRaises(Exception, lambda: 5 >= (2 * x + 3 * y >= 1))

        c1 = 2 * x + 3 * y <= (2 <= 4)  # 2x + 3y <= float(True)
        self.assertEqualExpr(c1, [x, y], [2, 3], None, [-highspy.kHighsInf, 1])

        # failure, non-linear terms
        self.assertRaises(Exception, lambda: 2 * x * 3 * y)

        # failure, order matters when having variables on both sides of inequality
        self.assertRaises(Exception, lambda: (4 * x <= 2 * x + 3 * y) <= 5)
        self.assertRaises(Exception, lambda: 4 * x <= (2 * x + 3 * y <= 5))
        self.assertRaises(Exception, lambda: (4 * x <= 2 * x + 3 * y) <= 5 * x)

        c1 = 5 * x <= 2 * x + 3 * y <= 5 * x
        self.assertEqualExpr(c1, [x, y, x], [2, 3, -5], None, [0, 0])

        # test various combinations with different inequalities
        self.assertRaises(Exception, lambda: (2 * x + 3 * y == 3 * y) == 3)
        self.assertRaises(Exception, lambda: 2 * x + 3 * y == (3 * y == 3))
        self.assertRaises(Exception, lambda: 2 * x + 3 * y == (3 * y <= 3))
        self.assertRaises(Exception, lambda: 2 * x + 3 * y == (3 * y >= 3))

        c1 = 2 * x + 3 * y == x
        self.assertEqualExpr(c1, [x, y, x], [2, 3, -1], None, [0, 0])

        c1 = 2 * x + 3 * y == 5
        self.assertEqualExpr(c1, [x, y], [2, 3], None, [5, 5])

        c1 = 5 == 2 * x + 3 * y
        self.assertEqualExpr(c1, [x, y], [2, 3], None, [5, 5])

        # 2*x + 3*y == 4.5
        c1 = 2 * x + 3 * y + 0.5 == 5
        self.assertEqualExpr(c1, [x, y], [2, 3], None, [4.5, 4.5])
        h.addConstr(c1)
        self.assertAlmostEqual((h.getLp().row_lower_[0], h.getLp().row_upper_[0]), (4.5, 4.5))

    def test_add_multiple_variables(self):
        # test basic functionality
        h = highspy.Highs()
        x = h.addVariables(2)
        self.assertEqual(h.numVariables, 2)

        # test multiple dimensions
        h = highspy.Highs()
        x = h.addVariables(2, 3, 4, out_array=False)
        self.assertEqual(h.numVariables, 2 * 3 * 4)
        self.assertEqual(isinstance(x, dict), True)

        # test multiple dimensions array
        h = highspy.Highs()
        x = h.addVariables(2, 3, 4)
        self.assertEqual(h.numVariables, 2 * 3 * 4)
        self.assertEqual(x.shape, (2, 3, 4))

        # test binary variables with objective and names
        h = highspy.Highs()
        x = h.addBinaries(20, obj=range(20), name_prefix="t_")
        self.assertEqual(h.numVariables, 20)
        self.assertEqual(h.getLp().col_names_[0], "t_0")
        self.assertEqual(h.getLp().col_names_[19], "t_19")
        self.assertEqual(h.getLp().col_cost_[0], 0)
        self.assertEqual(h.getLp().col_cost_[19], 19)

        # test prefix item with indices, not variable offset
        h = highspy.Highs()
        x = h.addVariables("a", "b", "c", name_prefix="t_")  # ('a', 'b', 'c')
        self.assertEqual(h.numVariables, 1)
        self.assertEqual(x["a", "b", "c"].name, "t_('a', 'b', 'c')")

        # Testing different ways of adding variables
        # Some are unlikely to be used, but this is expected behaviour
        N = 0
        h = highspy.Highs()
        x1 = h.addVariables(("a", "b", "c"), obj={"b": 20, "c": 10, "a": 50})
        N += 3
        self.assertEqual(h.numVariables, N)
        lp = h.getLp()
        self.assertEqual(lp.col_cost_[0], 50)
        self.assertEqual(lp.col_cost_[1], 20)
        self.assertEqual(lp.col_cost_[2], 10)

        x2 = h.addVariables(["a", "b", "c", "d"])  # 'a', 'b', 'c', 'd'
        N += 4
        self.assertEqual(h.numVariables, N)

        x3 = h.addVariables("abc")  # 'a', 'b', 'c'
        N += 3
        self.assertEqual(h.numVariables, N)

        x4 = h.addVariables(("ab", "b", "c", "d"))  # 'ab', 'b', 'c', 'd'
        N += 4
        self.assertEqual(h.numVariables, N)

        x5 = h.addVariables("ab", "b", "c", "d")  # ('a', 'b', 'c', 'd'), ('b', 'b', 'c', 'd')
        N += 2
        self.assertEqual(h.numVariables, N)
        self.assertTrue(("a", "b", "c", "d") in x5.keys())
        self.assertTrue(("b", "b", "c", "d") in x5.keys())

        x6 = h.addVariables(5, "a", 2, "b", "c")  # range(5), 'a', range(2), 'b', 'c'
        N += 5 * 2
        self.assertEqual(h.numVariables, N)

        x7 = h.addVariables([5, "a", 2, "b", "c"])  # 5, 'a', 2, 'b', 'c'
        N += 5
        self.assertEqual(h.numVariables, N)

        x8 = h.addVariables([(20, 1), (1, 2), (2, 6)], ub=[3, 2, 1], name_prefix="t")  # (20, 1), (1,2), (2,6)
        N += 3
        self.assertEqual(h.numVariables, N)

        x9 = h.addBinaries((20, 1), (1, 2), (2, 6))  # product((20, 1), (1,2), (2,6))) = (20, 1, 2), ..., (1, 2, 6)
        N += 8
        self.assertEqual(h.numVariables, N)

    def test_add_single_constraints(self):
        h = highspy.Highs()
        (x, y) = h.addVariables(2)
        added_constraints = h.addConstrs([2 * x + 3 * y <= 5])
        self.assertEqual(len(added_constraints), 1)
        self.assertEqual(h.numConstrs, 1)

    def test_add_multiple_constraints(self):
        # test manual constraints
        h = highspy.Highs()
        (x, y, z) = h.addVariables(3)
        added_constraints = h.addConstrs([x + y <= 5, 2 * x + 3 * y >= 10, x - z == 2])
        self.assertEqual(len(added_constraints), 3)
        self.assertEqual(h.numConstrs, 3)

        # test list comprehension constraints
        h = highspy.Highs()
        x = h.addVariables(5)
        h.addConstr(sum(x) == 1)
        self.assertEqual(h.numConstrs, 1)

        h.addConstrs(x[i] + x[j] <= 1 for i in range(5) for j in range(5))
        self.assertEqual(h.numConstrs, 26)

        # test names and sequence types (list, tuple, nothing)
        h = highspy.Highs()
        (x1, x2, x3) = h.addVariables(3)

        h.addConstrs((x2 - x1 >= 2), name_prefix="a")
        self.assertEqual(h.numConstrs, 1)
        h.addConstrs((x2 - x1 >= 2))
        self.assertEqual(h.numConstrs, 2)

        h.addConstrs(x2 - x1 >= 2, name_prefix="b")
        self.assertEqual(h.numConstrs, 3)
        h.addConstrs(x2 - x1 >= 2)
        self.assertEqual(h.numConstrs, 4)

        h.addConstrs([x2 - x1 >= 2], name_prefix="c")
        self.assertEqual(h.numConstrs, 5)
        h.addConstrs([x2 - x1 >= 2])
        self.assertEqual(h.numConstrs, 6)

    # r/w basis tests below works on unix but not windows?
    def test_write_basis_before_running(self):
        if platform == "linux" or platform == "darwin":
            h = self.get_basic_model()
            with tempfile.NamedTemporaryFile() as f:
                h.writeBasis(f.name)
                contents = f.read()
                self.assertEqual(contents, b"HiGHS v1\nNone\n")

    def test_write_basis_after_running(self):
        if platform == "linux" or platform == "darwin":
            h = self.get_basic_model()
            h.run()
            with tempfile.NamedTemporaryFile() as f:
                h.writeBasis(f.name)
                contents = f.read()
                self.assertEqual(contents, b"HiGHS v1\nValid\n# Columns 2\n1 1 \n# Rows 2\n0 0 \n")

    def test_read_basis(self):
        if platform == "linux" or platform == "darwin":
            # Read basis from one run model into an unrun model
            expected_status_before = highspy.HighsBasisStatus.kLower
            expected_status_after = highspy.HighsBasisStatus.kBasic

            h1 = self.get_basic_model()
            self.assertEqual(h1.getBasis().col_status[0], expected_status_before)
            h1.run()
            self.assertEqual(h1.getBasis().col_status[0], expected_status_after)

            h2 = self.get_basic_model()
            self.assertEqual(h2.getBasis().col_status[0], expected_status_before)

            with tempfile.NamedTemporaryFile() as f:
                h1.writeBasis(f.name)
                h2.readBasis(f.name)
                self.assertEqual(h2.getBasis().col_status[0], expected_status_after)

    def test_solve(self):
        """Test the solve method to ensure it runs the solver."""
        h = highspy.Highs()
        h.silent()
        x = h.addBinary(obj=1)
        h.setMaximize()
        h.solve()
        self.assertEqual(h.getSolution().col_value[0], 1)
        self.assertEqual(h.getObjectiveSense()[1], highspy.highs.ObjSense.kMaximize)

        h.setMinimize()
        h.optimize()
        self.assertEqual(h.getSolution().col_value[0], 0)
        self.assertEqual(h.getObjectiveSense()[1], highspy.highs.ObjSense.kMinimize)

    def test_minimize(self):
        """Test the minimize method with and without an objective."""
        h = highspy.Highs()
        h.silent()
        x, y, z = h.addBinaries(3, obj=-1)
        h.minimize()
        self.assertEqual(list(h.val([x, y, z])), [1, 1, 1])

        h.minimize(-x)
        self.assertEqual(h.val(x), 1)

        h.minimize(x)
        self.assertEqual(h.val(x), 0)

        h.minimize(x - y)
        self.assertEqual(h.val(x), 0)
        self.assertEqual(h.val(y), 1)

        self.assertRaises(Exception, h.minimize, x - y <= 5)
        self.assertRaises(Exception, h.minimize, 0 <= x - y <= 5)
        self.assertRaises(Exception, h.minimize, 4 <= x - y)
        self.assertRaises(Exception, h.minimize, 4)

    def test_maximize(self):
        """Test the maximize method with and without an objective."""
        h = highspy.Highs()
        h.silent()
        x, y, z = h.addBinaries(3, obj=1)
        h.maximize()
        self.assertEqual(list(h.val([x, y, z])), [1, 1, 1])

        h.maximize(x)
        self.assertEqual(h.val(x), 1)

        h.maximize(-x)
        self.assertEqual(h.val(x), 0)

        h.maximize(y - x)
        self.assertEqual(h.val(x), 0)
        self.assertEqual(h.val(y), 1)

        self.assertRaises(Exception, h.maximize, x - y <= 5)
        self.assertRaises(Exception, h.maximize, 0 <= x - y <= 5)
        self.assertRaises(Exception, h.maximize, 4 <= x - y)
        self.assertRaises(Exception, h.maximize, 4)

    def test_get_expr(self):
        h = self.get_basic_model()

        expr = h.getExpr(0)  # -x + y >= 2
        self.assertEqualExpr(expr, [0, 1], [-1, 1], None, [2, highspy.kHighsInf])

        expr = h.getExpr(1)  # x + y >= 0
        self.assertEqualExpr(expr, [0, 1], [1, 1], None, [0, highspy.kHighsInf])

        c = h.getConstrs()
        self.assertEqualExpr(c[0].expr(), [0, 1], [-1, 1], None, [2, highspy.kHighsInf])
        self.assertEqualExpr(c[1].expr(), [0, 1], [1, 1], None, [0, highspy.kHighsInf])

        self.assertRaises(Exception, lambda: h.getExpr(2))

    def test_add_variables(self):
        """Test adding multiple variables to the model."""
        h = highspy.Highs()

        var = h.addVariables()
        self.assertEqual(var, None)

        keys = ["a", "b", "c"]
        v = h.addVariables(keys)
        self.assertTrue(isinstance(v, dict))
        self.assertEqual(int(v["a"]), 0)

        h.addConstr(v["a"] + v["b"] + v["c"] == 1)
        h.maximize(v["a"])
        self.assertEqual(h.val(v), {"a": 1.0, "b": 0, "c": 0})
        self.assertEqual(list(map(int, h.getVariables())), [0, 1, 2])

        # provided len(parameters) != len(indices)
        self.assertRaises(Exception, lambda: h.addVariables(5, type=[highspy.HighsVarType.kContinuous, highspy.HighsVarType.kInteger]))
        self.assertRaises(Exception, lambda: h.addVariables(5, name=["a", "b"]))

        # some parameters not valid
        self.assertRaises(Exception, lambda: h.addVariables(5, type=[highspy.HighsVarType.kContinuous, None]))
        self.assertRaises(Exception, lambda: h.addVariables(5, name=["a", None]))
        self.assertRaises(Exception, lambda: h.addVariables(5, out_array=["a"]))
        self.assertRaises(Exception, lambda: h.addVariables(5, name_prefix=["a"]))

        # correct usage
        y = h.addVariables(5, type=[highspy.HighsVarType.kContinuous] * 5)
        self.assertEqual(list(map(int, y)), [3, 4, 5, 6, 7])

        y2 = h.addVariables(
            ["a", "b", "c"], type={"a": highspy.HighsVarType.kContinuous, "b": highspy.HighsVarType.kInteger, "c": highspy.HighsVarType.kContinuous}
        )
        self.assertEqual({k: int(v) for k, v in y2.items()}, {"a": 8, "b": 9, "c": 10})

        y3 = h.addVariables(2, name=["a", "b"])
        self.assertEqual(list(map(int, y3)), [11, 12])

    def test_delete_variable(self):
        h = highspy.Highs()

        keys = ["a", "b", "c"]
        D = h.addVariables(keys)
        X = h.addVariables(5)

        self.assertEqual({k: int(v) for k, v in D.items()}, {"a": 0, "b": 1, "c": 2})
        self.assertEqual(list(map(int, X)), [3, 4, 5, 6, 7])

        # delete variable and update collections
        h.deleteVariable(D["b"], D, X)

        self.assertEqual(h.numVariables, 7)
        self.assertEqual({k: int(v) for k, v in D.items()}, {"a": 0, "b": -1, "c": 1})
        self.assertEqual(list(map(int, X)), [2, 3, 4, 5, 6])

        # delete variable and update collections
        h.deleteVariable(X[3], D, X)

        self.assertEqual(h.numVariables, 6)
        self.assertEqual({k: int(v) for k, v in D.items()}, {"a": 0, "b": -1, "c": 1})
        self.assertEqual(list(map(int, X)), [2, 3, 4, -1, 5])

        # delete variable and update collections
        h.deleteVariable(X[2], D, *X)

        self.assertEqual(h.numVariables, 5)
        self.assertEqual({k: int(v) for k, v in D.items()}, {"a": 0, "b": -1, "c": 1})
        self.assertEqual(list(map(int, X)), [2, 3, -1, -1, 4])

    def test_remove_constraint(self):
        h = highspy.Highs()

        keys = ["a", "b", "c"]
        D = h.addVariables(keys)
        X = h.addVariables(5)

        c1 = h.addConstr(D["a"] + D["b"] + D["c"] == 1)
        cX = h.addConstrs((x >= 1 for x in X))
        c2 = h.addConstr(qsum(X) == 1)
        c3 = h.addConstr(qsum(D.values()) == 1)
        hash_test = {c1: "c1", c2: "c2", c3: "c3"}

        self.assertEqual(h.numConstrs, 8)
        self.assertEqual([int(c1)] + list(map(int, cX)) + [int(c2), int(c3)], [0, 1, 2, 3, 4, 5, 6, 7])

        # delete constr and update collections
        h.removeConstr(cX[2], c1, cX, c2, c3)

        self.assertEqual(h.numConstrs, 7)
        self.assertEqual([int(c1)] + list(map(int, cX)) + [int(c2), int(c3)], [0, 1, 2, -1, 3, 4, 5, 6])

        # delete variable and update collections
        h.removeConstr(c1, c1, cX)

        self.assertEqual(h.numConstrs, 6)
        self.assertEqual([int(c1)] + list(map(int, cX)) + [int(c2), int(c3)], [-1, 0, 1, -1, 2, 3, 5, 6])

        # add dictionary of constraints with highs_var keys
        cM = h.addConstrs({x: x >= 1 for x in X})
        self.assertEqual(h.numConstrs, 11)
        self.assertEqual([int(c1)] + list(map(int, cX)) + [int(c2), int(c3)], [-1, 0, 1, -1, 2, 3, 5, 6])
        self.assertEqual({int(k): int(c) for k, c in cM.items()}, {3: 6, 4: 7, 5: 8, 6: 9, 7: 10})

        # add dictionary of constraints with string keys
        cD = h.addConstrs({k: x >= 1 for k, x in D.items()})
        self.assertEqual(h.numConstrs, 14)
        self.assertEqual([int(c1)] + list(map(int, cX)) + [int(c2), int(c3)], [-1, 0, 1, -1, 2, 3, 5, 6])
        self.assertEqual({int(k): int(c) for k, c in cM.items()}, {3: 6, 4: 7, 5: 8, 6: 9, 7: 10})
        self.assertEqual({k: int(c) for k, c in cD.items()}, {"a": 11, "b": 12, "c": 13})

        # delete constraint from dictionary
        h.removeConstr(cM[X[2]], c1, cX, c2, c3, cM, cD)
        self.assertEqual(h.numConstrs, 13)
        self.assertEqual([int(c1)] + list(map(int, cX)) + [int(c2), int(c3)], [-1, 0, 1, -1, 2, 3, 5, 6])  # unchanged
        self.assertEqual({int(k): int(c) for k, c in cM.items()}, {3: 6, 4: 7, 5: -1, 6: 8, 7: 9})
        self.assertEqual({k: int(c) for k, c in cD.items()}, {"a": 10, "b": 11, "c": 12})

        # delete constraint from dictionary
        h.removeConstr(cD["b"], c1, cX, c2, c3, cM, cD)
        self.assertEqual(h.numConstrs, 12)
        self.assertEqual([int(c1)] + list(map(int, cX)) + [int(c2), int(c3)], [-1, 0, 1, -1, 2, 3, 5, 6])  # unchanged
        self.assertEqual({int(k): int(c) for k, c in cM.items()}, {3: 6, 4: 7, 5: -1, 6: 8, 7: 9})  # unchanged
        self.assertEqual({k: int(c) for k, c in cD.items()}, {"a": 10, "b": -1, "c": 11})

        # delete non-existent constraint
        self.assertRaises(Exception, lambda: h.removeConstr(cD["b"], c1, cX, c2, c3, cM, cD))

    def test_qsum(self):
        """Test summation."""
        h = highspy.Highs()
        X = h.addVariables(10)

        # qsum
        expr = qsum(X)
        self.assertEqualExpr(expr, X, [1] * 10)

        expr = qsum((2 * x for x in X))
        self.assertEqualExpr(expr, X, [2] * 10)

        expr = qsum((qsum(X) for k in range(11))).simplify()
        self.assertEqualExpr(expr, X, [11] * 10)

        # sum
        expr = sum(X)
        self.assertEqualExpr(expr, X, [1] * 10, 0)

        expr = sum((2 * x for x in X))
        self.assertEqualExpr(expr, X, [2] * 10, 0)

        # init sum with empty expression
        expr = sum((2 * x for x in X), highs_linear_expression())
        self.assertEqualExpr(expr, X, [2] * 10)

        # manual sum
        expr = h.expr()
        for x in X:
            expr += x
        self.assertEqualExpr(expr, X, [1] * 10)

    def test_user_interrupts(self):
        N = 8
        h = highspy.Highs()
        h.silent()

        x = h.addBinaries(N, N)
        y = np.fliplr(x)

        h.addConstrs(x.sum(axis = 0) == 1)  # each col has exactly one queen
        h.addConstrs(x.sum(axis = 1) == 1)  # each row has exactly one queen

        h.addConstrs(x.diagonal(k).sum() <= 1 for k in range(-N + 1, N))  # each diagonal has at most one queen
        h.addConstrs(y.diagonal(k).sum() <= 1 for k in range(-N + 1, N))  # each 'reverse' diagonal has at most one queen

        h.HandleUserInterrupt = True
        t = h.startSolve()
        self.assertRaises(Exception, lambda: h.startSolve())
        h.cancelSolve()
        h.wait()

        h = self.get_basic_model()
        h.HandleKeyboardInterrupt = True
        self.assertEqual(h.HandleKeyboardInterrupt, True)
        self.assertEqual(h.HandleUserInterrupt, True)

        h.solve()
        h.minimize()
        h.maximize()
        h.minimize()

        h.joinSolve(h.startSolve())
        h.joinSolve(h.startSolve(), 0)

        # replace wait function with Ctrl+C signal
        highspy.highs.Highs.wait = lambda self, t: signal.raise_signal(signal.SIGINT)
        h.HandleKeyboardInterrupt = False

        self.assertEqual(h.HandleKeyboardInterrupt, False)
        self.assertEqual(h.HandleUserInterrupt, False)

        h.joinSolve(h.startSolve(), 0)
        h.startSolve()
        h.joinSolve(None, 0)

        with self.assertRaises(SystemExit):
            h.startSolve()
            h.joinSolve(None, 5)
            unittest.main(exit=False)

    def test_callbacks(self):
        N = 8
        h = highspy.Highs()
        h.silent(False)

        x = h.addBinaries(N, N)
        y = np.fliplr(x)

        h.addConstrs(h.qsum(x[i, :]) == 1 for i in range(N))  # each row has exactly one queen
        h.addConstrs(h.qsum(x[:, j]) == 1 for j in range(N))  # each col has exactly one queen

        h.addConstrs(h.qsum(x.diagonal(k)) <= 1 for k in range(-N + 1, N))  # each diagonal has at most one queen
        h.addConstrs(h.qsum(y.diagonal(k)) <= 1 for k in range(-N + 1, N))  # each 'reverse' diagonal has at most one queen

        do_nothing = lambda e: None

        check_callback = []
        chk_callback = lambda e: check_callback.append(e)

        check_solution = []
        chk_solution = lambda e: check_solution.append(e.val(x))

        # check callback is added and called
        h.cbLogging += chk_callback
        h.cbMipSolution += chk_solution
        self.assertEqual(len(h.cbLogging.callbacks), 1)
        self.assertEqual(len(h.cbMipSolution.callbacks), 1)
        h.solve()
        self.assertNotEqual(len(check_callback), 0)
        self.assertNotEqual(len(check_solution), 0)
        check_callback.clear()
        check_solution.clear()

        # check callback is removed and not called
        h.cbLogging -= chk_callback
        h.cbMipSolution -= chk_solution
        self.assertEqual(len(h.cbLogging.callbacks), 0)
        self.assertEqual(len(h.cbMipSolution.callbacks), 0)
        h.solve()
        self.assertEqual(len(check_callback), 0)
        self.assertEqual(len(check_solution), 0)

        h.disableCallbacks()
        self.assertRaises(Exception, lambda: h.cbLogging.subscribe(do_nothing))
        h.enableCallbacks()

        h.cbLogging += chk_callback
        h.cbSimplexInterrupt += do_nothing
        h.cbIpmInterrupt += do_nothing
        h.cbMipSolution += do_nothing
        h.cbMipImprovingSolution += do_nothing
        h.cbMipLogging += do_nothing
        h.cbMipInterrupt += do_nothing
        h.cbMipGetCutPool += do_nothing
        h.cbMipDefineLazyConstraints += do_nothing

        self.assertEqual(len(h.cbLogging.callbacks), 1)

        # check callback is disabled and not called
        h.disableCallbacks()
        self.assertEqual(len(h.cbLogging.callbacks), 1)
        h.solve()
        self.assertEqual(len(check_callback), 0)
        check_callback.clear()

        # check callback is enabled and called
        h.enableCallbacks()
        self.assertEqual(len(h.cbLogging.callbacks), 1)
        h.solve()
        self.assertNotEqual(len(check_callback), 0)
        check_callback.clear()

        h.cbLogging -= chk_callback
        h.cbSimplexInterrupt -= do_nothing
        h.cbIpmInterrupt -= do_nothing
        h.cbMipSolution -= do_nothing
        h.cbMipImprovingSolution -= do_nothing
        h.cbMipLogging -= do_nothing
        h.cbMipInterrupt -= do_nothing
        h.cbMipGetCutPool -= do_nothing
        h.cbMipDefineLazyConstraints -= do_nothing

        self.assertEqual(len(h.cbLogging.callbacks), 0)
        h.cbLogging.subscribe(do_nothing)
        self.assertEqual(len(h.cbLogging.callbacks), 1)
        h.cbLogging.unsubscribe(do_nothing)
        self.assertEqual(len(h.cbLogging.callbacks), 0)
        h.cbLogging.unsubscribe(do_nothing)  # does not throw error if not exists
        self.assertEqual(len(h.cbLogging.callbacks), 0)

        h.cbLogging.subscribe(do_nothing, "test")
        self.assertEqual(len(h.cbLogging.callbacks), 1)
        h.cbLogging.unsubscribe_by_data("test")
        self.assertEqual(len(h.cbLogging.callbacks), 0)

        h.cbLogging += do_nothing
        h.cbLogging += do_nothing
        h.cbLogging += do_nothing
        self.assertEqual(len(h.cbLogging.callbacks), 3)
        h.cbLogging.clear()
        self.assertEqual(len(h.cbLogging.callbacks), 0)

        self.assertRaises(Exception, lambda: h.__setattr__("cbLogging", None))
        self.assertRaises(Exception, lambda: h.__setattr__("cbSimplexInterrupt", None))
        self.assertRaises(Exception, lambda: h.__setattr__("cbIpmInterrupt", None))
        self.assertRaises(Exception, lambda: h.__setattr__("cbMipSolution", None))
        self.assertRaises(Exception, lambda: h.__setattr__("cbMipImprovingSolution", None))
        self.assertRaises(Exception, lambda: h.__setattr__("cbMipLogging", None))
        self.assertRaises(Exception, lambda: h.__setattr__("cbMipInterrupt", None))
        self.assertRaises(Exception, lambda: h.__setattr__("cbMipGetCutPool", None))
        self.assertRaises(Exception, lambda: h.__setattr__("cbMipDefineLazyConstraints", None))


class TestHighsLinearExpressionPy(unittest.TestCase):
    def setUp(self):
        self.h = highspy.Highs()
        self.h.silent()

        self.x = self.h.addVariables(10)

    def assertEqualExpr(self, expr, idxs, vals, constant=None, bounds=None):
        self.assertEqual(list(map(int, expr.idxs)), list(map(int, idxs)), "variable index")
        self.assertEqual(expr.vals, vals, "variable values")
        self.assertEqual(expr.constant, constant, "constant")
        self.assertEqual(expr.bounds, (bounds[0], bounds[1]) if bounds is not None else None, "bounds")

    def test_init_empty(self):
        # Test initialization with no arguments
        expr = highspy.highs.highs_linear_expression()
        self.assertEqualExpr(expr, [], [])

        expr = self.h.expr()
        self.assertEqualExpr(expr, [], [])
        self.assertRaises(Exception, lambda: self.h.expr([]))
        self.assertRaises(Exception, lambda: self.h.expr(self.h))

    def test_init_var(self):
        # Test initialization with a highs_var
        expr = highspy.highs.highs_linear_expression(self.x[0])
        self.assertEqualExpr(expr, [self.x[0]], [1.0])

        expr = 1.0 * self.x[0]
        self.assertEqualExpr(expr, [self.x[0]], [1.0])

        expr = self.h.expr(self.x[0])
        self.assertEqualExpr(expr, [self.x[0]], [1.0])

    def test_init_const(self):
        # Test initialization with a constant
        expr = highspy.highs.highs_linear_expression(5)
        self.assertEqualExpr(expr, [], [], 5)

        expr = self.h.expr(5)
        self.assertEqualExpr(expr, [], [], 5)

    def test_mutable(self):
        x, y, z = self.x[0:3]

        expr = x + 2 * y
        expr2 = expr  # reference to expr
        expr3 = expr.copy()  # copy of expr
        self.assertEqualExpr(expr, [x, y], [1, 2])
        self.assertEqualExpr(expr2, [x, y], [1, 2])
        self.assertEqualExpr(expr3, [x, y], [1, 2])

        expr += z
        self.assertEqualExpr(expr, [x, y, z], [1, 2, 1])
        self.assertEqualExpr(expr2, [x, y, z], [1, 2, 1])
        self.assertEqualExpr(expr3, [x, y], [1, 2])

        expr *= 2
        self.assertEqualExpr(expr, [x, y, z], [2, 4, 2])
        self.assertEqualExpr(expr2, [x, y, z], [2, 4, 2])
        self.assertEqualExpr(expr3, [x, y], [1, 2])

        # test simplify
        expr += x + y + z
        expr += x + y + z
        expr += x + y + z
        self.assertEqualExpr(expr, [x, y, z] + [x, y, z] * 3, [2, 4, 2] + [1, 1, 1] * 3)

        expr4 = expr.simplify()
        self.assertEqualExpr(expr, [x, y, z] + [x, y, z] * 3, [2, 4, 2] + [1, 1, 1] * 3)
        self.assertEqualExpr(expr4, [x, y, z], [5, 7, 5])

        # test edge cases
        e1 = x + z <= 3
        e2 = 2 * y

        # addition
        self.assertRaises(Exception, lambda: e1 + (e2 + 3))  # cannot add if one has bounds and the other has constant
        self.assertRaises(Exception, lambda: e1 + 5)  # cannot add constant to expr with bounds
        self.assertRaises(Exception, lambda: e1 + [])  # unknown type

        expr = e1 + (1 <= (e2 + 4) <= 2)
        self.assertEqualExpr(expr, [x, z, y], [1, 1, 2], None, [-self.h.inf, 1])

        expr = e2.copy()
        expr += e1
        self.assertEqualExpr(expr, [y, x, z], [2, 1, 1], None, [-self.h.inf, 3])

        # subtract
        self.assertRaises(Exception, lambda: e1 - (e2 + 3))  # cannot add if one has bounds and the other has constant
        self.assertRaises(Exception, lambda: e1 - 5)  # cannot add constant to expr with bounds
        self.assertRaises(Exception, lambda: e1 - [])  # unknown type

        expr = e1 - (1 <= (e2 + 4) <= 2)  # (-inf <= x + z <= 3) + (2 <= -2y <= 3)
        self.assertEqualExpr(expr, [x, z, y], [1, 1, -2], None, [-self.h.inf, 6])

        expr = e2.copy()
        expr -= e1
        self.assertEqualExpr(expr, [y, x, z], [2, -1, -1], None, [-3, self.h.inf])

    def test_immutable(self):
        x, y, z = self.x[0:3]

        expr = (1 * x - 2 * y).simplify()
        self.assertEqualExpr(expr, [x, y], [1, -2])

        expr2 = expr + z
        self.assertEqualExpr(expr, [x, y], [1, -2])
        self.assertEqualExpr(expr2, [x, y, z], [1, -2, 1])

        expr2 = expr + 5
        self.assertEqualExpr(expr, [x, y], [1, -2])
        self.assertEqualExpr(expr2, [x, y], [1, -2], 5)

        expr2 = expr <= 3
        self.assertEqualExpr(expr, [x, y], [1, -2])
        self.assertEqualExpr(expr2, [x, y], [1, -2], None, [-self.h.inf, 3])

        expr2 = 0 <= expr <= 3
        self.assertEqualExpr(expr, [x, y], [1, -2])
        self.assertEqualExpr(expr2, [x, y], [1, -2], None, [0, 3])

        expr2 = 0 <= expr <= 3
        expr2 += z
        self.assertEqualExpr(expr, [x, y], [1, -2])
        self.assertEqualExpr(expr2, [x, y, z], [1, -2, 1], None, [0, 3])

    def test_negation(self):
        # Test negation of a highs_var
        expr = -self.x[0]
        self.assertEqualExpr(expr, [self.x[0]], [-1])

        # Test negation of a highs_linear_expression
        x, y = self.x[0:2]
        expr = x - 2 * y
        self.assertEqualExpr(expr, [x, y], [1, -2])
        negr = -expr
        self.assertEqualExpr(negr, [x, y], [-1, 2])

        # Test negation of a highs_linear_expression
        expr = -self.h.qsum(self.x)
        self.assertEqualExpr(expr, list(map(int, self.x)), [-1] * len(self.x))

    def test_equality(self):
        x, y = self.x[0:2]

        expr = x == y
        self.assertEqualExpr(expr, [x, y], [1, -1], None, [0, 0])

        expr = x + y == [1, 2]
        self.assertEqualExpr(expr, [x, y], [1, 1], None, [1, 2])
        self.assertRaises(Exception, lambda: x == [x, y])
        self.assertRaises(Exception, lambda: x == self.h)
        self.assertRaises(Exception, lambda: x != self.h)
        self.assertRaises(Exception, lambda: x + y != y)
        self.assertRaises(Exception, lambda: x + 5 != self.h)

    def test_le_inequality(self):
        x, y = self.x[0:2]

        # Test inequality of two highs_linear_expressions
        expr = x <= y
        self.assertEqualExpr(expr, [x, y], [1, -1], None, [-self.h.inf, 0])
        self.assertRaises(Exception, lambda: x <= self.h)

    def test_ge_inequality(self):
        x, y = self.x[0:2]

        # Test inequality of two highs_linear_expressions
        expr = x >= y  # y - x <= 0
        self.assertEqualExpr(expr, [y, x], [1, -1], None, [-self.h.inf, 0])
        self.assertRaises(Exception, lambda: x >= self.h)
        self.assertRaises(Exception, lambda: x + 4 >= self.h)
        self.assertRaises(Exception, lambda: x + 4 >= (y <= 2))

    def test_chain_inequality(self):
        x, y, z = self.x[0:3]

        # test basic chain inequality
        expr = 2 <= x + y <= 6
        self.assertEqualExpr(expr, [x, y], [1, 1], None, [2, 6])

        expr = 2 <= (x + y) <= 6
        self.assertEqualExpr(expr, [x, y], [1, 1], None, [2, 6])

        expr = 2 <= x <= 6
        self.assertEqualExpr(expr, [x], [1], None, [2, 6])

        # advanced chain use cases
        expr = y <= 6 + x <= y  # -6 <= x - y <= -6
        self.assertEqualExpr(expr, [x, y], [1, -1], None, [-6, -6])

        expr = x <= 6 <= x  # x == 6
        self.assertEqualExpr(expr, [x], [1], None, [6, 6])

        expr = x + y <= 6 <= x + y  # 6 <= x + y <= 6
        self.assertEqualExpr(expr, [x, y], [1, 1], None, [6, 6])

        expr = x + y + 1 <= 6 <= x + y + 2  # 4 <= x + y <= 5
        self.assertEqualExpr(expr, [x, y], [1, 1], None, [4, 5])

        # test chain ordering with a constant
        t1 = x + y
        t2 = x - y
        self.assertEqualExpr(t1, [x, y], [1, 1])
        self.assertEqualExpr(t2, [x, y], [1, -1])

        expr = (t1 + t2).simplify()  # 2x + 0y
        self.assertEqualExpr(expr, [x, y], [2, 0])

        t3 = 2 <= expr <= 4  # 2 <= 2x + 0y <= 4
        self.assertEqualExpr(t3, [x, y], [2, 0], None, [2, 4])

        vx, vl = [x, y, x, y], [1, 1, 1, -1]

        t3 = 2 <= (t1 + t2) <= 4  # 2 <= 2x + 0y <= 4
        self.assertEqualExpr(t3, vx, vl, None, [2, 4])

        t3 = 2 <= t1 + t2 <= 4  # 2 <= 2x + 0y <= 4
        self.assertEqualExpr(t3, vx, vl, None, [2, 4])

        t3 = (t1 + t2) <= 4  # -inf <= 2x + 0y <= 4
        self.assertEqualExpr(t3, vx, vl, None, [-self.h.inf, 4])

        t3 = t1 + t2 <= 4  # -inf <= 2x + 0y <= 4
        self.assertEqualExpr(t3, vx, vl, None, [-self.h.inf, 4])

        t3 = 2 <= t1 + t2  # 2 <= 2x + 0y <= inf
        self.assertEqualExpr(t3, vx, vl, None, [2, self.h.inf])

        t3 = t1 + t2 + 5  # 2x + 0y + 5
        self.assertEqualExpr(t3, vx, vl, 5)

        t3 = 5 + t1 + t2  # 2x + 0y + 5
        self.assertEqualExpr(t3, vx, vl, 5)

        t3 = 5 <= t1 + t2 + 5  # 0 <= 2x + 0y <= inf
        self.assertEqualExpr(t3, vx, vl, None, [0, self.h.inf])

        t3 = 5 + t1 + t2 <= 4  # -inf <= 2x + 0y <= -1
        self.assertEqualExpr(t3, vx, vl, None, [-self.h.inf, -1])

        t3 = 2 <= 5 + t1 + t2  # -3 <= 2x + 0y <= inf
        self.assertEqualExpr(t3, vx, vl, None, [-3, self.h.inf])

        t3 = 2 <= 5 + t1 + t2 <= 6  # -3 <= 2x + 0y <= 1
        self.assertEqualExpr(t3, vx, vl, None, [-3, 1])

        # test chain with variables on both sides
        t3 = 5 + x - x <= y <= 10  # 5 <= y <= 10
        self.assertEqualExpr(t3, [y], [1], None, [5, 10])

        t3 = 5 + 2 * x - 2 * x <= y <= 10  # 5 <= y <= 10
        self.assertEqualExpr(t3, [y], [1], None, [5, 10])

        t3 = 5 <= y <= 10 + 2 * x - 2 * x  # 5 <= y <= 10
        self.assertEqualExpr(t3, [y], [1], None, [5, 10])

        t3 = 5 - x + x <= y <= 10 + y - y  # 5 <= y <= 10
        self.assertEqualExpr(t3, [y], [1], None, [5, 10])

        t3 = 10 >= y >= 5 + x - x  # 5 <= y <= 10
        self.assertEqualExpr(t3, [y], [1], None, [5, 10])

        t3 = 10 >= y >= 5 + 2 * x - 2 * x  # 5 <= y <= 10
        self.assertEqualExpr(t3, [y], [1], None, [5, 10])

        t3 = 10 + 2 * x - 2 * x >= y >= 5  # 5 <= y <= 10
        self.assertEqualExpr(t3, [y], [1], None, [5, 10])

        t3 = 10 + y - y >= y >= 5 - x + x  # 5 <= y <= 10
        self.assertEqualExpr(t3, [y], [1], None, [5, 10])

        vx, vl, nl = list(self.x), [1] * len(self.x), [-1] * len(self.x)
        t3 = qsum(self.x) <= 10  # -inf <= sum(x) <= 10
        self.assertEqualExpr(t3, vx, vl, None, [-self.h.inf, 10])

        t3 = qsum(self.x) <= y  # -inf <= sum(x) - y <= 0
        self.assertEqualExpr(t3, vx + [y], vl + [-1], None, [-self.h.inf, 0])

        t3 = y <= qsum(self.x) <= y  # sum(x) == 0
        self.assertEqualExpr(t3, vx + [y], vl + [-1], None, [0, 0])

        t3 = y >= qsum(self.x) >= y  # sum(x) == 0
        self.assertEqualExpr(t3, vx + [y], vl + [-1], None, [0, 0])

        t3 = qsum(self.x) == y  # sum(x) - y == 0
        self.assertEqualExpr(t3, vx + [y], vl + [-1], None, [0, 0])

        t3 = y == qsum(self.x)  # sum(x) - y == 0
        self.assertEqualExpr(t3, vx + [y], vl + [-1], None, [0, 0])

        t3 = y + 1 <= qsum(self.x) <= y + 6  # 1 <= sum(x) - y <= 6
        self.assertEqualExpr(t3, vx + [y], vl + [-1], None, [1, 6])

        t3 = y + 6 >= qsum(self.x) >= y + 1  # 1 <= sum(x) - y <= 6
        self.assertEqualExpr(t3, vx + [y], vl + [-1], None, [1, 6])

        t3 = qsum(self.x) >= y >= qsum(self.x)  # sum(x) - y == 0
        self.assertEqualExpr(t3, vx + [y], vl + [-1], None, [0, 0])

        t3 = qsum(self.x) <= y <= qsum(self.x)  # sum(x) - y == 0
        self.assertEqualExpr(t3, vx + [y], vl + [-1], None, [0, 0])

        t3 = qsum(self.x) + 5 >= y >= qsum(self.x) + 2  # -5 <= sum(x) - y <= -2
        self.assertEqualExpr(t3, vx + [y], vl + [-1], None, [-5, -2])

        t3 = qsum(self.x) + 2 <= y <= qsum(self.x) + 5  # -5 <= sum(x) - y <= -2
        self.assertEqualExpr(t3, vx + [y], vl + [-1], None, [-5, -2])

        t3 = (qsum(self.x) == 1 + qsum(self.x)).simplify()  # [] == 1
        self.assertEqualExpr(t3, vx, [0] * len(vx), None, [1, 1])

        self.assertRaises(Exception, lambda: x <= y <= 1)
        self.assertRaises(Exception, lambda: (qsum(self.x) + 5 >= y) >= qsum(self.x) + 2)
        self.assertRaises(Exception, lambda: qsum(self.x) + 2 <= (y <= qsum(self.x) + 5))
        self.assertRaises(Exception, lambda: qsum(self.x) + 2 <= y <= 5)
        self.assertRaises(Exception, lambda: 2 <= y <= 5 + qsum(self.x))

    def test_order_priority(self):
        x, y, z = self.x[:3]

        # prefer more variables
        self.assertEqualExpr(x <= y + z, [y, z, x], [1, 1, -1], None, [0, self.h.inf])
        self.assertEqualExpr(y + z >= x, [y, z, x], [1, 1, -1], None, [0, self.h.inf])
        self.assertEqualExpr(x >= y + z, [y, z, x], [1, 1, -1], None, [-self.h.inf, 0])
        self.assertEqualExpr(y + z <= x, [y, z, x], [1, 1, -1], None, [-self.h.inf, 0])
        self.assertEqualExpr(y + z == x, [y, z, x], [1, 1, -1], None, [0, 0])
        self.assertEqualExpr(x == y + z, [y, z, x], [1, 1, -1], None, [0, 0])
        self.assertEqualExpr(x + y == 2 * x + 2 * y + z, [x, y, z, x, y], [2, 2, 1, -1, -1], None, [0, 0])

        # prefer constant
        self.assertEqualExpr(y <= x + 2, [y, x], [1, -1], None, [-self.h.inf, 2])
        self.assertEqualExpr(x + 2 >= y, [y, x], [1, -1], None, [-self.h.inf, 2])
        self.assertEqualExpr(x + 2 <= y, [y, x], [1, -1], None, [2, self.h.inf])
        self.assertEqualExpr(y >= x + 2, [y, x], [1, -1], None, [2, self.h.inf])
        self.assertEqualExpr(y + 2 == x, [x, y], [1, -1], None, [2, 2])
        self.assertEqualExpr(x == y + 2, [x, y], [1, -1], None, [2, 2])

        # prefer 'left'
        self.assertEqualExpr(x + 2 <= y + 3, [x, y], [1, -1], None, [-self.h.inf, 1])
        self.assertEqualExpr(y + 3 >= x + 2, [x, y], [1, -1], None, [-self.h.inf, 1])
        self.assertEqualExpr(x + 2 >= y + 3, [y, x], [1, -1], None, [-self.h.inf, -1])
        self.assertEqualExpr(y + 3 <= x + 2, [y, x], [1, -1], None, [-self.h.inf, -1])
        self.assertEqualExpr(x + 2 == y + 3, [x, y], [1, -1], None, [1, 1])
        self.assertEqualExpr(y + 3 == x + 2, [y, x], [1, -1], None, [-1, -1])

        self.assertEqualExpr(6 <= x + y <= 8, [x, y], [1, 1], None, [6, 8])
        self.assertEqualExpr(6 + x <= y <= 8 + x, [y, x], [1, -1], None, [6, 8])
        self.assertEqualExpr(6 + x <= y + 2 <= 8 + x, [y, x], [1, -1], None, [4, 6])
        self.assertEqualExpr(x <= 6 <= x, [x], [1], None, [6, 6])
        self.assertEqualExpr(x <= y + z <= x + 5, [y, z, x], [1, 1, -1], None, [0, 5])

    def test_chain_inequality_hacks(self):
        x, y = self.x[0:2]

        # Test hacks around chain inequality
        # These don't need to be supported, good to check if logic changes
        t = x + y
        self.assertEqualExpr(t, [x, y], [1, 1])

        bool(t <= 10)
        expr = 5 <= t
        self.assertEqualExpr(expr, [x, y], [1, 1], None, [5, 10])
        expr = 5 <= t
        self.assertEqualExpr(expr, [x, y], [1, 1], None, [5, self.h.inf])

        bool(t <= 10)
        expr = t <= 18
        self.assertEqualExpr(expr, [x, y], [1, 1], None, [-self.h.inf, 18])

        bool(t <= 10)
        expr = t >= 5
        self.assertEqualExpr(expr, [x, y], [1, 1], None, [5, 10])
        expr = t >= 5
        self.assertEqualExpr(expr, [x, y], [1, 1], None, [5, self.h.inf])

        expr = bool(t <= 10) and t >= 5
        self.assertEqualExpr(expr, [x, y], [1, 1], None, [5, 10])

        # test hack if expression is modified
        bool(t <= 10)
        t += y
        expr = 5 <= t
        self.assertEqualExpr(expr, [x, y, y], [1, 1, 1], None, [5, self.h.inf])

        # test hack if using 'equal' temporary expr
        bool(x + y <= 10)
        expr = 5 <= x + y
        self.assertEqualExpr(expr, [x, y], [1, 1], None, [5, 10])
        expr = 5 <= x + y
        self.assertEqualExpr(expr, [x, y], [1, 1], None, [5, self.h.inf])

        # test hack if using 'not-equal' temporary expr
        bool(x + 3 * y <= 10)
        expr = 5 <= x + y
        self.assertEqualExpr(expr, [x, y], [1, 1], None, [5, self.h.inf])

    def test_bounds_already_set(self):
        x, y = self.x[0:2]
        self.assertRaises(Exception, lambda: (1 <= 2 * x + 3 * y) <= 5)
        self.assertRaises(Exception, lambda: (x <= 4) <= 5)
        self.assertRaises(Exception, lambda: 2 <= (x <= 4) <= 5)
        self.assertRaises(Exception, lambda: 2 <= (4 <= x))
        self.assertRaises(Exception, lambda: 2 <= (x <= 4))
        self.assertRaises(Exception, lambda: 4 >= (x >= 2))

        # Cannot a constant if bounds already exist
        self.assertRaises(Exception, lambda: (x + y <= 3) + 5)
        self.assertEqualExpr((x + y <= 3) + (self.h.expr() == 5), [x, y], [1, 1], None, [-self.h.inf, 8])

        # Cannot a expr with constant if bounds already exist
        self.assertRaises(Exception, lambda: (x + y <= 3) + (x + 5))
        self.assertEqualExpr((x + y <= 3) + (x == -5), [x, y, x], [1, 1, 1], None, [-self.h.inf, -2])

    def test_addition(self):
        x, y = self.x[0:2]

        # Test addition of two variables
        expr = x + y
        self.assertEqualExpr(expr, [x, y], [1, 1])

        # Test addition of two exprs (with bounds)
        e1 = 2 <= x - y <= 4
        e2 = 2 <= 2 * x - y <= 4
        self.assertEqualExpr(e1, [x, y], [1, -1], None, [2, 4])
        self.assertEqualExpr(e2, [x, y], [2, -1], None, [2, 4])

        expr = e1 + e2
        self.assertEqualExpr(expr, [x, y, x, y], [1, -1, 2, -1], None, [4, 8])

        expr = (e1 + e2).simplify()
        self.assertEqualExpr(expr, [x, y], [3, -2], None, [4, 8])

        # Test multiplication/addition of two exprs (with bounds)
        expr = (2 * e1 + 3 * e2).simplify()
        self.assertEqualExpr(expr, [x, y], [8, -5], None, [10, 20])

    def test_subtraction(self):
        x, y = self.x[0:2]

        # Test subtraction of two highs_linear_expressions
        expr = x - y
        self.assertEqualExpr(expr, [x, y], [1, -1])

        # Test subtraction of two exprs (with bounds)
        e1 = 2 <= x - y <= 4
        e2 = 2 <= 2 * x - y <= 4
        self.assertEqualExpr(e1, [x, y], [1, -1], None, [2, 4])
        self.assertEqualExpr(e2, [x, y], [2, -1], None, [2, 4])

        expr = e1 + (-1.0 * e2)
        self.assertEqualExpr(expr, [x, y, x, y], [1, -1, -2, 1], None, [-2, 2])

        expr = e1 + (-e2)
        self.assertEqualExpr(expr, [x, y, x, y], [1, -1, -2, 1], None, [-2, 2])

        expr = e1 - e2
        self.assertEqualExpr(expr, [x, y, x, y], [1, -1, -2, 1], None, [-2, 2])

        expr = (e1 - e2).simplify()
        self.assertEqualExpr(expr, [x, y], [-1, 0], None, [-2, 2])

        # Test multiplication/subtraction of two exprs (with bounds)
        expr = (2 * e1 - e2).simplify()
        self.assertEqualExpr(expr, [x, y], [0, -1], None, [0, 6])

        # test rsub
        expr = 5 - (x + y)
        self.assertEqualExpr(expr, [x, y], [-1, -1], 5)

    def test_multiply(self):
        x, y = self.x[0:2]

        # basic tests
        e1 = x * 3
        self.assertEqualExpr(e1, [x], [3])
        e1 = 3 * x
        self.assertEqualExpr(e1, [x], [3])

        e1 = x - y
        self.assertEqualExpr(e1, [x, y], [1, -1])

        e2 = 2 * e1
        self.assertEqualExpr(e2, [x, y], [2, -2])

        e1 *= 3
        self.assertEqualExpr(e1, [x, y], [3, -3])

        e2 = -1 * e1
        self.assertEqualExpr(e2, [x, y], [-3, 3])

        e2 = 0 * e1
        self.assertEqualExpr(e2, [x, y], [0, 0])

        e2 = highs_linear_expression(-1) * e1
        self.assertEqualExpr(e2, [x, y], [-3, 3])

        self.assertRaises(Exception, lambda: e1 * x)
        self.assertRaises(Exception, lambda: e1 * highs_linear_expression(x))
        self.assertRaises(Exception, lambda: e1 * highs_linear_expression(x - x))

        # Test with constant
        e1 = x - y + 4
        self.assertEqualExpr(e1, [x, y], [1, -1], 4)

        e2 = e1 * 2.5
        self.assertEqualExpr(e2, [x, y], [2.5, -2.5], 10)

        e1 *= 2.5
        self.assertEqualExpr(e1, [x, y], [2.5, -2.5], 10)

        e2 = -1.0 * e1
        self.assertEqualExpr(e2, [x, y], [-2.5, 2.5], -10)

        e1 *= -1.0
        self.assertEqualExpr(e1, [x, y], [-2.5, 2.5], -10)

        e1 *= highs_linear_expression(-1.0)
        self.assertEqualExpr(e1, [x, y], [2.5, -2.5], 10)

        # Test with bounds
        e1 = x - 2 * y == [1, 4]
        self.assertEqualExpr(e1, [x, y], [1, -2], None, [1, 4])

        e2 = e1 * 2.5
        self.assertEqualExpr(e2, [x, y], [2.5, -5], None, [2.5, 10])

        e1 *= 2.5
        self.assertEqualExpr(e1, [x, y], [2.5, -5], None, [2.5, 10])

        e2 = -1.0 * e1
        self.assertEqualExpr(e2, [x, y], [-2.5, 5], None, [-10, -2.5])

        e1 *= -1.0
        self.assertEqualExpr(e1, [x, y], [-2.5, 5], None, [-10, -2.5])

        e1 *= highs_linear_expression(-1.0)
        self.assertEqualExpr(e1, [x, y], [2.5, -5], None, [2.5, 10])

        e1 = highs_linear_expression(-1.0)
        e1 *= x + 4
        self.assertEqualExpr(e1, [x], [-1], -4)

        e1 = highs_linear_expression(-1.0)
        e1 *= x - 4 <= 3
        self.assertEqualExpr(e1, [x], [-1], None, [-7, self.h.inf])

        e1 = highs_linear_expression(1.0)
        e1 *= x - 4 <= 3
        self.assertEqualExpr(e1, [x], [1], None, [-self.h.inf, 7])

    def test_simplify(self):
        x, y, z = self.x[0:3]

        # basics
        expr = x + y + z
        self.assertEqualExpr(expr, [x, y, z], [1, 1, 1])
        expr = expr.simplify()
        self.assertEqualExpr(expr, [x, y, z], [1, 1, 1])

        expr = z + x + y  # simplify reorders variables
        self.assertEqualExpr(expr, [z, x, y], [1, 1, 1])
        expr = expr.simplify()
        self.assertEqualExpr(expr, [x, y, z], [1, 1, 1])

        expr = x + x + x + x + x
        self.assertEqualExpr(expr, [x] * 5, [1] * 5)
        expr = expr.simplify()
        self.assertEqualExpr(expr, [x], [5])

        expr = x - x + x - x + x
        self.assertEqualExpr(expr, [x] * 5, [1, -1, 1, -1, 1])
        expr = expr.simplify()
        self.assertEqualExpr(expr, [x], [1])

        expr = -x
        expr += x
        expr *= -5
        self.assertEqualExpr(expr, [x, x], [5, -5])
        expr = expr.simplify()
        self.assertEqualExpr(expr, [x], [0])

        # with constant
        expr = x + y + 1 + y
        self.assertEqualExpr(expr, [x, y, y], [1, 1, 1], 1)
        expr = expr.simplify()
        self.assertEqualExpr(expr, [x, y], [1, 2], 1)

        # with bounds
        expr = 0 <= x + y + 1 + y <= 5
        self.assertEqualExpr(expr, [x, y, y], [1, 1, 1], None, [-1, 4])
        expr = expr.simplify()
        self.assertEqualExpr(expr, [x, y], [1, 2], None, [-1, 4])

    def test_evaluate(self):
        h = self.h
        x, y = h.addVariables(2, lb=-h.inf)

        h.addConstrs(-x + y >= 2, x + y >= 0)
        h.minimize(y)

        self.assertAlmostEqual(h.val(x), -1)
        self.assertAlmostEqual(h.val(y), 1)
        self.assertAlmostEqual(h.val(x + y), 0)
        self.assertAlmostEqual(h.val(y - x), 2)
        self.assertEqual(h.val(y - x >= 2), True)
        self.assertEqual(h.val(y - x <= 1), False)

    def test_repr_str(self):
        x, y = self.x[0:2]
        self.assertEqual(repr(x), "highs_var(0)")
        self.assertEqual(str(x), "highs_var(0)")

        c = self.h.addConstr(x + y <= 5)
        self.assertEqual(repr(c), "highs_cons(0)")
        self.assertEqual(str(c), "highs_cons(0)")

        expr = c.expr()
        self.assertEqual(repr(expr), "-inf <= 1.0_v0  1.0_v1 <= 5.0")
        self.assertEqual(str(expr), "-inf <= 1.0_v0  1.0_v1 <= 5.0")

        expr = x + y
        self.assertEqual(repr(expr), "1.0_v0  1.0_v1")
        self.assertEqual(str(expr), "1.0_v0  1.0_v1")

        expr = x + y + x + x + y
        self.assertEqual(repr(expr), "1.0_v0  1.0_v1  1.0_v0  1.0_v0  1.0_v1")
        self.assertEqual(str(expr), "3.0_v0  2.0_v1")

        expr = y + x + 5
        self.assertEqual(repr(expr), "1.0_v1  1.0_v0  5.0")
        self.assertEqual(str(expr), "1.0_v0  1.0_v1  5.0")

        expr = y + x == 5
        self.assertEqual(repr(expr), "1.0_v1  1.0_v0 == 5.0")
        self.assertEqual(str(expr), "1.0_v0  1.0_v1 == 5.0")

        expr = y + x <= 5
        self.assertEqual(repr(expr), "-inf <= 1.0_v1  1.0_v0 <= 5.0")
        self.assertEqual(str(expr), "-inf <= 1.0_v0  1.0_v1 <= 5.0")
