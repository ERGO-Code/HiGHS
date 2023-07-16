import tempfile
import unittest
import highspy
import numpy as np
from io import StringIO


class TestHighsPy(unittest.TestCase):
    def get_basic_model(self):
        """
        min y
        s.t.
        -x + y >= 2
        x + y >= 0
        """
        inf = highspy.kHighsInf
        h = highspy.Highs()
        h.setOptionValue('output_flag', False)
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
        lp.num_col_ = 2;
        lp.num_row_ = 3;
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
    
    def get_infeasible_model(self):
        inf = highspy.kHighsInf
        lp = highspy.HighsLp()
        lp.num_col_ = 2;
        lp.num_row_ = 2;
        lp.col_cost_ = np.array([10, 15], dtype=np.double)
        lp.col_lower_ = np.array([0, 0], dtype=np.double)
        lp.col_upper_ = np.array([inf, inf], dtype=np.double)
        lp.row_lower_ = np.array([3, 1], dtype=np.double)
        lp.row_upper_ = np.array([3, 1], dtype=np.double)
        lp.a_matrix_.start_ = np.array([0, 2, 4])
        lp.a_matrix_.index_ = np.array([0, 1, 0, 1])
        lp.a_matrix_.value_ = np.array([2, 1, 1, 3], dtype=np.double)
        lp.offset_ = 0;
        h = highspy.Highs()
        h.setOptionValue('output_flag', False)
        status = h.passModel(lp)
        self.assertEqual(status, highspy.HighsStatus.kOk)
        h.setOptionValue('presolve', 'off')
        return h
    
    def test_version(self):
        h = self.get_basic_model()
        self.assertEqual(h.version(), "1.5.3")
        self.assertEqual(h.versionMajor(), 1)
        self.assertEqual(h.versionMinor(), 5)
        self.assertEqual(h.versionPatch(), 3)

    def test_basics(self):
        h = self.get_basic_model()
        h.passColName(0, 'Col0')
        h.passColName(1, 'Col1')
        h.passRowName(0, 'Row0')
        h.passRowName(1, 'Row1')
#        h.setOptionValue('output_flag', True)
        h.writeModel('')
        h.setOptionValue('output_flag', False)
        self.assertEqual(h.setOptionValue('presolve', 'off'), highspy.HighsStatus.kOk)
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
        h.addRows(1, np.array([0], dtype=np.double), np.array([inf]), 2,
                  np.array([0]), np.array([0, 1]), np.array([-1, 1], dtype=np.double))
        h.run()
        sol = h.getSolution()
        self.assertAlmostEqual(sol.col_value[0], 0)
        self.assertAlmostEqual(sol.col_value[1], 0)

        # change the upper bound of x to -5
        h.changeColsBounds(1, np.array([0]), np.array([-inf], dtype=np.double),
                           np.array([-5], dtype=np.double))
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
        [status, offset] = h.getObjectiveOffset();
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

        [status, output_flag] = h.getOptionValue('output_flag')
        [status, solver] = h.getOptionValue('solver')
        [status, primal_feasibility_tolerance] = h.getOptionValue('primal_feasibility_tolerance')
        [status, simplex_update_limit] = h.getOptionValue('simplex_update_limit')
        self.assertEqual(output_flag, True);
        self.assertEqual(solver, 'choose');
        self.assertEqual(primal_feasibility_tolerance, 1e-7);
        self.assertEqual(simplex_update_limit, 5000);
        # Illegal name
        option_value = h.getOptionValue('simplex_limit')
        self.assertEqual(option_value[0], highspy.HighsStatus.kError)
        
        # test bool option
        [status, type] = h.getOptionType('output_flag')
        self.assertEqual(type, highspy.HighsOptionType.kBool)

        h.setOptionValue('output_flag', True)
        [status, value] = h.getOptionValue('output_flag')
        self.assertTrue(value)
        h.setOptionValue('output_flag', False)
        [status, value] = h.getOptionValue('output_flag')
        self.assertFalse(value)

        # test string option
        [status, type] = h.getOptionType('presolve')
        self.assertEqual(type, highspy.HighsOptionType.kString)
        h.setOptionValue('presolve', 'off')
        [status, value] = h.getOptionValue('presolve')
        self.assertEqual(value, 'off')
        h.setOptionValue('presolve', 'on')
        [status, value] = h.getOptionValue('presolve')
        self.assertEqual(value, 'on')

        # test int option
        [status, type] = h.getOptionType('threads')
        self.assertEqual(type, highspy.HighsOptionType.kInt)
        h.setOptionValue('threads', 1)
        [status, value] = h.getOptionValue('threads')
        self.assertEqual(value, 1)
        h.setOptionValue('threads', 2)
        [status, value] = h.getOptionValue('threads')
        self.assertEqual(value, 2)

        # test double option
        [status, type] = h.getOptionType('time_limit')
        self.assertEqual(type, highspy.HighsOptionType.kDouble)
        h.setOptionValue('time_limit', 1.7)
        [status, value] = h.getOptionValue('time_limit')
        self.assertAlmostEqual(value, 1.7)
        h.setOptionValue('time_limit', 2.7)
        [status, value] = h.getOptionValue('time_limit')
        self.assertAlmostEqual(value, 2.7)

    def test_clear(self):
        h = self.get_basic_model()
        self.assertEqual(h.getNumCol(), 2)
        self.assertEqual(h.getNumRow(), 2)
        self.assertEqual(h.getNumNz(), 4)

        [status, orig_feas_tol] = h.getOptionValue('primal_feasibility_tolerance')
        new_feas_tol = orig_feas_tol + 1
        h.setOptionValue('primal_feasibility_tolerance', new_feas_tol)
        [status, value] = h.getOptionValue('primal_feasibility_tolerance')
        self.assertAlmostEqual(value, new_feas_tol)
        h.clear()
        self.assertEqual(h.getNumCol(), 0)
        self.assertEqual(h.getNumRow(), 0)
        self.assertEqual(h.getNumNz(), 0)
        [status, value] = h.getOptionValue('primal_feasibility_tolerance')
        self.assertAlmostEqual(value, orig_feas_tol)

        h = self.get_basic_model()
        h.setOptionValue('primal_feasibility_tolerance', new_feas_tol)
        [status, value] = h.getOptionValue('primal_feasibility_tolerance')
        self.assertAlmostEqual(value, new_feas_tol)
        h.clearModel()
        self.assertEqual(h.getNumCol(), 0)
        self.assertEqual(h.getNumRow(), 0)
        self.assertEqual(h.getNumNz(), 0)
        [status, value] = h.getOptionValue('primal_feasibility_tolerance')
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
        [status, orig_feas_tol] = h.getOptionValue('primal_feasibility_tolerance')
        new_feas_tol = orig_feas_tol + 1
        h.setOptionValue('primal_feasibility_tolerance', new_feas_tol)
        [status, value] = h.getOptionValue('primal_feasibility_tolerance')
        self.assertAlmostEqual(value, new_feas_tol)
        h.resetOptions()
        [status, value] = h.getOptionValue('primal_feasibility_tolerance')
        self.assertAlmostEqual(value, orig_feas_tol)

    def test_ranging(self):
        inf = highspy.kHighsInf
        h = self.get_basic_model()
        # Cost ranging
        #c0 2 -1 1 0
        #c1 0 0 inf inf
        #
        ## Bound ranging
        ## Columns
        #c0 1 -inf inf 1
        #c1 1 1 inf 1
        ## Rows
        #r0 -inf -inf inf inf
        #r1 -inf -inf inf inf
        h.run()
        [status, ranging] = h.getRanging()
        self.assertEqual(ranging.col_cost_dn.objective_[0], 2);
        self.assertEqual(ranging.col_cost_dn.value_[0], -1);
        self.assertEqual(ranging.col_cost_up.value_[0], 1);
        self.assertEqual(ranging.col_cost_up.objective_[0], 0);
        self.assertEqual(ranging.col_cost_dn.objective_[1], 0);
        self.assertEqual(ranging.col_cost_dn.value_[1], 0);
        self.assertEqual(ranging.col_cost_up.value_[1], inf);
        self.assertEqual(ranging.col_cost_up.objective_[1], inf);
#
        self.assertEqual(ranging.col_bound_dn.objective_[0], 1);
        self.assertEqual(ranging.col_bound_dn.value_[0], -inf);
        self.assertEqual(ranging.col_bound_up.value_[0], inf);
        self.assertEqual(ranging.col_bound_up.objective_[0], 1);
        self.assertEqual(ranging.col_bound_dn.objective_[1], 1);
        self.assertEqual(ranging.col_bound_dn.value_[1], 1);
        self.assertEqual(ranging.col_bound_up.value_[1], inf);
        self.assertEqual(ranging.col_bound_up.objective_[1], 1);
#
        self.assertEqual(ranging.row_bound_dn.objective_[0], -inf);
        self.assertEqual(ranging.row_bound_dn.value_[0], -inf);
        self.assertEqual(ranging.row_bound_up.value_[0], inf);
        self.assertEqual(ranging.row_bound_up.objective_[0], inf);
        self.assertEqual(ranging.row_bound_dn.objective_[1], -inf);
        self.assertEqual(ranging.row_bound_dn.value_[1], -inf);
        self.assertEqual(ranging.row_bound_up.value_[1], inf);
        self.assertEqual(ranging.row_bound_up.objective_[1], inf);
        
    def test_write_basis_before_running(self):
        h = self.get_basic_model()
        with tempfile.NamedTemporaryFile() as f:
            h.writeBasis(f.name)
            contents = f.read()
            self.assertEqual(contents, b'HiGHS v1\nNone\n')
        
    def test_write_basis_after_running(self):
        h = self.get_basic_model()
        h.run()
        with tempfile.NamedTemporaryFile() as f:
            h.writeBasis(f.name)
            contents = f.read()
            self.assertEqual(
                contents, b'HiGHS v1\nValid\n# Columns 2\n1 1 \n# Rows 2\n0 0 \n'
            )

    def test_read_basis(self):
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
