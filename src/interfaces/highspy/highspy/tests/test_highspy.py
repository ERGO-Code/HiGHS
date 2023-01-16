import unittest
import highspy
import numpy as np
from pyomo.common.tee import capture_output
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
        self.assertEqual(h.version(), "1.5.0")
        self.assertEqual(h.versionMajor(), 1)
        self.assertEqual(h.versionMinor(), 5)
        self.assertEqual(h.versionPatch(), 0)

    def test_basics(self):
        h = self.get_basic_model()
        self.assertEqual(h.setOptionValue('presolve', 'off'), highspy.HighsStatus.kOk)
#        h.setOptionValue('output_flag', True)
        h.run()
        [status, valid, integral, feasible] = h.assessPrimalSolution()
        self.assertEqual(status, highspy.HighsStatus.kOk)
        self.assertEqual(valid, True)
        self.assertEqual(integral, True)
        self.assertEqual(feasible, True)
        
        # Info can be obtained from the class instance, specific call
        # and, in the case of objective_function_value,
        # h.getObjectiveValue()
        info = h.getInfo()
        objective_function_value0 = info.objective_function_value
        self.assertAlmostEqual(objective_function_value0, 1)
        [status, objective_function_value1] = h.getDoubleInfoValue("objective_function_value")
        self.assertAlmostEqual(objective_function_value0, objective_function_value1)
        self.assertAlmostEqual(h.getObjectiveValue(), objective_function_value0)

        simplex_iteration_count0 = info.simplex_iteration_count
        self.assertAlmostEqual(simplex_iteration_count0, 2)
        [status, simplex_iteration_count1] = h.getIntInfoValue("simplex_iteration_count")
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
        [status, mip_node_count1] = h.getDoubleInfoValue("mip_node_count")
        self.assertEqual(status, highspy.HighsStatus.kError)
        [status, mip_node_count1] = h.getInt64InfoValue("mip_node_count")
        self.assertEqual(status, highspy.HighsStatus.kOk)
        self.assertAlmostEqual(mip_node_count0, mip_node_count1)

    def test_options(self):
        h = highspy.Highs()
        # test bool option
        [status, type] = h.getOptionType('output_flag')
        self.assertEqual(type, highspy.HighsOptionType.kBool)

        h.setOptionValue('output_flag', True)
        [status, value] = h.getBoolOptionValue('output_flag')
        self.assertTrue(value)
        h.setOptionValue('output_flag', False)
        [status, value] = h.getBoolOptionValue('output_flag')
        self.assertFalse(value)

        # test string option
        [status, type] = h.getOptionType('presolve')
        self.assertEqual(type, highspy.HighsOptionType.kString)
        h.setOptionValue('presolve', 'off')
        [status, value] = h.getStringOptionValue('presolve')
        self.assertEqual(value, 'off')
        h.setOptionValue('presolve', 'on')
        [status, value] = h.getStringOptionValue('presolve')
        self.assertEqual(value, 'on')

        # test int option
        [status, type] = h.getOptionType('threads')
        self.assertEqual(type, highspy.HighsOptionType.kInt)
        h.setOptionValue('threads', 1)
        [status, value] = h.getIntOptionValue('threads')
        self.assertEqual(value, 1)
        h.setOptionValue('threads', 2)
        [status, value] = h.getIntOptionValue('threads')
        self.assertEqual(value, 2)

        # test double option
        [status, type] = h.getOptionType('time_limit')
        self.assertEqual(type, highspy.HighsOptionType.kDouble)
        h.setOptionValue('time_limit', 1.7)
        [status, value] = h.getDoubleOptionValue('time_limit')
        self.assertAlmostEqual(value, 1.7)
        h.setOptionValue('time_limit', 2.7)
        [status, value] = h.getDoubleOptionValue('time_limit')
        self.assertAlmostEqual(value, 2.7)

    def test_clear(self):
        h = self.get_basic_model()
        self.assertEqual(h.getNumCol(), 2)
        self.assertEqual(h.getNumRow(), 2)
        self.assertEqual(h.getNumNz(), 4)

        [status, orig_feas_tol] = h.getDoubleOptionValue('primal_feasibility_tolerance')
        new_feas_tol = orig_feas_tol + 1
        h.setOptionValue('primal_feasibility_tolerance', new_feas_tol)
        [status, value] = h.getDoubleOptionValue('primal_feasibility_tolerance')
        self.assertAlmostEqual(value, new_feas_tol)
        h.clear()
        self.assertEqual(h.getNumCol(), 0)
        self.assertEqual(h.getNumRow(), 0)
        self.assertEqual(h.getNumNz(), 0)
        [status, value] = h.getDoubleOptionValue('primal_feasibility_tolerance')
        self.assertAlmostEqual(value, orig_feas_tol)

        h = self.get_basic_model()
        h.setOptionValue('primal_feasibility_tolerance', new_feas_tol)
        [status, value] = h.getDoubleOptionValue('primal_feasibility_tolerance')
        self.assertAlmostEqual(value, new_feas_tol)
        h.clearModel()
        self.assertEqual(h.getNumCol(), 0)
        self.assertEqual(h.getNumRow(), 0)
        self.assertEqual(h.getNumNz(), 0)
        [status, value] = h.getDoubleOptionValue('primal_feasibility_tolerance')
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
        [status, orig_feas_tol] = h.getDoubleOptionValue('primal_feasibility_tolerance')
        new_feas_tol = orig_feas_tol + 1
        h.setOptionValue('primal_feasibility_tolerance', new_feas_tol)
        [status, value] = h.getDoubleOptionValue('primal_feasibility_tolerance')
        self.assertAlmostEqual(value, new_feas_tol)
        h.resetOptions()
        [status, value] = h.getDoubleOptionValue('primal_feasibility_tolerance')
        self.assertAlmostEqual(value, orig_feas_tol)

    def test_dual_ray(self):
        h = self.get_infeasible_model()
        h.run()
        h.setOptionValue('output_flag', True)
        # Check that there is a dual ray
        [status, has_dual_ray] = h.getDualRay()
        self.assertTrue(has_dual_ray)
        
        num_row = h.getLp().num_row_
        values = np.empty(num_row, dtype=np.double)
        h.getDualRay(values)
        self.assertAlmostEqual(values[0], 0.5)
        self.assertAlmostEqual(values[1], -1)
        
        v = h.foo(values)
        self.assertEqual(v, 1) #Should be -1
        v = h.foo()
        self.assertEqual(v, 1)
 
    def test_check_solution_feasibility(self):
        h = self.get_basic_model()
        [status, valid, integral, feasible] = h.assessPrimalSolution()
        self.assertEqual(status, highspy.HighsStatus.kError)
        h.run()
        [status, valid, integral, feasible] = h.assessPrimalSolution()
        self.assertEqual(valid, True)
        self.assertEqual(integral, True)
        self.assertEqual(feasible, True)

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

        
    def test_basic_solves(self):
        h = self.get_basic_model()
        h.run()
        num_row = h.getLp().num_row_
        basic_variables = np.empty(num_row, dtype=np.intc)
        h.getBasicVariables(basic_variables)
        self.assertEqual(basic_variables[0], 1)
        self.assertEqual(basic_variables[1], 0)
        row = 0
        row_vector = np.empty(num_row, dtype=np.double)
        row_num_nz = np.empty(1, dtype=np.intc)
        row_indices = np.empty(num_row, dtype=np.intc)
        h.getBasisInverseRow(row, row_vector, row_num_nz, row_indices)#None, None)

    def test_log_callback(self):
        h = self.get_basic_model()
        h.setOptionValue('output_flag', True)

        class Foo(object):
            def __str__(self):
                return 'an instance of Foo'

            def __repr__(self):
                return self.__str__()

        def log_callback(log_type, message, data):
            print('got a log message: ', log_type, data, message)

        h.setLogCallback(log_callback, Foo())
        out = StringIO()
        with capture_output(out) as t:
            h.run()
        out = out.getvalue()
        self.assertIn('got a log message:  HighsLogType.kInfo an instance of Foo Presolving model', out)

