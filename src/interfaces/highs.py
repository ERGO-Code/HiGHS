#  ___________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright 2017 National Technology and Engineering Solutions of Sandia, LLC
#  Under the terms of Contract DE-NA0003525 with National Technology and
#  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain
#  rights in this software.
#  This software is distributed under the 3-clause BSD License.
#  ___________________________________________________________________________

import logging
import re
import sys
import csv
import subprocess

from pyomo.common.tempfiles import TempfileManager

from pyomo.common import Executable
from pyomo.common.collections import Bunch
from pyomo.core.base.component import name
from pyomo.opt import SolverFactory, OptSolver, ProblemFormat, ResultsFormat, SolverResults, TerminationCondition, SolutionStatus, ProblemSense
from pyomo.opt.base.solvers import _extract_version
from pyomo.opt.solver import SystemCallSolver

logger = logging.getLogger('pyomo.solvers')

_highs_version = None


def configure_highs():
    global _highs_version
    if _highs_version is not None:
        return
    _highs_version = _extract_version("")
    if not Executable("highs"):
        return
    result = subprocess.run([Executable('highs').path()],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=1, universal_newlines=True)
    if not result.returncode:
        _highs_version = _extract_version(result.stdout)


@SolverFactory.register('highs', doc='The GLPK LP/MIP solver')
class GLPK(OptSolver):
    """The GLPK LP/MIP solver"""

    def __new__(cls, *args, **kwds):
        configure_highs()
        try:
            mode = kwds['solver_io']
            if mode is None:
                mode = 'lp'
            del kwds['solver_io']
        except KeyError:
            mode = 'lp'
        
        if mode  == 'lp':           
            return SolverFactory('_highs_shell', **kwds)
        elif mode == 'mps':
            raise ValueError(
                'This model does not support')
        elif mode == 'python':
            raise ValueError(
               'This model does not support')
        else:
            raise ValueError(
                'This model does not support')


@SolverFactory.register('_highs_shell',doc='Shell interface to the GNU Linear Programming Kit')
class highsSHELL(SystemCallSolver):
    """Shell interface to the highs LP/MIP solver"""

    def __init__(self, **kwargs):
        configure_highs()
        #
        # Call base constructor.
        #
        kwargs['type'] = 'highs'
        SystemCallSolver.__init__(self, **kwargs)

        self._rawfile = None
        self._log_file = None
        #
        # Valid problem formats, and valid results for each format.
        #
        self._valid_problem_formats = [ProblemFormat.cpxlp,
                                       ProblemFormat.mps, ]
        self._valid_result_formats = {
            ProblemFormat.cpxlp: ResultsFormat.soln,
            ProblemFormat.mps:   ResultsFormat.soln,
        }
        self.set_problem_format(ProblemFormat.cpxlp)

        # Note: Undefined capabilities default to 'None'.
        self._capabilities = Bunch()
        self._capabilities.linear = True
        self._capabilities.integer = True

    def _default_results_format(self, prob_format):
        return ResultsFormat.soln

    def _default_executable(self):
        executable = Executable('highs')
        if not executable:
            msg = ("Could not locate the 'highs' executable, which is "
                   "required for solver '%s'")
            logger.warning(msg % self.name)
            self.enable = False
            return None
        return executable.path()

    def _get_version(self):
        """
        Returns a tuple describing the solver executable version.
        """
        if _highs_version is None:
            return _extract_version('')
        return _highs_version

    def create_command_line(self, executable, problem_files):
        self._log_file = 'Highs.log'
        #
        # Define solution file.
        #
        self._lpfile = problem_files[0]
        self._rawfile = TempfileManager.create_tempfile(suffix='.highs.raw')
        self._soln_file = self._rawfile
        #
        # Define command line.
        #
        cmd = [executable]
        if self._timer:
            cmd.insert(0, self._timer)
        for key in self.options:
            opt = self.options[key]
            if opt is None or (isinstance(opt, str) and opt.strip() == ''):
                # Handle the case for options that must be
                # specified without a value.
                cmd.append("--%s" % key)
            else:
                cmd.extend(["--%s" % key, str(opt)])

        if self._timelimit is not None and self._timelimit > 0.0:
            cmd.extend(['--time_limit', str(self._timelimit)])

        cmd.extend(['--solution_file', self._rawfile])
        cmd.extend(['--model_file', problem_files[0]])

        return Bunch(cmd=cmd, log_file=self._log_file, env=None)

    def process_logfile(self):
        """
        Process logfile
        """
        results = SolverResults()

        # For the lazy programmer, handle long variable names.
        prob = results.problem
        solv = results.solver
        solv.termination_condition = TerminationCondition.unknown
        stats = results.solver.statistics
        bbound = stats.branch_and_bound

        prob.upper_bound = float('inf')
        prob.lower_bound = float('-inf')
        bbound.number_of_created_subproblems = 0
        bbound.number_of_bounded_subproblems = 0
        with open(self._log_file, 'r') as output:
            for line in output:
                toks = line.split()
                if 'tree is empty' in line:
                    bbound.number_of_created_subproblems = toks[-1][:-1]
                    bbound.number_of_bounded_subproblems = toks[-1][:-1]
                elif len(toks) == 4 and toks[1] == 'status':
                    if toks[3] == 'Optimal':
                        solv.termination_condition = TerminationCondition.optimal
                elif len(toks) == 9 and toks[0] == 'LP':
                    prob.number_of_constraints = int(toks[3])
                    prob.number_of_nonzeros = int(toks[7])
                    prob.number_of_variables = int(toks[5])
                elif len(toks) >= 1 and toks[0] == 'Objective':
                    self.obj_val = float(toks[3])
        return results

    def _highs_get_solution_status(self, status):
        # if GLP_FEAS     == status: return SolutionStatus.feasible
        # elif GLP_INFEAS == status: return SolutionStatus.infeasible
        # elif GLP_NOFEAS == status: return SolutionStatus.infeasible
        # elif GLP_UNDEF  == status: return SolutionStatus.other
        # elif GLP_OPT    == status: return SolutionStatus.optimal
        return SolutionStatus.optimal

    def _process_lp_file(self, row, reader, results, variable_names, constraint_names):

        varname_dict = {} 

        cid = 0
        while True:
            row = next(reader)
            if len(row) == 0:
                continue
            if row[0] == 's.t.':
                break
            _, vname = row
            cid += 1
            variable_names[cid] = (vname)
            varname_dict[vname] = True

        rid = 0
        while True:
            row = next(reader)
            if len(row) == 0:
                continue
            if row[0] == 'bounds':
                break
            if ':' in row[0] and 'ONE_VAR_CONSTANT' not in row[0]:
                cname = row[0]
                rid += 1
                constraint_names[rid] = (cname[:-1])
        while True:
            row = next(reader)
            
            if len(row) == 0:
                continue
            if row[0] == 'end':
                break
            vname = row[5]
            if vname not in varname_dict:
                cid += 1
                variable_names[cid] = (vname)    
        self.max_cid = cid
        self.max_rid = rid

    def process_soln_file(self, results):

        pdata = self._lpfile
        psoln = self._rawfile

        prob = results.problem
        solv = results.solver

        prob.name = 'unknown'   # will ostensibly get updated.

        # Step 1: Make use of the highs's machine parseable format (--wglp) to
        #    collect variable and constraint names.
        glp_line_count = ' -- File not yet opened'

        # The trick for getting the variable names correctly matched to their
        # values is the note that the --wglp option outputs them in the same
        # order as the --write output.
        # Note that documentation for these formats is available from the highs
        # documentation of 'glp_read_prob' and 'glp_write_sol.

        variable_names = dict()    # cols
        constraint_names = dict()  # rows
        obj_name = 'objective'

        with open(pdata, 'r') as csvfile:

            reader = csv.reader(csvfile, delimiter=' ')

            row = next(reader)

            try:
                row = next(reader)
                row = next(reader)
                if len(row) == 1 and (row[0] == 'min' or row[0] == 'max'):
                    prob.sense = 'min' == row[0] and ProblemSense.minimize or ProblemSense.maximize
                row = next(reader)
                obj_name = row[0][:-1]

                self._process_lp_file(row, reader, results,
                                     variable_names, constraint_names)


            except Exception:
                print("ERROR: " + str(sys.exc_info()[1]))
                msg = "Error parsing solution data file, line %d" % reader.line_num
                raise ValueError(msg)
        # Step 2: Make use of the highs's machine parseable format (--write) to
        #    collect solution variable and constraint values.
        with open(psoln, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=' ')
            row = next(reader)
            try:
                row = next(reader)
                while (row[0] != 'Columns'):
                    row = next(reader)

                self._process_soln_bas(
                    row, reader, results, obj_name, variable_names, constraint_names)

            except Exception:
                print("ERROR: " + str(sys.exc_info()[1]))
                msg = "Error parsing solution data file, line %d" % reader.line_num
                raise ValueError(msg)

    def _process_soln_bas(self, row, reader, results, obj_name, variable_names, constraint_names):
        """
        Process a basic solution
        """

        solv = results.solver

        if True:
            soln = results.solution.add()
            soln.status = SolutionStatus.feasible
            solv.termination_condition = TerminationCondition.optimal

            # TODO: Should we have a gap value for LP solves?
            soln.gap = 0.0
            results.problem.lower_bound = self.obj_val
            results.problem.upper_bound = self.obj_val

            # I'd like to choose the correct answer rather than just doing
            # something like commenting the obj_name line.  The point is that
            # we ostensibly could or should make use of the user's choice in
            # objective name.  In that vein I'd like to set the objective value
            # to the objective name.  This would make parsing on the user end
            # less 'arbitrary', as in the yaml key 'f'.  Weird.
            soln.objective[obj_name] = {'Value': self.obj_val}

            extract_duals = False
            extract_reduced_costs = False
            for suffix in self._suffixes:
                if re.match(suffix, "dual"):
                    extract_duals = True
                elif re.match(suffix, "rc"):
                    extract_reduced_costs = True

            range_duals = {}

            cid = 0
            while True:
                row = next(reader)
                cid += 1
                if cid > self.max_cid:
                    break
                # NOTE: we are not using the column status (cst) value right now
                cprim, cdual, rtype = row

                vname = variable_names[int(cid)]

                cprim = float(cprim)
                if extract_reduced_costs is False:
                    soln.variable[vname] = {"Value": cprim}
                else:
                    soln.variable[vname] = {"Value": cprim, "Rc": float(cdual)}
            rid = 0
            row = next(reader)

            while True:
                row = next(reader)
                rid += 1
                if rid > self.max_rid:
                    break
                if not extract_duals:
                    continue

                # NOTE: we are not using the row status (rst) value right now
                rprim, rdual, rtype = row

                cname = constraint_names[int(rid)]

                rdual = float(rdual)
                if cname.startswith('c_'):
                    soln.constraint[cname] = {"Dual": rdual}
                elif cname.startswith('r_l_'):
                    range_duals.setdefault(cname[4:], [0, 0])[0] = rdual
                elif cname.startswith('r_u_'):
                    range_duals.setdefault(cname[4:], [0, 0])[1] = rdual

            # For the range constraints, supply only the dual with the largest
            # magnitude (at least one should always be numerically zero).
            scon = soln.Constraint
            for key, (ld, ud) in range_duals.items():
                if abs(ld) > abs(ud):
                    scon['r_l_'+key] = {"Dual": ld}
                else:
                    scon['r_l_'+key] = {"Dual": ud}      # Use the same key
