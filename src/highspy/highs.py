from ._core import (
    # enum classes
    ObjSense,
    MatrixFormat,
    HessianFormat,
    SolutionStatus,
    BasisValidity,
    HighsModelStatus,
    HighsPresolveStatus,
    HighsBasisStatus,
    HighsVarType,
    HighsOptionType,
    HighsInfoType,
    HighsStatus,
    HighsLogType,
    # classes
    HighsSparseMatrix,
    HighsLp,
    HighsHessian,
    HighsModel,
    HighsInfo,
    HighsOptions,
    _Highs,
    # structs
    HighsSolution,
    HighsObjectiveSolution,
    HighsBasis,
    HighsRangingRecord,
    HighsRanging,
    # constants
    kHighsInf,
    kHighsIInf,
)

from collections.abc import Mapping
from itertools import product
from decimal import Decimal
from threading import local
import numpy as np

class Highs(_Highs):
    """
    HiGHS solver interface
    """
    def __init__(self):
        super().__init__()

    # Silence logging
    def silent(self, turn_off_output=True):
        """
        Disables solver output to the console.
        """
        super().setOptionValue("output_flag", not turn_off_output)
    
    # solve
    def solve(self):
        """Runs the solver on the current problem.

        Returns:
            A HighsStatus object containing the solve status.
        """
        return super().run()

    def optimize(self):
        """
        Alias for the solve method.
        """
        return super().run()

    # reset the objective and sense, then solve
    def minimize(self, obj=None):
        """Solves a minimization of the objective and optionally updates the costs.

        Args:
            obj: An optional highs_linear_expression representing the new objective function.

        Raises:
            Exception: If obj is an inequality or not a highs_linear_expression.

        Returns:
            A HighsStatus object containing the solve status after minimization.
        """
        if obj is not None:
            # if we have a single variable, wrap it in a linear expression
            if isinstance(obj, highs_var) == True:
                obj = highs_linear_expression(obj)

            elif isinstance(obj, highs_linear_expression) == False or obj.bounds != None:
                raise Exception('Objective cannot be an inequality') 

            # reset objective
            super().changeColsCost(self.numVariables, range(self.numVariables), [0]*self.numVariables)

            # if we have duplicate variables, add the vals
            vars,vals = obj.unique_elements()
            super().changeColsCost(len(vars), vars, vals)
            super().changeObjectiveOffset(obj.constant or 0.0)

        super().changeObjectiveSense(ObjSense.kMinimize)
        return super().run()

    # reset the objective and sense, then solve
    def maximize(self, obj=None):
        """Solves a maximization of the objective and optionally updates the costs.

        Args:
            obj: An optional highs_linear_expression representing the new objective function.

        Raises:
            Exception: If obj is an inequality or not a highs_linear_expression.

        Returns:
            A HighsStatus object containing the solve status after maximization.
        """
        if obj is not None:
            # if we have a single variable, wrap it in a linear expression
            if isinstance(obj, highs_var) == True:
                obj = highs_linear_expression(obj)

            elif isinstance(obj, highs_linear_expression) == False or obj.bounds != None:
                raise Exception('Objective cannot be an inequality') 

            # reset objective
            super().changeColsCost(self.numVariables, range(self.numVariables), [0]*self.numVariables)

            # if we have duplicate variables, add the vals
            vars,vals = obj.unique_elements()
            super().changeColsCost(len(vars), vars, vals)
            super().changeObjectiveOffset(obj.constant or 0.0)

        super().changeObjectiveSense(ObjSense.kMaximize)
        return super().run()

    def internal_get_value(self, array_values, index_collection):
        """
        Internal method to get the value of an index from an array of values. Could be value or dual, variable or constraint.
        """
        if isinstance(index_collection, (int, highs_var, highs_cons)):
            return array_values[int(index_collection)]

        elif isinstance(index_collection, highs_linear_expression):
            return index_collection.evaluate(array_values)

        elif isinstance(index_collection, Mapping):
            return {k: self.internal_get_value(array_values, v) for k,v in index_collection.items()}

        else:
            return np.asarray([self.internal_get_value(array_values, v) for v in index_collection])

    def val(self, var):
        """Gets the value of a variable/index or expression in the solution.

        Args:
            var: A highs_var/index or highs_linear_expression object representing the variable.

        Returns:
            The value of the variable in the solution.
        """
        return self.internal_get_value(super().getSolution().col_value, var)

    def vals(self, vars):
        """Gets the values of multiple variables in the solution.

        Args:
            vars: A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.

        Returns:
            If vars is a Mapping, returns a dict where keys are the same keys from the input vars and values are the solution values of the corresponding variables. If vars is an iterable, returns a list of solution values for the variables.
        """
        return self.internal_get_value(super().getSolution().col_value, vars)

    def variableName(self, var):
        """Retrieves the name of a specific variable.

        Args:
            var: A highs_var object representing the variable.

        Raises:
            Exception: If the variable name cannot be found.

        Returns:
            The name of the specified variable.
        """
        [status, name] = super().getColName(int(var))
        failed = status != HighsStatus.kOk
        if failed:
            raise Exception('Variable name not found') 
        return name

    def variableNames(self, vars):
        """Retrieves the names of multiple variables.

        Args:
            vars: An iterable of highs_var objects or a mapping where keys are identifiers and values are highs_var objects.

        Raises:
            Exception: If any variable name cannot be found.

        Returns:
            If vars is a mapping, returns a dict where keys are the same keys from the input vars and values are the names of the corresponding variables. 
            If vars is an iterable, returns a list of names for the specified variables.
        """
        if isinstance(vars, Mapping):
            return { key: v.name for key, v in vars.items() }
        else:
            return [v.name for v in vars]

    def allVariableNames(self):
        """Retrieves the names of all variables in the model.

        Returns:
            A list of strings representing the names of all variables.
        """
        return super().getLp().col_names_

    def variableValue(self, var):
        """Retrieves the value of a specific variable in the solution.

        Args:
            var: A highs_var object representing the variable.

        Returns:
            The value of the specified variable in the solution.
        """
        return self.val(var)

    def variableValues(self, vars):
        """Retrieves the values of multiple variables in the solution.

        Args:
            vars: A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.

        Returns:
            If vars is a Mapping, returns a dict where keys are the same keys from the input vars and values are the solution values of the corresponding variables. If vars is an iterable, returns a list of solution values for the variables.
        """
        return self.vals(vars)


    def allVariableValues(self):
        """Retrieves the values of all variables in the solution.

        Returns:
            A list of values for all variables in the solution.
        """
        return super().getSolution().col_value

    def variableDual(self, var):
        """Retrieves the dual value of a specific variable/index or expression in the solution.

        Args:
            var: A highs_var object representing the variable.

        Returns:
            The dual value of the specified variable in the solution.
        """
        return self.internal_get_value(super().getSolution().col_dual, var)

    def variableDuals(self, vars):
        """Retrieves the dual values of multiple variables in the solution.

        Args:
            vars: A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.

        Returns:
            If vars is a Mapping, returns a dict where keys are the same keys from the input vars and values are the dual values of the corresponding variables. If vars is an iterable, returns a list of dual values for the variables.
        """
        return self.internal_get_value(super().getSolution().col_dual, vars)


    def allVariableDuals(self):
        """Retrieves the dual values of all variables in the solution.

        Returns:
            A list of dual values for all variables in the solution.
        """
        return super().getSolution().col_dual

    def constrValue(self, con):
        """Retrieves the value of a specific constraint in the solution.

        Args:
            con: A highs_con object representing the constraint.

        Returns:
            The value of the specified constraint in the solution.
        """
        return self.internal_get_value(super().getSolution().row_value, con)

    def constrValues(self, cons):
        """Retrieves the values of multiple constraints in the solution.

        Args:
            cons: A collection of highs_con objects representing the constraints. Can be a Mapping (e.g., dict) where keys are constraint names and values are highs_con objects, or an iterable of highs_con objects.

        Returns:
            If cons is a Mapping, returns a dict where keys are the same keys from the input cons and values are the solution values of the corresponding constraints. If cons is an iterable, returns a list of solution values for the constraints.
        """
        return self.internal_get_value(super().getSolution().row_value, cons)

    def allConstrValues(self):
        """Retrieves the values of all constraints in the solution.

        Returns:
            A list of values for all constraints in the solution.
        """
        return super().getSolution().row_value
    
    def constrDual(self, con):
        """Retrieves the dual value of a specific constraint in the solution.

        Args:
            con: A highs_con object representing the constraint.

        Returns:
            The dual value of the specified constraint in the solution.
        """
        return self.internal_get_value(super().getSolution().row_dual, con)
    
    def constrDuals(self, cons):
        """Retrieves the dual values of multiple constraints in the solution.

        Args:
            cons: A collection of highs_con objects representing the constraints. Can be a Mapping (e.g., dict) where keys are constraint names and values are highs_con objects, or an iterable of highs_con objects.

        Returns:
            If cons is a Mapping, returns a dict where keys are the same keys from the input cons and values are the dual values of the corresponding constraints. If cons is an iterable, returns a list of dual values for the constraints.
        """
        return self.internal_get_value(super().getSolution().row_dual, cons)

    def allConstrDuals(self):
        """Retrieves the dual values of all constraints in the solution.

        Returns:
            A list of dual values for all constraints in the solution.
        """
        return super().getSolution().row_dual

    def addVariable(self, lb = 0, ub = kHighsInf, obj = 0, type=HighsVarType.kContinuous, name = None):
        """Adds a variable to the model.

        Args:
            lb: Lower bound of the variable (default is 0).
            ub: Upper bound of the variable (default is infinity).
            obj: Objective coefficient of the variable (default is 0).
            type: Type of the variable (continuous, integer; default is continuous).
            name: Optional name for the variable.

        Returns:
            A highs_var object representing the added variable.
        """
        status = super().addCol(obj, lb, ub, 0, [], [])
        
        if status != HighsStatus.kOk:
            raise Exception("Failed to add variable to the model.")
        
        var = highs_var(self.numVariables - 1, self)

        if type != HighsVarType.kContinuous:
            super().changeColIntegrality(var.index, type)
            
        if name != None:
            super().passColName(var.index, name)
   
        return var

    def addVariables(self, *nvars, **kwargs):
        """Adds multiple variables to the model.

        Args:
            *args: A sequence of variables to be added. Can be a collection of scalars or indices (or mix).
            
            **kwargs: Optional keyword arguments.  Can be scalars, arrays, or mappings.
                lb: Lower bound of the variables (default is 0).
                ub: Upper bound of the variables (default is infinity).
                obj: Objective coefficient of the variables (default is 0).
                type: Type of the variables (continuous, integer; default is continuous).
                name: A collection of names for the variables (list or mapping).
                name_prefix: Prefix for the variable names.  Constructed name will be name_prefix + index.
                out_array: Return an array of highs_var objects instead of a dictionary.

        Returns:
            A highs_var collection (array or dictionary) representing the added variables.
        """        
        if len(nvars) == 0:
            return None

        lb = kwargs.get('lb', 0)
        ub = kwargs.get('ub', kHighsInf)
        obj = kwargs.get('obj', 0)
        vartype = kwargs.get('type', HighsVarType.kContinuous)
        name_prefix = kwargs.get('name_prefix', None)
        name = kwargs.get('name', None)
        out_array = kwargs.get('out_array', all([isinstance(n, int) for n in nvars]))

        shape = list(map(int, nvars)) if all([isinstance(n, int) for n in nvars]) else 1
        nvars = [range(n) if isinstance(n, int) else n for n in nvars]
        indices = list(nvars[0] if len(nvars) == 1 else product(*nvars))  # unpack tuple if needed
        N = len(indices)

        # parameter can be scalar, array, or mapping lookup (i.e., dictionary, custom class, etc.)
        # scalar: repeat for all N, array: use as is, lookup: convert to array using indices
        Rf = lambda x: np.fromiter((x[i] for i in indices), np.float64) if isinstance(x, Mapping) else (x if hasattr(x, "__getitem__") else np.full(N, x, dtype=np.float64))
        Ri = lambda x: np.fromiter((x[i] for i in indices), np.int32) if isinstance(x, Mapping) else (x if hasattr(x, "__getitem__") else np.full(N, x, dtype=np.int32))

        start_idx = self.numVariables        
        idx = np.arange(start_idx, start_idx + N, dtype=np.int32)
        status = super().addCols(N, Rf(obj), Rf(lb), Rf(ub), 0, [], [], [])
        
        if status != HighsStatus.kOk:
            raise Exception("Failed to add columns to the model.")            

        # only set integrality if we have non-continuous variables
        if vartype != HighsVarType.kContinuous:
            super().changeColsIntegrality(N, idx, Ri(vartype))
               
        if name or name_prefix:
            names = name or [f"{name_prefix}{i}" for i in indices]

            for i,n in zip(idx, names):
                super().passColName(int(i), str(n))

        return np.asarray([highs_var(i, self) for i in idx], dtype=highs_var).reshape(shape) if out_array == True else {index: highs_var(i, self) for index,i in zip(indices, idx)}

    def addIntegrals(self, *nvars, **kwargs):
        """
        Alias for the addVariables method, for integer variables.
        """
        kwargs.setdefault('type', HighsVarType.kInteger)
        return self.addVariables(*nvars, **kwargs)

    def addBinaries(self, *nvars, **kwargs):
        """
        Alias for the addVariables method, for binary variables.
        """
        kwargs.setdefault('lb', 0)
        kwargs.setdefault('ub', 1)
        kwargs.setdefault('type', HighsVarType.kInteger)

        return self.addVariables(*nvars, **kwargs)

    def addIntegral(self, lb = 0, ub = kHighsInf, obj = 0, name = None):
        """
        Alias for the addVariable method, for integer variables.
        """
        return self.addVariable(lb, ub, obj, HighsVarType.kInteger, name)

    def addBinary(self, obj = 0, name = None):
        """
        Alias for the addVariable method, for binary variables.
        """
        return self.addVariable(0, 1, obj, HighsVarType.kInteger, name)

    def deleteVariable(self, var_or_index, *args):
        """Deletes a variable from the model and updates the indices of subsequent variables in provided collections.

        Args:
            var_or_index: A highs_var object or an index representing the variable to be deleted.
            *args: Optional collections (lists, dicts, etc.) of highs_var objects whose indices need to be updated.
        """
        # Determine the index of the variable to delete
        index = int(var_or_index)

        # Delete the variable from the model if it exists
        if index < self.numVariables:
            super().deleteVars(1, [index])

        # Update the indices of variables in the provided collections
        for collection in args:
            if isinstance(collection, dict):
                # Update indices in a dictionary of variables
                for key, var in collection.items():
                    if var.index > index:
                        var.index -= 1
                    elif var.index == index:
                        var.index = -1

            elif hasattr(collection, '__iter__'):
                # Update indices in an iterable of variables
                for var in collection:
                    if var.index > index:
                        var.index -= 1
                    elif var.index == index:
                        var.index = -1

            # If the collection is a single highs_var object, check and update if necessary
            elif isinstance(collection, highs_var) and collection.index > index:
                collection.index -= 1
            elif isinstance(collection, highs_var) and collection.index == index:
                collection.index = -1

    def getVariables(self):
        """Retrieves all variables in the model.

        Returns:
            A list of highs_var objects, each representing a variable in the model.
        """
        return [highs_var(i, self) for i in range(self.numVariables)]

    @property
    def inf(self):
        """Represents infinity in the context of the solver.

        Returns:
            The value used to represent infinity.
        """
        return kHighsInf

    @property
    def numVariables(self):
        """Gets the number of variables in the model.

        Returns:
            The number of variables.
        """
        return super().getNumCol()

    @property
    def numConstrs(self):
        """Gets the number of constraints in the model.

        Returns:
            The number of constraints.
        """
        return super().getNumRow()

    #
    # add constraints
    #
    def addConstr(self, expr, name=None):
        """Adds a constraint to the model.

        Args:
            expr: A highs_linear_expression to be added.            
            name: Optional name of the constraint.

        Returns:
            A highs_con object representing the added constraint.
        """  
        con = self.__addRow(expr, self.numConstrs)

        if name != None:
            super().passRowName(con.index, name)

        return con

    def addConstrs(self, *args, **kwargs):
        """Adds multiple constraints to the model.

        Args:
            *args: A sequence of highs_linear_expression to be added.
            
            **kwargs: Optional keyword arguments.
                name_prefix: Prefix for the constraint names.  Constructed name will be name_prefix + index.
                name: A collection of names for the constraints (list or mapping).

        Returns:
            A highs_con collection array representing the added constraints.
        """  
        name_prefix = kwargs.get('name_prefix', None)
        name = kwargs.get('name', None)
        generator = args

        # unpack generator if needed
        if len(args) == 1 and hasattr(args[0], "__iter__") == True:
            generator = args[0]

        initial_rows = self.numConstrs

        try:
            if isinstance(generator, Mapping) == False:
                cons = [self.__addRow(expr, initial_rows + count) for count, expr in enumerate(generator)]
            else:
                cons = {key: self.__addRow(expr, initial_rows + count) for count, (key, expr) in enumerate(generator.items())}

            # TODO: Mapping support with constraing names can be improved, e.g., by allowing a name collection to be passed
            if name or name_prefix:
                names = name or [f"{name_prefix}{n}" for n in range(self.numConstrs - initial_rows)]

                for c,n in zip(range(initial_rows, self.numConstrs), names):
                    super().passRowName(int(c), str(n))

        except Exception as e:
            # rollback model if error - remove any constraints that were added
            status = super().deleteRows(self.numConstrs - initial_rows, np.arange(initial_rows, self.numConstrs, dtype=np.int32))

            if status != HighsStatus.kOk:
                raise Exception("Failed to rollback model after failure.  Model might be in a undeterminate state.")
            else:
                raise e

        return cons

    def __addRow(self, expr, idx):
        """
        Internal method to add a constraint to the model.
        """
        if expr.bounds != None:
            vars,vals = expr.unique_elements() # if we have duplicate variables, add the vals
            super().addRow(expr.bounds[0], expr.bounds[1], len(vars), vars, vals)
            return highs_cons(idx, self)
        else:
            raise Exception("Constraint bounds must be set via comparison (>=,==,<=).")

    def expr(self, optional=None):
        """Creates a new highs_linear_expression object.

        Returns:
            A highs_linear_expression object.
        """
        return highs_linear_expression(optional)

    def getExpr(self, cons):
        """Retrieves the highs_linear_expression of a constraint.

        Args:
            cons: A highs_con object or index representing the constraint.

        Returns:
            A highs_linear_expression object representing the expression of the constraint.
        """
        status, lb, ub, nnz = super().getRow(int(cons))

        if status != HighsStatus.kOk:
            raise Exception("Error retrieving constraint expression.")

        status, idx, val = super().getRowEntries(int(cons))

        if status != HighsStatus.kOk:
            raise Exception("Error retrieving constraint expression entries.")

        expr = highs_linear_expression()
        expr.bounds = [lb, ub]
        expr.vars = list(idx)
        expr.vals = list(val)
        return expr

    def chgCoeff(self, cons, var, val):
        """Changes the coefficient of a variable in a constraint.

        Args:
            cons: A highs_con object representing the constraint.
            var: A highs_var object representing the variable.
            val: The new coefficient value for the variable in the constraint.
        """
        super().changeCoeff(int(cons), int(var), val)

    def getConstrs(self):
        """Retrieves all constraints in the model.

        Returns:
            A list of highs_cons objects, each representing a constraint in the model.
        """
        return [highs_cons(i, self) for i in range(self.numConstrs)]

    def removeConstr(self, cons_or_index, *args):
        """Removes a constraint from the model and updates the indices of subsequent constraints in provided collections.

        Args:
            cons_or_index: A highs_cons object or an index representing the constraint to be removed.
            *args: Optional collections (lists, dicts, etc.) of highs_cons objects whose indices need to be updated after the removal.
        """
        # Determine the index of the constraint to delete
        index = int(cons_or_index)

        # Delete the variable from the model if it exists
        if index < self.numConstrs:
           status = super().deleteRows(1, [index])

           if status != HighsStatus.kOk:
               raise Exception("Failed to delete constraint from the model.")

        # Update the indices of constraints in the provided collections
        for collection in args:
            if isinstance(collection, dict):
                # Update indices in a dictionary of constraints
                for key, con in collection.items():
                    if con.index > index:
                        con.index -= 1
                    elif con.index == index:
                        con.index = -1

            elif hasattr(collection, '__iter__'):
                # Update indices in an iterable of constraints
                for con in collection:
                    if con.index > index:
                        con.index -= 1
                    elif con.index == index:
                        con.index = -1
            # If the collection is a single highs_cons object, check and update if necessary
            elif isinstance(collection, highs_cons) and collection.index > index:
                collection.index -= 1
            elif isinstance(collection, highs_cons) and collection.index == index:
                collection.index = -1

    def setMinimize(self):
        """
        Sets the objective sense of the model to minimization.
        """
        super().changeObjectiveSense(ObjSense.kMinimize)

    def setMaximize(self):
        """
        Sets the objective sense of the model to maximization.
        """
        super().changeObjectiveSense(ObjSense.kMaximize)

    def setInteger(self, var_or_collection):
        """Sets a variable/collection to integer.

        Args:
            var_or_collection: A highs_var object/collection representing the variable to be set as integer.
        """
        if hasattr(var_or_collection, '__iter__'):
            idx = np.fromiter(map(int, var_or_collection), dtype=np.int32)
            super().changeColsIntegrality(len(idx), idx, np.full(len(idx), HighsVarType.kInteger, dtype=np.int32))
        else:
            super().changeColIntegrality(int(var_or_collection), HighsVarType.kInteger)

    def setContinuous(self, var_or_collection):
        """Sets a variable/collection to continuous.

        Args:
            var_or_collection: A highs_var object/collection representing the variable to be set as continuous.
        """
        if hasattr(var_or_collection, '__iter__'):
            idx = np.fromiter(map(int, var_or_collection), dtype=np.int32)
            super().changeColsIntegrality(len(idx), idx, np.full(len(idx), HighsVarType.kContinuous, dtype=np.int32))
        else:
            super().changeColIntegrality(int(var_or_collection), HighsVarType.kContinuous)

    @staticmethod
    def qsum(items, initial=None):
        """Performs a faster sum for highs_linear_expressions.
        
        Args:
            items: A collection of highs_linear_expressions or highs_vars to be summed.
        """
        expr = highs_linear_expression(initial)

        for v in items:
            expr += v

        return expr


## The following classes keep track of variables
## It is currently quite basic and may fail in complex scenarios

# highs variable
class highs_var(object):
    """
    Basic constraint builder for HiGHS
    """
    __slots__ = ['index', 'highs']

    def __init__(self, i, highs):
        self.index = i
        self.highs = highs

    def __repr__(self):
        return f"highs_var({self.index})"

    @property
    def name(self):
        status, name = self.highs.getColName(self.index)
        
        if status != HighsStatus.kOk:
            raise Exception("Error retrieving variable name.")
        return name

    @name.setter
    def name(self, value):
        self.highs.passColName(self.index, value)

    def __int__(self):
        return int(self.index)

    def __hash__(self):
        return int(self.index)

    def __neg__(self):
        expr = highs_linear_expression(self)
        expr *= -1.0
        return expr
    
    def __add__(self, other):
        expr = highs_linear_expression(self)
        expr += other
        return expr

    def __radd__(self, other):
        expr = highs_linear_expression(self)
        expr += other
        return expr

    def __mul__(self, other):
        expr = highs_linear_expression(self)
        expr *= other
        return expr

    def __rmul__(self, other):
        expr = highs_linear_expression(self)
        expr *= other
        return expr

    def __rsub__(self, other):
        expr = highs_linear_expression(other)
        expr -= self
        return expr

    def __sub__(self, other):
        expr = highs_linear_expression(self)
        expr -= other
        return expr

    # self <= other
    def __le__(self, other):
        if isinstance(other, highs_linear_expression):
            return other.__ge__(self)
        else:
            return highs_linear_expression(self).__le__(other)

    # self == other
    def __eq__(self, other):
        if isinstance(other, highs_linear_expression):
            return other.__eq__(self)
        else:
            return highs_linear_expression(self).__eq__(other)

    # self != other
    def __ne__(self, other):
        if other == None:
            return True
        else:
            raise Exception('Invalid comparison.')

    # self >= other
    def __ge__(self, other):
        if isinstance(other, highs_linear_expression):
            return other.__le__(self)
        else:
            return highs_linear_expression(self).__ge__(other)

# highs constraint
class highs_cons(object):
    """
    Basic constraint for HiGHS
    """
    __slots__ = ['index', 'highs']
    
    def __init__(self, i, highs):
        self.index = i
        self.highs = highs

    def __repr__(self):
        return f"highs_cons({self.index})"

    def __int__(self):
        return int(self.index)

    def __hash__(self):
        return int(self.index)

    def expr(self):
        """
        Retrieves the expression of the constraint.
        """
        return self.highs.getExpr(self)

    @property
    def name(self):
        status, name = self.highs.getRowName(self.index)
        
        if status != HighsStatus.kOk:
            raise Exception("Error retrieving constraint name.")
        return name

    @name.setter
    def name(self, value):
        self.highs.passRowName(self.index, value)
   

# highs constraint builder
# 
# Note: we only allow LHS/RHS to be set once via comparisons (>=,==,<=), otherwise it gets confusing!
# e.g, for (0 <= x <= 1) <= 2, should this be:
#     0 <= x <= 1 (i.e., tighter bound, x <= 1 and x <= 2 implies x <= 1),
# or, 0 <= x <= 2 (i.e., last comparision overrides previous)?
#
# Throwing an error makes it obvious to the user. Note we can still set via addition,
# e.g., (x <= 1) + (y <= 5) to get x + y <= 6
#
# 
# For chained comparisons, we throw an error if we have mismatched variables in the "bounds",
# e.g., x <= y <= 1.  This requires 2 constraints, x <= 1 and x <= y, and is not supported, whereas
#   x + 2 <= y <= x + 5 is supported, since we can easily rewrite this as 2 <= y - x <= 5.
#
# As such, trivial chaining is supported, e.g., lb <= expr <= ub, where lb and ub are numeric.
#
# Note: when dealing with inequalities (>=,==,<=) we need to decide which variables to move the LHS/RHS.
# For consistency, we:
# 1. Move variables to side with the most variables, e.g., x <= y + z  ->     0 <= y + z - x <= inf
#                                                      y + z >= x      ->     0 <= y + z - x <= inf
#                                                          x >= y + z  ->  -inf <= y + z - x <= 0
#                                                      y + z <= x      ->  -inf <= y + z - x <= 0 
#                                                      y + z == x      ->     y + z - x == 0
#                                                          x == y + z  ->     y + z - x == 0
#
# 2. If equal, move to side without a constant,      e.g., y <= x + 2  ->  -inf <= y - x <= 2
#                                                      x + 2 >= y      ->  -inf <= y - x <= 2
#                                                      x + 2 <= y      ->     2 <= y - x <= inf
#                                                          y >= x + 2  ->     2 <= y - x <= inf
#                                                      y + 2 == x      ->     x - y == 2
#                                                          x == y + 2  ->     x - y == 2
#
# 3. If both have constants, move to the 'left',     e.g., x + 2 <= y + 3  ->  -inf <= x - y <=  1
#                                                          y + 3 >= x + 2  ->  -inf <= x - y <=  1
#                                                          x + 2 >= y + 3  ->  -inf <= y - x <= -1
#                                                          y + 3 <= x + 2  ->  -inf <= y - x <= -1
#                                                          x + 2 == y + 3  -> x - y ==  1
#                                                          y + 3 == x + 2  -> y - x == -1
# Hence:
#           6 <= x + y <= 8      ->  6 <= x + y     <= 8  #(1)
#       6 + x <= y     <= 8 + x  ->  6 <= y - x     <= 8  #(2)
#       6 + x <= y + 2 <= 8 + x  ->  4 <= y - x     <= 6  #(3)
#           x <= 6     <= x      ->  6 <= x         <= 6  #(1)
#           x <= y + z <= x + 5  ->  0 <= y + z - x <= 5  #(1)
class highs_linear_expression(object):
    """
    Basic constraint builder for HiGHS
    """
    __slots__ = ['vars', 'vals', 'bounds', 'constant']

    def __init__(self, other=None):
        self.constant = None  # constant is only valid when bounds are None
        self.bounds = None    # bounds are only valid when constant is None

        if other is None:
            self.vars = []
            self.vals = []

        elif isinstance(other, highs_linear_expression):
            self.vars = list(other.vars)
            self.vals = list(other.vals)
            self.constant = other.constant
            self.bounds = list(other.bounds) if other.bounds != None else None

        elif isinstance(other, highs_var):
            self.vars = [other.index]
            self.vals = [1.0]
        
        elif isinstance(other, (int, float, Decimal)):
            self.vars = []
            self.vals = []
            self.constant = other

        else:
            raise Exception('Invalid type for highs_linear_expression')

    def simplify(self):
        """
        Simplifies the linear expression by combining duplicate variables.
        """
        copy = highs_linear_expression()
        copy.vars, copy.vals = (v.tolist() for v in self.unique_elements())
        copy.bounds = list(self.bounds) if self.bounds != None else None
        copy.constant = self.constant
        return copy

    def copy(self):
        """
        Creates a copy of the linear expression.
        """
        return highs_linear_expression(self)

    def evaluate(self, values):
        """
        Evaluates the linear expression given a solution array (values).
        """
        result = sum(v * values[c] for c,v in zip(self.vars, self.vals)) + (self.constant or 0.0)
        return result if self.bounds == None else (self.bounds[0] <= result <= self.bounds[1])

    def __repr__(self):
        # display duplicate variables
        v = str.join("  ", [f"{c}_v{x}" for x,c in zip(self.vars,self.vals)])

        if self.bounds == None:
            return f"{v}" + (f"  {self.constant}" if self.constant != None else '')
        elif self.bounds[0] == self.bounds[1]:
            return f"{v} == {self.bounds[0] - (self.constant or 0.0)}"
        else:
            return f"{self.bounds[0]} <= {v} <= {self.bounds[1]}"

    def __str__(self):
        # display unique variables (values are totaled)
        vars,vals = self.unique_elements()
        v = str.join("  ", [f"{c}_v{x}" for x,c in zip(vars,vals)])

        if self.bounds == None:
            return f"{v}" + (f"  {self.constant}" if self.constant != None else '')
        elif self.bounds[0] == self.bounds[1]:
            return f"{v} == {self.bounds[0] - (self.constant or 0.0)}"
        else:
            return f"{self.bounds[0]} <= {v} <= {self.bounds[1]}"

    # self != other
    def __ne__(self, other):
        if other == None:
            return True
        else:
            raise Exception('Invalid comparison.')

    # self == other
    def __eq__(self, other):
        if self.bounds != None:
            raise Exception('Bounds have already been set.')

        # self == c
        elif isinstance(other, (int, float, Decimal)):
            copy = highs_linear_expression(self)
            copy.bounds = [float(other) - (self.constant or 0.0), float(other) - (self.constant or 0.0)]
            copy.constant = None
            return copy

        # self == other
        elif isinstance(other, highs_linear_expression):
            if other.bounds != None:
                raise Exception('Bounds have already been set.')

            copy = highs_linear_expression()

            # prefer most vars, constant, 'left'
            if len(other.vars) > len(self.vars) or (len(other.vars) == len(self.vars) and other.constant == None and self.constant != None):
                copy.vars = other.vars + self.vars
                copy.vals = other.vals + [-v for v in self.vals]
                copy.bounds = [(self.constant or 0.0) - (other.constant or 0.0), (self.constant or 0.0) - (other.constant or 0.0)]
            else:
                copy.vars = self.vars + other.vars
                copy.vals = self.vals + [-v for v in other.vals]
                copy.bounds = [(other.constant or 0.0) - (self.constant or 0.0), (other.constant or 0.0) - (self.constant or 0.0)]

            return copy

        # self == x
        elif isinstance(other, highs_var):
            copy = highs_linear_expression()

            if len(self.vars) == 0 or len(self.vars) == 1 and self.constant != None:
                copy.vars = [other.index] + self.vars
                copy.vals = [1.0] + [-v for v in self.vals]
                copy.bounds = [(self.constant or 0.0), (self.constant or 0.0)]
            else:
                copy.vars = self.vars + [other.index]
                copy.vals = self.vals + [-1.0]
                copy.bounds = [-(self.constant or 0.0), -(self.constant or 0.0)]

            return copy

        # support expr == [lb, ub] --> lb <= expr <= ub
        elif hasattr(other, "__getitem__") and hasattr(other, "__len__") and len(other) == 2:
            if not (isinstance(other[0], (int, float, Decimal)) and isinstance(other[1], (int, float, Decimal))):
                raise Exception('Provided bounds were not valid numbers.')

            copy = highs_linear_expression(self)
            copy.bounds = [float(other[0]) - (copy.constant or 0.0), float(other[1]) - (copy.constant or 0.0)]
            copy.constant = None
            return copy

        else:
            raise Exception('Unknown comparison.')

    # self <= other
    def __le__(self, other):
        if self.bounds != None:
            raise Exception('Bounds have already been set.')

        elif self.__is_active_chain():
            other = other if isinstance(other, highs_linear_expression) else highs_linear_expression(other)
            order = self.__get_chain(other, False)  # [self, LHS, inner, RHS, other], ignores None values

            # inner <= (self == RHS) <= other
            # LHS <= (self == inner) <= other
            if self.__is_equal_except_bounds(order[2]) == True:
                return self.__compose_chain(*order[1:])

            # self <= (other == inner) <= RHS
            # self <= (other == LHS) <= inner
            elif other.__is_equal_except_bounds(order[1]):
                return self.__compose_chain(*order[:-1])

        self.__reset_chain(self, other, None)
            
        # self <= c
        if isinstance(other, (int, float, Decimal)):
            copy = highs_linear_expression(self)
            copy.bounds = [-kHighsInf, float(other) - (copy.constant or 0.0)]
            copy.constant = None
            return copy

        # self <= other
        elif isinstance(other, highs_linear_expression):
            if other.bounds == None:
                copy = highs_linear_expression()

                # prefer most vars, constant, 'left'
                if len(other.vars) > len(self.vars) or (len(other.vars) == len(self.vars) and other.constant == None and self.constant != None):
                    copy.vars = other.vars + self.vars
                    copy.vals = other.vals + [-v for v in self.vals]
                    copy.bounds = [(self.constant or 0.0) - (other.constant or 0.0), kHighsInf]
                else:
                    copy.vars = self.vars + other.vars
                    copy.vals = self.vals + [-v for v in other.vals]
                    copy.bounds = [-kHighsInf, (other.constant or 0.0) - (self.constant or 0.0)]

                return copy
            else:
                raise Exception('Bounds have already been set.')

        # self <= x
        elif isinstance(other, highs_var):
            copy = highs_linear_expression()

            if len(self.vars) == 0 or len(self.vars) == 1 and self.constant != None:
                copy.vars = [other.index] + self.vars
                copy.vals = [1.0] + [-v for v in self.vals]
                copy.bounds = [(self.constant or 0.0), kHighsInf]
            else:
                copy.vars = self.vars + [other.index]
                copy.vals = self.vals + [-1.0]
                copy.bounds = [-kHighsInf, -(self.constant or 0.0)]

            return copy

        else:
            raise Exception('Unknown comparison.')

    # other <= self
    def __ge__(self, other):
        if self.bounds != None:
            raise Exception('Bounds have already been set.')

        elif self.__is_active_chain():
            other = other if isinstance(other, highs_linear_expression) else highs_linear_expression(other)
            order = self.__get_chain(other, True) # [other, LHS, inner, RHS, self], ignores None values

            # other <= (self == LHS) <= inner
            # other <= (self == inner) <= RHS
            if self.__is_equal_except_bounds(order[1]) == True:
                return self.__compose_chain(*order[:-1])

            # LHS <= (other == inner) <= self
            # inner <= (other == RHS) <= self
            elif other.__is_equal_except_bounds(order[2]):
                return self.__compose_chain(*order[1:])

        self.__reset_chain(None, other, self)

        # c <= self
        if isinstance(other, (int, float, Decimal)):
            copy = highs_linear_expression(self)
            copy.bounds = [float(other) - (self.constant or 0.0), kHighsInf]
            copy.constant = None
            return copy

        # other <= self
        elif isinstance(other, highs_linear_expression):
            if other.bounds == None:
                copy = highs_linear_expression()

                # prefer most vars, constant, 'left'
                if len(self.vars) > len(other.vars) or (len(self.vars) == len(other.vars) and self.constant == None and other.constant != None):
                    copy.vars = self.vars + other.vars
                    copy.vals = self.vals + [-v for v in other.vals]
                    copy.bounds = [(other.constant or 0.0) - (self.constant or 0.0), kHighsInf]
                else:
                    copy.vars = other.vars + self.vars
                    copy.vals = other.vals + [-v for v in self.vals]
                    copy.bounds = [-kHighsInf, (self.constant or 0.0) - (other.constant or 0.0)]

                return copy
            else:
                raise Exception('Bounds have already been set.')

        # x <= self
        elif isinstance(other, highs_var):
            copy = highs_linear_expression()

            if len(self.vars) > 1:
                copy.vars = self.vars + [other.index]
                copy.vals = self.vals + [-1.0]
                copy.bounds = [-(self.constant or 0.0), kHighsInf]
            else:
                copy.vars = [other.index] + self.vars
                copy.vals = [1.0] + [-v for v in self.vals]
                copy.bounds = [-kHighsInf, (self.constant or 0.0)]

            return copy

        else:
            raise Exception('Unknown comparison.')

    def __radd__(self, other):
        return self + other

    # (LHS <= self <= RHS) + (LHS <= other <= RHS)
    def __add__(self, other):
        copy = highs_linear_expression(self)
        copy += other
        return copy

    def __neg__(self):
        return -1.0 * self

    def __rmul__(self, other):
        return self * other

    def __mul__(self, other):
        copy = highs_linear_expression(self)
        copy *= other
        return copy

    # other - self
    def __rsub__(self, other):
        copy = highs_linear_expression(other)
        copy -= self
        return copy

    def __sub__(self, other):
        copy = highs_linear_expression(self)
        copy -= other
        return copy

    def unique_elements(self):
        """
        Collects unique variables and sums their corresponding values.  Keeps all values (including zeros).
        """
        # sort by groups for fast unique
        groups = np.asarray(self.vars, dtype=np.int32)
        order = np.argsort(groups, kind='stable')
        groups = groups[order]

        # get unique groups
        index = np.ones(len(groups), dtype=bool)
        index[:-1] = groups[1:] != groups[:-1]

        if index.all():
            values = np.asarray(self.vals, dtype=np.float64)
            values = values[order]
            return groups, values
        else:
            values = np.array(self.vals)
            values = np.cumsum(values[order])
            values = values[index]
            groups = groups[index]

            # calculate the correct sum (diff of cumsum)
            values[1:] = values[1:] - values[:-1]
            return groups, values

    def reduced_elements(self):
        """
        Similar to unique_elements, except keeps only non-zero values
        """
        vx, vl = self.unique_elements()
        zx = np.nonzero(vl)
        return vx[zx], vl[zx]

    ##
    ## mutable functions
    ##

    # (LHS <= self <= RHS) += (LHS <= other <= RHS)
    def __iadd__(self, other):
        if isinstance(other, highs_var):
            self.vars.append(other.index)
            self.vals.append(1.0)
            return self

        elif isinstance(other, highs_linear_expression):
            if (self.constant != None and other.bounds != None or self.bounds != None and other.constant != None):
                raise Exception('''Cannot add a bounded constraint to a constraint with a constant, i.e., (lb <= expr1 <= ub) + (expr2 + c). 
                    Unsure of your intent. Did you want: lb + c <= expr1 + expr2 <= ub + c?  Try: (lb <= expr1 <= ub) + (expr2 == c) instead.''')

            self.vars.extend(other.vars)
            self.vals.extend(other.vals)

            if self.constant != None or other.constant != None:
                self.constant = (self.constant or 0.0) + (other.constant or 0.0)

            # (l1 <= expr1 <= u1) + (l2 <= expr2 <= u2)  -->  l1 + l2 <= expr1 + expr2 <= u1 + u2
            if self.bounds != None and other.bounds != None:
                self.bounds[0] += other.bounds[0]
                self.bounds[1] += other.bounds[1]

            # (expr1) + (lb <= expr2 <= ub)  -->  lb <= expr1 + expr2 <= ub
            elif self.bounds == None and other.bounds != None:
                self.bounds = list(other.bounds)

            return self

        elif isinstance(other, (int, float, Decimal)):
            if self.bounds != None:
                raise Exception('''Cannot add a constant to a bounded constraint, i.e., (lb <= expr <= ub) + c. 
                    Unsure of your intent. Did you want: lb + c <= expr <= ub + c?  Try: (lb <= expr <= ub) + (highs_linear_expression() == c) instead.''')

            self.constant = float(other) + (self.constant or 0.0)
            return self

        else:
            raise Exception('Unexpected parameters.')

    def __isub__(self, other):
        if isinstance(other, highs_var):
            self.vars.append(other.index)
            self.vals.append(-1.0)
            return self

        elif isinstance(other, highs_linear_expression):
            if (self.constant != None and other.bounds != None or self.bounds != None and other.constant != None):
                raise Exception('''Cannot subtract a bounded constraint to a constraint with a constant, i.e., (lb <= expr1 <= ub) - (expr2 + c). 
                    Unsure of your intent. Did you want: lb - c <= expr1 - expr2 <= ub - c?  Try: (lb <= expr1 <= ub) - (expr2 == c) instead.''')

            self.vars.extend(other.vars)
            self.vals.extend([-v for v in other.vals])

            if self.constant != None or other.constant != None:
                self.constant = (self.constant or 0.0) - (other.constant or 0.0)

            # (l1 <= expr1 <= u1) - (l2 <= expr2 <= u2)  -->  l1 - u2 <= expr1 - expr2 <= u1 - l2
            if self.bounds != None and other.bounds != None:
                self.bounds[0] -= other.bounds[1]
                self.bounds[1] -= other.bounds[0]

            # expr1 - (lb <= expr2 <= ub)  -->  -ub <= expr1 - expr2 <= -lb
            elif self.bounds == None and other.bounds != None:
                self.bounds = [-other.bounds[1], -other.bounds[0]]

            return self

        elif isinstance(other, (int, float, Decimal)):
            if self.bounds != None:
                raise Exception('''Cannot subtract a constant to a bounded constraint, i.e., (lb <= expr <= ub) - c. 
                    Unsure of your intent. Did you want: lb - c <= expr <= ub - c?  Try: (lb <= expr <= ub) - (highs_linear_expression() == c) instead.''')

            self.constant = (self.constant or 0.0) - other
            return self

        else:
            raise Exception('Unexpected parameters.')

    def __imul__(self, other):
        if isinstance(other, (int, float, Decimal)):
            scale = float(other)

        # other is a constant expression, so treat as a scalar
        elif isinstance(other, highs_linear_expression) and other.vars == [] and other.constant != None:
            scale = float(other.constant)
        else:
            scale = None

        if scale != None:
            self.vals = [scale * v for v in self.vals]

            if self.constant != None:
                self.constant *= scale

            if self.bounds != None:
                # negative scale reverses bounds
                if scale >= 0:
                    self.bounds = [scale * self.bounds[0], scale * self.bounds[1]]
                else:
                    self.bounds = [scale * self.bounds[1], scale * self.bounds[0]]

            return self

        elif isinstance(other, highs_var):
            raise Exception('Only linear expressions are allowed.')

        # self only has a constant, so treat as a scalar
        elif self.vars == [] and self.constant != None:
            scale = float(self.constant)

            self.vars = other.vars
            self.vals = [scale * v for v in other.vals]

            if other.constant != None:
                self.constant = other.constant * scale
            else:
                self.constant = None

            if other.bounds != None:
                # negative scale reverses bounds
                if scale >= 0:
                    self.bounds = [scale * other.bounds[0], scale * other.bounds[1]]
                else:
                    self.bounds = [scale * other.bounds[1], scale * other.bounds[0]]

            return self
        else:
            raise Exception('Unexpected parameters.')


    # The following is needed to support chained comparison, i.e., lb <= expr <= ub. This is interpreted 
    # as '__bool__(lb <= expr) and (expr <= ub)'; returning (expr <= ub), since __bool__(lb <= expr) == True.
    # 
    # We essentially want to "rewrite" this as '(lb <= expr) <= ub', while keeping the expr instance immutable.
    # As a slight hack, we can use a shared (thread local) object to keep track of the chain.
    #
    # Whenever we perform an inequality, we first check if the current expression ('self') is part of a chain. 
    # If it is, we copy the inner '__chain' expression rather than 'self'.
    # 
    # This inner '__chain' is set by __bool__(lb <= expr) and is reset after the inequality is evaluated.
    #
    # Two potential issues:
    #   1. It is possible to manually construct this sequence, e.g., 
    #         tmp = x0 + x1
    #         bool(tmp <= 10)
    #         print(5 <= tmp) # outputs: 5 <= 1.0_v0 + 1.0_v1 <= 10
    #         print(5 <= tmp) # outputs: 5 <= 1.0_v0 + 1.0_v1 <= inf : chain is broken
    #  
    #      Note that:
    #         bool(tmp <= 10)
    #         tmp += x0 + x1  # changes tmp, so the chain is broken
    #         print(5 <= tmp) # outputs: 5 <= 2.0_v0 + 2.0_v1 <= inf
    #  
    #   2. The chain might "break" if run within a debugger (on same thread), i.e., "watched debugger expressions" 
    #      that evaluate any variant of highs_linear_expression inequalities.
    #
    # I believe these issues are low risk, the approach thread safe, and the performance/overhead is minimal.
    #
    __chain = local()

    # capture the chain
    def __bool__(self):
        highs_linear_expression.__chain.check = self

        # take copies of the chain to avoid issues with mutable expressions
        LHS = getattr(highs_linear_expression.__chain, 'left', None)
        EXR = getattr(highs_linear_expression.__chain, 'inner', None)
        RHS = getattr(highs_linear_expression.__chain, 'right', None)

        highs_linear_expression.__chain.left  = highs_linear_expression(LHS) if LHS != None else None
        highs_linear_expression.__chain.inner = highs_linear_expression(EXR) if EXR != None else None
        highs_linear_expression.__chain.right = highs_linear_expression(RHS) if RHS != None else None
        return True

    def __is_equal_except_bounds(self, other):
        return self.vars == other.vars and self.vals == other.vals and self.constant == other.constant

    def __is_active_chain(self):
        return getattr(highs_linear_expression.__chain, 'check', None) is not None

    def __reset_chain(self, LHS=None, inner=None, RHS=None):
        highs_linear_expression.__chain.check = None
        highs_linear_expression.__chain.left = LHS
        highs_linear_expression.__chain.inner = inner
        highs_linear_expression.__chain.right = RHS

    def __get_chain(self, other, is_ge_than):
        LHS = getattr(highs_linear_expression.__chain, 'left', None)
        RHS = getattr(highs_linear_expression.__chain, 'right', None)
        inner = getattr(highs_linear_expression.__chain, 'inner', None)

        # assume (LHS is None) ^ (RHS is None) == 1, i.e., only LHS or RHS is set
        assert((LHS is None) ^ (RHS is None) == 1)

        order = np.asarray([other, LHS, inner, RHS, self])
        order = order[order != None]

        if is_ge_than == False:  # swap order for: self <= other
            order[0], order[-1] = order[-1], order[0]

        return order

    # left <= self <= right
    def __compose_chain(self, left, inner, right):
        self.__reset_chain()
        assert (isinstance(inner, highs_linear_expression) == False or inner.bounds == None), 'Bounds already set in chain comparison.'

        # check to see if we have a valid chain, i.e., left.vars "==" right.vars
        # we can assume that left and right are both highs_linear_expression
        LHS_vars, LHS_vals = left.reduced_elements()
        RHS_vars, RHS_vals = right.reduced_elements()
        
        if np.array_equal(LHS_vars, RHS_vars) == False or np.array_equal(LHS_vals, RHS_vals) == False:
            raise Exception('Mismatched variables in chain comparison.')

        if len(LHS_vars) > len(inner.vars):
            copy = highs_linear_expression()
            copy.vars = LHS_vars.tolist() + inner.vars
            copy.vals = LHS_vals.tolist() + [-v for v in inner.vals]
            copy.bounds = [(inner.constant or 0.0) - (right.constant or 0.0), (inner.constant or 0.0) - (left.constant or 0.0)]
        else:
            copy = highs_linear_expression(inner)
            copy.vars.extend(LHS_vars)
            copy.vals.extend([-v for v in LHS_vals])
            copy.bounds = [(left.constant or 0.0) - (copy.constant or 0.0), (right.constant or 0.0) - (copy.constant or 0.0)]
            copy.constant = None

        return copy

@staticmethod
def qsum(items, initial=None):
    """Performs a faster sum for highs_linear_expressions.
        
    Args:
        items: A collection of highs_linear_expressions or highs_vars to be summed.
    """
    return Highs.qsum(items, initial)
