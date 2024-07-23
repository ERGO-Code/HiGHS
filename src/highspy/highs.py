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
from itertools import groupby, product
from operator import itemgetter
from decimal import Decimal

class Highs(_Highs):
    """HiGHS solver interface"""

    def __init__(self):
        super().__init__()

    # Silence logging
    def silent(self):
        """Disables solver output to the console."""
        super().setOptionValue("output_flag", False)
    
    # solve
    def solve(self):
        """Runs the solver on the current problem.

        Returns:
            A HighsStatus object containing the solve status.
        """
        return super().run()

    def optimize(self):
        """Alias for the solve method."""
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
        if obj != None:
            # if we have a single variable, wrap it in a linear expression
            if isinstance(obj, highs_var) == True:
                obj = highs_linear_expression(obj)

            if isinstance(obj, highs_linear_expression) == False or obj.LHS != -self.inf or obj.RHS != self.inf:
                raise Exception('Objective cannot be an inequality') 

            # reset objective
            super().changeColsCost(self.numVariables, range(self.numVariables), [0]*self.numVariables)

            # if we have duplicate variables, add the vals
            vars,vals = zip(*[(var, sum(v[1] for v in Vals)) for var, Vals in groupby(sorted(zip(obj.vars, obj.vals)), key=itemgetter(0))])
            super().changeColsCost(len(vars), vars, vals)
            super().changeObjectiveOffset(obj.constant)

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
        if obj != None:
            # if we have a single variable, wrap it in a linear expression
            if isinstance(obj, highs_var) == True:
                obj = highs_linear_expression(obj)

            if isinstance(obj, highs_linear_expression) == False or obj.LHS != -self.inf or obj.RHS != self.inf:
                raise Exception('Objective cannot be an inequality') 

            # reset objective
            super().changeColsCost(self.numVariables, range(self.numVariables), [0]*self.numVariables)

            # if we have duplicate variables, add the vals
            vars,vals = zip(*[(var, sum(v[1] for v in Vals)) for var, Vals in groupby(sorted(zip(obj.vars, obj.vals)), key=itemgetter(0))])
            super().changeColsCost(len(vars), vars, vals)
            super().changeObjectiveOffset(obj.constant)

        super().changeObjectiveSense(ObjSense.kMaximize)
        return super().run()

    def internal_get_value(self, var_index_collection, col_value):
        """Internal method to get the value of a variable in the solution. Could be value or dual."""
        if isinstance(var_index_collection, int):
            return col_value[var_index_collection]
        elif isinstance(var_index_collection, highs_var):
            return col_value[var_index_collection.index]
        elif isinstance(var_index_collection, Mapping):
            return {k: col_value[v.index] for k,v in var_index_collection.items()}
        else:
            return [col_value[v.index] for v in var_index_collection]

    def val(self, var):
        """Gets the value of a variable in the solution.

        Args:
            var: A highs_var object representing the variable.

        Returns:
            The value of the variable in the solution.
        """
        return super().getSolution().col_value[var.index]

    def vals(self, vars):
        """Gets the values of multiple variables in the solution.

        Args:
            vars: A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.

        Returns:
            If vars is a Mapping, returns a dict where keys are the same keys from the input vars and values are the solution values of the corresponding variables. If vars is an iterable, returns a list of solution values for the variables.
        """
        col_value = super().getSolution().col_value
        return {k: self.internal_get_value(v, col_value) for k,v in vars.items()} if isinstance(vars, Mapping) else [self.internal_get_value(v, col_value) for v in vars]

    def variableName(self, var):
        """Retrieves the name of a specific variable.

        Args:
            var: A highs_var object representing the variable.

        Raises:
            Exception: If the variable name cannot be found.

        Returns:
            The name of the specified variable.
        """
        [status, name] = super().getColName(var.index)
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
        return super().getSolution().col_value[var.index]

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
        """Retrieves the dual value of a specific variable in the solution.

        Args:
            var: A highs_var object representing the variable.

        Returns:
            The dual value of the specified variable in the solution.
        """
        return super().getSolution().col_dual[var.index]

    def variableDuals(self, vars):
        """Retrieves the dual values of multiple variables in the solution.

        Args:
            vars: A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.

        Returns:
            If vars is a Mapping, returns a dict where keys are the same keys from the input vars and values are the dual values of the corresponding variables. If vars is an iterable, returns a list of dual values for the variables.
        """
        col_dual = super().getSolution()
        return {k: self.internal_get_value(v, col_dual) for k,v in vars.items()} if isinstance(vars, Mapping) else [self.internal_get_value(v, col_dual) for v in vars]


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
        return super().getSolution().row_value[con.index]

    def constrValues(self, cons):
        """Retrieves the values of multiple constraints in the solution.

        Args:
            cons: A collection of highs_con objects representing the constraints. Can be a Mapping (e.g., dict) where keys are constraint names and values are highs_con objects, or an iterable of highs_con objects.

        Returns:
            If cons is a Mapping, returns a dict where keys are the same keys from the input cons and values are the solution values of the corresponding constraints. If cons is an iterable, returns a list of solution values for the constraints.
        """
        row_value = super().getSolution().row_value
        return {k: row_value[c.index] for k,c in cons.items()} if isinstance(cons, Mapping) else [row_value[c.index] for c in cons]


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
        return super().getSolution().row_dual[con.index]
    
    def constrDuals(self, cons):
        """Retrieves the dual values of multiple constraints in the solution.

        Args:
            cons: A collection of highs_con objects representing the constraints. Can be a Mapping (e.g., dict) where keys are constraint names and values are highs_con objects, or an iterable of highs_con objects.

        Returns:
            If cons is a Mapping, returns a dict where keys are the same keys from the input cons and values are the dual values of the corresponding constraints. If cons is an iterable, returns a list of dual values for the constraints.
        """
        row_dual = super().getSolution().row_dual
        return {k: row_dual[c.index] for k,c in cons.items()} if isinstance(cons, Mapping) else [row_dual[c.index] for c in cons]

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
        out_array = kwargs.get('out_array', len(nvars) == 1 and isinstance(nvars[0], int))

        nvars = [range(n) if isinstance(n, int) else n for n in nvars]
        indices = list(nvars[0] if len(nvars) == 1 else product(*nvars))  # unpack tuple if needed
        N = len(indices)

        # parameter can be scalar, array, or mapping lookup (i.e., dictionary, custom class, etc.)
        # scalar: repeat for all N, array: use as is, lookup: convert to array using indices
        R = lambda x: [x[i] for i in indices] if isinstance(x, Mapping) else (x if hasattr(x, "__getitem__") else [x] * N)

        start_idx = self.numVariables        
        idx = range(start_idx, start_idx + N)
        status = super().addCols(N, R(obj), R(lb), R(ub), 0, [], [], [])
        
        if status != HighsStatus.kOk:
            raise Exception("Failed to add columns to the model.")            

        # only set integrality if we have non-continuous variables
        if vartype != HighsVarType.kContinuous:
            super().changeColsIntegrality(N, idx, R(vartype))
               
        if name or name_prefix:
            names = name or [f"{name_prefix}{i}" for i in indices]

            for i,n in zip(idx, names):
                super().passColName(int(i), str(n))

        return [highs_var(i, self) for i in idx] if out_array == True else {index: highs_var(i, self) for index,i in zip(indices, idx)}

    def addIntegrals(self, *nvars, **kwargs):
        """Alias for the addVariables method, for integer variables."""
        kwargs.setdefault('type', HighsVarType.kInteger)
        return self.addVariables(*nvars, **kwargs)

    def addBinaries(self, *nvars, **kwargs):
        """Alias for the addVariables method, for binary variables."""
        kwargs.setdefault('lb', 0)
        kwargs.setdefault('ub', 1)
        kwargs.setdefault('type', HighsVarType.kInteger)

        return self.addVariables(*nvars, **kwargs)

    def addIntegral(self, lb = 0, ub = kHighsInf, obj = 0, name = None):
        """Alias for the addVariable method, for integer variables."""
        return self.addVariable(lb, ub, obj, HighsVarType.kInteger, name)

    def addBinary(self, obj = 0, name = None):
        """Alias for the addVariable method, for binary variables."""
        return self.addVariable(0, 1, obj, HighsVarType.kInteger, name)

    def deleteVariable(self, var_or_index, *args):
        """Deletes a variable from the model and updates the indices of subsequent variables in provided collections.

        Args:
            var_or_index: A highs_var object or an index representing the variable to be deleted.
            *args: Optional collections (lists, dicts, etc.) of highs_var objects whose indices need to be updated.
        """
        # Determine the index of the variable to delete
        index = var_or_index.index if isinstance(var_or_index, highs_var) else var_or_index

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
            elif hasattr(collection, '__iter__'):
                # Update indices in an iterable of variables
                for var in collection:
                    if var.index > index:
                        var.index -= 1
            # If the collection is a single highs_var object, check and update if necessary
            elif isinstance(collection, highs_var) and collection.index > index:
                collection.index -= 1

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
    def addConstr(self, cons, name=None):
        """Adds a constraint to the model.

        Args:
            cons: A highs_linear_expression to be added.            
            name: Optional name of the constraint.

        Returns:
            A highs_con object representing the added constraint.
        """  
        # if we have duplicate variables, add the vals
        vars,vals = zip(*[(var, sum(v[1] for v in Vals)) for var, Vals in groupby(sorted(zip(cons.vars, cons.vals)), key=itemgetter(0))])
        super().addRow(cons.LHS - cons.constant, cons.RHS - cons.constant, len(vars), vars, vals)
        con = highs_cons(self.numConstrs - 1, self)

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

        lower = []
        upper = []
        starts = [0]
        indices = []
        values = []
        nnz = 0;

        for cons in generator:
            # if we have duplicate variables, add the vals together
            vars,vals = zip(*[(var, sum(v[1] for v in Vals)) for var, Vals in groupby(sorted(zip(cons.vars, cons.vals)), key=itemgetter(0))])
            
            indices.extend(vars)
            values.extend(vals)
            nnz += len(vars)

            lower.append(cons.LHS - cons.constant)
            upper.append(cons.RHS - cons.constant)
            starts.append(nnz)

        new_rows = len(lower)
        super().addRows(new_rows, lower, upper, nnz, starts, indices, values);
        cons = [highs_cons(self.numConstrs - new_rows + n, self) for n in range(new_rows)] 

        if name or name_prefix:
            names = name or [f"{name_prefix}{n}" for n in range(new_rows)]

            for c,n in zip(cons, names):
                super().passRowName(int(c.index), str(n))

        return cons


    def chgCoeff(self, cons, var, val):
        """Changes the coefficient of a variable in a constraint.

        Args:
            cons: A highs_con object representing the constraint.
            var: A highs_var object representing the variable.
            val: The new coefficient value for the variable in the constraint.
        """
        super().changeCoeff(cons.index, var.index, val)

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
        index = cons_or_index.index if isinstance(cons_or_index, highs_cons) else cons_or_index

        # Delete the variable from the model if it exists
        if index < self.numConstrs:
           super().deleteRows(1, [index])

        # Update the indices of constraints in the provided collections
        for collection in args:
            if isinstance(collection, dict):
                # Update indices in a dictionary of constraints
                for key, con in collection.items():
                    if con.index > index:
                        con.index -= 1
            elif hasattr(collection, '__iter__'):
                # Update indices in an iterable of constraints
                for con in collection:
                    if con.index > index:
                        con.index -= 1
            # If the collection is a single highs_cons object, check and update if necessary
            elif isinstance(collection, highs_cons) and collection.index > index:
                collection.index -= 1


    def setMinimize(self):
        """Sets the objective sense of the model to minimization."""
        super().changeObjectiveSense(ObjSense.kMinimize)

    def setMaximize(self):
        """Sets the objective sense of the model to maximization."""
        super().changeObjectiveSense(ObjSense.kMaximize)

    def setInteger(self, var):
        """Sets a variable's type to integer.

        Args:
            var: A highs_var object representing the variable to be set as integer.
        """
        super().changeColIntegrality(var.index, HighsVarType.kInteger)

    def setContinuous(self, var):
        """Sets a variable's type to continuous.

        Args:
            var: A highs_var object representing the variable to be set as continuous.
        """
        super().changeColIntegrality(var.index, HighsVarType.kContinuous)

## The following classes keep track of variables
## It is currently quite basic and may fail in complex scenarios

# highs variable
class highs_var(object):
    """Basic constraint builder for HiGHS"""
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

    def __hash__(self):
        return self.index

    def __neg__(self):
        return -1.0 * highs_linear_expression(self)
    
    def __le__(self, other):
        return highs_linear_expression(self) <= other

    def __eq__(self, other):
        return highs_linear_expression(self) == other
    
    def __ge__(self, other):
        return highs_linear_expression(self) >= other

    def __add__(self, other):
        return highs_linear_expression(self) + other

    def __radd__(self, other):
        return highs_linear_expression(self) + other

    def __mul__(self, other):
        return highs_linear_expression(self) * other

    def __rmul__(self, other):
        return highs_linear_expression(self) * other

    def __rsub__(self, other):
        return -1.0 * highs_linear_expression(self) + other

    def __sub__(self, other):
        return highs_linear_expression(self) - other

# highs constraint
class highs_cons(object):
    """Basic constraint for HiGHS"""
    __slots__ = ['index', 'highs']
    
    def __init__(self, i, highs):
        self.index = i
        self.highs = highs

    def __repr__(self):
        return f"highs_cons({self.index})"

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
class highs_linear_expression(object):
    """Basic constraint builder for HiGHS"""
    __slots__ = ['vars', 'vals', 'LHS', 'RHS', 'constant']

    def __init__(self, other=None):
        self.constant = 0
        self.LHS = -kHighsInf
        self.RHS = kHighsInf

        if isinstance(other, highs_linear_expression):
            self.vars = list(other.vars)
            self.vals = list(other.vals)
            self.constant = other.constant
            self.LHS = other.LHS
            self.RHS = other.RHS

        elif isinstance(other, highs_var):
            self.vars = [other.index]
            self.vals = [1.0]
        else:
            self.vars = []
            self.vals = []

    def __neg__(self):
        return -1.0 * self

    # (LHS <= self <= RHS) <= (other.LHS <= other <= other.RHS)
    def __le__(self, other):
        if isinstance(other, highs_linear_expression):
            if self.LHS != -kHighsInf and self.RHS != kHighsInf and len(other.vars) > 0 or other.LHS != -kHighsInf:
                raise Exception('Cannot construct constraint with variables as bounds.')

            # move variables from other to self
            self.vars.extend(other.vars)
            self.vals.extend([-1.0 * v for v in other.vals])
            self.constant -= other.constant
            self.RHS = 0
            return self

        elif isinstance(other, highs_var):
            return NotImplemented

        elif isinstance(other, (int, float, Decimal)):
            self.RHS = min(self.RHS, other)
            return self

        else:
            return NotImplemented

   # (LHS <= self <= RHS) == (other.LHS <= other <= other.RHS)
    def __eq__(self, other):
        if isinstance(other, highs_linear_expression):
            if self.LHS != -kHighsInf and len(other.vars) > 0 or other.LHS != -kHighsInf:
                raise Exception('Cannot construct constraint with variables as bounds.')

            # move variables from other to self
            self.vars.extend(other.vars)
            self.vals.extend([-1.0 * v for v in other.vals])
            self.constant -= other.constant
            self.LHS = 0
            self.RHS = 0
            return self

        elif isinstance(other, highs_var):
            return NotImplemented

        elif isinstance(other, (int, float, Decimal)):
            if self.LHS != -kHighsInf or self.RHS != kHighsInf:
                raise Exception('Logic error in constraint equality.')

            self.LHS = other
            self.RHS = other
            return self

        else:
            return NotImplemented

    # (other.LHS <= other <= other.RHS) <= (LHS <= self <= RHS)
    def __ge__(self, other):
        if isinstance(other, highs_linear_expression):
            return other <= self

        elif isinstance(other, highs_var):
            return NotImplemented

        elif isinstance(other, (int, float, Decimal)):
            self.LHS = max(self.LHS, other)
            return self

        else:
            return NotImplemented

    def __radd__(self, other):
        return self + other

    # (LHS <= self <= RHS) + (LHS <= other <= RHS)
    def __add__(self, other):
        if isinstance(other, highs_linear_expression):
            self.vars.extend(other.vars)
            self.vals.extend(other.vals)
            self.constant += other.constant
            self.LHS = max(self.LHS, other.LHS)
            self.RHS = min(self.RHS, other.RHS)
            return self

        elif isinstance(other, highs_var):
            self.vars.append(other.index)
            self.vals.append(1.0)
            return self

        elif isinstance(other, (int, float, Decimal)):
            self.constant += other
            return self

        else:
            return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __mul__(self, other):
        result = highs_linear_expression(self)

        if isinstance(other, (int, float, Decimal)):
            result.vals = [float(other) * v for v in self.vals]
            result.constant *= float(other)
            return result
        elif isinstance(other, highs_var):
            raise Exception('Only linear expressions are allowed.')
        else:
            return NotImplemented

    def __rsub__(self, other):
        return other + -1.0 * self

    def __sub__(self, other):
        if isinstance(other, highs_linear_expression):
            return self + (-1.0 * other)
        elif isinstance(other, highs_var):
            return self + (-1.0 * highs_linear_expression(other))
        elif isinstance(other, (int, float, Decimal)):
            return self + (-1.0 * other)
        else:
            return NotImplemented
