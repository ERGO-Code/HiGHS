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


from itertools import groupby
from operator import itemgetter
from decimal import Decimal

class Highs(_Highs):
    """HiGHS solver interface"""
    __slots__ = ['_batch', '_vars', '_cons']

    def __init__(self):
        super().__init__()
        
        self._batch = highs_batch(self)
        self._vars = []
        self._cons = []

    # Silence logging
    def silent(self):
        super().setOptionValue("output_flag", False)
    
    # solve
    def solve(self):
        return super().run()

    # reset the objective and sense, then solve
    def minimize(self, obj=None):
        if obj != None:
            # if we have a single variable, wrap it in a linear expression
            if isinstance(obj, highs_var) == True:
                obj = highs_linear_expression(obj)

            if isinstance(obj, highs_linear_expression) == False or obj.LHS != -self.inf or obj.RHS != self.inf:
                raise Exception('Objective cannot be an inequality') 

            # reset objective
            self.update()
            super().changeColsCost(self.numVars, range(self.numVars), [0]*self.numVars)

            # if we have duplicate variables, add the vals
            vars,vals = zip(*[(var, sum(v[1] for v in Vals)) for var, Vals in groupby(sorted(zip(obj.vars, obj.vals)), key=itemgetter(0))])
            super().changeColsCost(len(vars), vars, vals)
            super().changeObjectiveOffset(obj.constant)

        super().changeObjectiveSense(ObjSense.kMinimize)
        return super().run()

    # reset the objective and sense, then solve
    def maximize(self, obj=None):
        if obj != None:
            # if we have a single variable, wrap it in a linear expression
            if isinstance(obj, highs_var) == True:
                obj = highs_linear_expression(obj)

            if isinstance(obj, highs_linear_expression) == False or obj.LHS != -self.inf or obj.RHS != self.inf:
                raise Exception('Objective cannot be an inequality') 

            # reset objective
            self.update()
            super().changeColsCost(self.numVars, range(self.numVars), [0]*self.numVars)

            # if we have duplicate variables, add the vals
            vars,vals = zip(*[(var, sum(v[1] for v in Vals)) for var, Vals in groupby(sorted(zip(obj.vars, obj.vals)), key=itemgetter(0))])
            super().changeColsCost(len(vars), vars, vals)
            super().changeObjectiveOffset(obj.constant)

        super().changeObjectiveSense(ObjSense.kMaximize)
        return super().run()

    
    # update variables
    def update(self):
        current_batch_size = len(self._batch.obj)
        if current_batch_size > 0:
            names = [self._batch.name[i] for i in range(current_batch_size)]
            super().addVars(int(current_batch_size), self._batch.lb, self._batch.ub)
            super().changeColsCost(current_batch_size, self._batch.idx, self._batch.obj)

            # only set integrality if we have non-continuous variables
            if any([t != HighsVarType.kContinuous for t in self._batch.type]):
                super().changeColsIntegrality(current_batch_size, self._batch.idx, self._batch.type)

            for i in range(current_batch_size):
                super().passColName(int(self._batch.idx[i]), str(names[i]))
        self._batch = highs_batch(self)

    def val(self, var):
        return super().getSolution().col_value[var.index]

    def vals(self, vars):
        sol = super().getSolution()
        return [sol.col_value[v.index] for v in vars]

    def varName(self, var):
        [status, name] = super().getColName(var.index)
        failed = status != HighsStatus.kOk
        if failed:
            raise Exception('Variable name not found') 
        return name

    def varNames(self, vars):
        names = list()
        for v in vars:
            [status, name] = super().getColName(v.index)
            failed = status != HighsStatus.kOk
            if failed:
                raise Exception('Variable name not found')
            names.append(name)
        return names

    def allVarNames(self):
        return super().getLp().col_names_

    def varValue(self, var):
        return super().getSolution().col_value[var.index]

    def varValues(self, vars):
        col_value = super().getSolution().col_value
        return [col_value[v.index] for v in vars]

    def allVarValues(self):
        return super().getSolution().col_value

    def varDual(self, var):
        return super().getSolution().col_dual[var.index]

    def varDuals(self, vars):
        col_dual = super().getSolution()
        return [col_dual[v.index] for v in vars]

    def allVarDuals(self):
        return super().getSolution().col_dual

    def constrValue(self, constr_name):
        status_index = super().getRowByName(constr_name)
        failed = status_index[0] != HighsStatus.kOk
        if failed:
            raise Exception('Constraint name not found') 
        return super().getSolution().row_value[status_index[1]]

    def constrValues(self, constr_names):
        row_value = super().getSolution().row_value
        index = list()
        for name in constr_names:
            status_index = super().getRowByName(name)
            failed = status_index[0] != HighsStatus.kOk
            if failed:
                raise Exception('Constraint name not found')
            index.append(status_index[1])
        return [row_value[index[v]] for v in range(len(index))]

    def allConstrValues(self):
        return super().getSolution().row_value

    def constrDual(self, constr_name):
        status_index = super().getRowByName(constr_name)
        failed = status_index[0] != HighsStatus.kOk
        if failed:
            raise Exception('Constraint name not found') 
        return super().getSolution().row_dual[status_index[1]]

    def constrDuals(self, constr_names):
        row_dual = super().getSolution().row_dual
        index = list()
        for name in constr_names:
            status_index = super().getRowByName(name)
            failed = status_index[0] != HighsStatus.kOk
            if failed:
                raise Exception('Constraint name not found')
            index.append(status_index[1])
        return [row_dual[index[v]] for v in range(len(index))]

    def allConstrDuals(self):
        return super().getSolution().row_dual

    #
    # add variable & useful constants
    #
    # Change the name of addVar to addVariable to prevent shadowing of
    # highspy binding to Highs::addVar
    def addVariable(self, lb = 0, ub = kHighsInf, obj = 0, type=HighsVarType.kContinuous, name = None):
        var = self._batch.add(obj, lb, ub, type, name, self)
        self._vars.append(var)
        # No longer acumulate a batch of variables so that addVariable
        # behaves like Highs::addVar and highspy bindings modifying
        # column data and adding rows can be used
        self.update()
        return var

    def addIntegral(self, lb = 0, ub = kHighsInf, obj = 0, name = None):
        return self.addVariable(lb, ub, obj, HighsVarType.kInteger, name)

    def addBinary(self, obj = 0, name = None):
        return self.addVariable(0, 1, obj, HighsVarType.kInteger, name)

    # Change the name of removeVar to deleteVariable
    def deleteVariable(self, var):
        for i in self._vars[var.index+1:]:
            i.index -= 1

        del self._vars[var.index]
    
        # only delete from model if it exists
        if var.index < self.numVars:
            super().deleteVars(1, [var.index])

    # Change the name of getVars to getVariables
    def getVariables(self):
        return self._vars

    @property
    def inf(self):
        return kHighsInf

    @property
    def numVars(self):
        return super().getNumCol()

    @property
    def numConstrs(self):
        return super().getNumRow()

    #
    # add constraints
    #
    def addConstr(self, cons, name=None):
        self.update()

        # if we have duplicate variables, add the vals
        vars,vals = zip(*[(var, sum(v[1] for v in Vals)) for var, Vals in groupby(sorted(zip(cons.vars, cons.vals)), key=itemgetter(0))])
        super().addRow(cons.LHS - cons.constant, cons.RHS - cons.constant, len(vars), vars, vals)

        cons = highs_cons(self.numConstrs - 1, self, name)
        self._cons.append(cons)
        return cons

    def chgCoeff(self, cons, var, val):
        super().changeCoeff(cons.index, var.index, val)

    def getConstrs(self):
        return self._cons

    def removeConstr(self, cons):
        for i in self._cons[cons.index+1:]:
            i.index -= 1

        del self._cons[cons.index]
        super().deleteRows(1, [cons.index])

    # set to minimization
    def setMinimize(self):
        super().changeObjectiveSense(ObjSense.kMinimize)

    # set to maximization
    def setMaximize(self):
        super().changeObjectiveSense(ObjSense.kMaximize)

    # Set to integer
    def setInteger(self, var):
        super().changeColIntegrality(var.index, HighsVarType.kInteger)

    # Set to continuous
    def setContinuous(self, var):
        super().changeColIntegrality(var.index, HighsVarType.kContinuous)

## The following classes keep track of variables
## It is currently quite basic and may fail in complex scenarios

# highs variable
class highs_var(object):
    """Basic constraint builder for HiGHS"""
    __slots__ = ['index', '_varName', 'highs']

    def __init__(self, i, highs, name=None):
        self.index = i
        self.highs = highs
        self.name = f"__v{i}" if name == None else name

    def __repr__(self):
        return f"{self.name}"

    @property
    def name(self):
        if self.index < self.highs.numVars and self.highs.numVars > 0:
            return self.highs.getLp().col_names_[self.index]
        else:
            return self._varName

    @name.setter
    def name(self, value):
        if value == None or len(value) == 0:
            raise Exception('Name cannot be empty')

        self._varName = value
        if self.index < self.highs.numVars and self.highs.numVars > 0:
            self.highs.passColName(self.index, self._varName)

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
    __slots__ = ['index', '_constrName', 'highs']
    
    def __init__(self, i, highs, name):
        self.index = i
        self.highs = highs
        self.name = f"__c{i}" if name == None else name

    def __repr__(self):
        return f"{self.name}"

    @property
    def name(self):
        return self._constrName

    @name.setter
    def name(self, value):
        if value == None or len(value) == 0:
            raise Exception('Name cannot be empty')

        self._constrName = value
        self.highs.passRowName(self.index, self._constrName)
   

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

# used to batch add new variables
class highs_batch(object):
    """Batch constraint builder for HiGHS"""
    __slots__ = ['obj', 'lb', 'ub', 'type', 'name', 'highs', 'idx']

    def __init__(self, highs):
        self.highs = highs

        self.obj = []
        self.lb = []
        self.ub = []
        self.type = []
        self.idx = []
        self.name = []

    def add(self, obj, lb, ub, type, name, solver):
        self.obj.append(obj)
        self.lb.append(lb)
        self.ub.append(ub)
        self.type.append(type)
        self.name.append(name)

        newIndex = self.highs.numVars + len(self.obj)-1
        self.idx.append(newIndex)
        return highs_var(newIndex, solver, name)
