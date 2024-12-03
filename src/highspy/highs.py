from __future__ import annotations
import sys
import numpy as np
from numbers import Integral
from itertools import product
from threading import Thread, local, RLock, Lock
from typing import Optional, Any, overload, Callable, Sequence, Mapping, Iterable, SupportsIndex, cast, Union

from ._core import (
    ObjSense,
    HighsVarType,
    HighsStatus,
    cb,  # type: ignore
    _Highs,  # type: ignore
    readonly_ptr_wrapper_double,
    kHighsInf,
)

# backwards typing support information for HighspyArray
np_version = tuple(map(int, np.__version__.split('.')))
if sys.version_info >= (3, 9) and np_version >= (1,22,0):
    ndarray_object_type = np.ndarray[Any, np.dtype[np.object_]]
else:
    ndarray_object_type = np.ndarray

class Highs(_Highs):
    """
    HiGHS solver interface
    """

    __handle_keyboard_interrupt: bool = False
    __handle_user_interrupt: bool = False
    __solver_should_stop: bool = False
    __solver_stopped: RLock = RLock()
    __solver_started: Lock = Lock()
    __solver_status: Optional[HighsStatus] = None

    def __init__(self):
        super().__init__()
        self.callbacks = [HighsCallback(cb.HighsCallbackType(_), self) for _ in range(int(cb.HighsCallbackType.kCallbackMax) + 1)]
        self.enableCallbacks()

    # Silence logging
    def silent(self, turn_off_output: bool = True):
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
        if not self.HandleKeyboardInterrupt:
            return super().run()
        else:
            return self.joinSolve(self.startSolve())

    def startSolve(self):
        """
        Starts the solver in a separate thread.  Useful for handling KeyboardInterrupts.
        Do not attempt to modify the model while the solver is running.

        Returns:
            A Thread object representing the solver thread.
        """
        if not self.is_solver_running():
            self.__solver_started.acquire()
            self.__solver_should_stop = False
            self.__solver_status = None

            t = Thread(target=Highs.__solve, args=(self,), daemon=True)
            t.start()

            # wait for solver thread to start to avoid synchronization issues
            try:
                self.__solver_started.acquire(True)
            finally:
                self.__solver_started.release()
                return t
        else:
            raise Exception("Solver is already running.")

    def is_solver_running(self):
        is_running = True
        try:
            # try to acquire lock, if we can't, solver is already running
            is_running = not self.__solver_stopped.acquire(False)
            return is_running
        finally:
            if not is_running:
                self.__solver_stopped.release()

    # internal solve method for use with threads
    # will set the status of the solver when finished and release the shared lock
    def __solve(self):
        try:
            self.__solver_stopped.acquire(True)
            self.__solver_started.release()  # allow main thread to continue
            self.__solver_status = super().run()

            # avoid potential deadlock in Windows
            # can remove once HiGHS is updated to handle this internally
            _Highs.resetGlobalScheduler(False)
        finally:
            self.__solver_stopped.release()

    def joinSolve(self, solver_thread: Optional[Thread] = None, interrupt_limit: int = 5):
        """
        Waits for the solver to finish. If solver_thread is provided, it will handle KeyboardInterrupts.

        Args:
            solver_thread: A Thread object representing the solver thread (optional).
            interrupt_limit: The number of times to allow KeyboardInterrupt before forcing termination (optional).

        Returns:
            A HighsStatus object containing the solve status.
        """
        if solver_thread is not None and interrupt_limit <= 0:
            result = (False, None)

            try:
                while not result[0]:
                    result = self.wait(0.1)
                return result[1]

            except KeyboardInterrupt:
                print("KeyboardInterrupt: Waiting for HiGHS to finish...")
                self.cancelSolve()

        elif interrupt_limit > 0:
            result = (False, None)

            for count in range(interrupt_limit):
                try:
                    while not result[0]:
                        result = self.wait(0.1)
                    return result[1]

                except KeyboardInterrupt:
                    print(f"Ctrl-C pressed {count+1} times: Waiting for HiGHS to finish. ({interrupt_limit} times to force termination)")
                    self.cancelSolve()

            # if we reach this point, we should force termination
            print("Forcing termination...")
            exit(1)

        try:
            # wait for shared lock, i.e., solver to finish
            self.__solver_stopped.acquire(True)
        except KeyboardInterrupt:
            pass  # ignore additional KeyboardInterrupt here
        finally:
            self.__solver_stopped.release()

        return self.__solver_status

    def wait(self, timeout: float = -1.0):
        result = False, None

        try:
            result = (
                self.__solver_stopped.acquire(True, timeout=timeout),
                self.__solver_status,
            )
            return result
        finally:
            if result[0]:
                self.__solver_status = None  # reset status
                self.__solver_stopped.release()

    def optimize(self):
        """
        Alias for the solve method.
        """
        return self.solve()

    # reset the objective and sense, then solve
    def minimize(self, obj: Optional[Union[highs_var, highs_linear_expression]] = None):
        """
        Solves a minimization of the objective and optionally updates the costs.

        Args:
            obj: An optional highs_linear_expression representing the new objective function.

        Raises:
            Exception: If obj is an inequality or not a highs_linear_expression.

        Returns:
            A HighsStatus object containing the solve status after minimization.
        """
        if obj is not None:
            # if we have a single variable, wrap it in a linear expression
            expr = highs_linear_expression(obj) if isinstance(obj, highs_var) else obj

            if expr.bounds is not None:
                raise Exception("Objective cannot be an inequality")

            # reset objective
            super().changeColsCost(
                self.numVariables,
                np.arange(self.numVariables, dtype=np.int32),
                np.full(self.numVariables, 0, dtype=np.float64),
            )

            # if we have duplicate variables, add the vals
            idxs, vals = expr.unique_elements()
            super().changeColsCost(len(idxs), idxs, vals)
            super().changeObjectiveOffset(expr.constant or 0.0)

        super().changeObjectiveSense(ObjSense.kMinimize)
        return self.solve()

    # reset the objective and sense, then solve
    def maximize(self, obj: Optional[Union[highs_var, highs_linear_expression]] = None):
        """
        Solves a maximization of the objective and optionally updates the costs.

        Args:
            obj: An optional highs_linear_expression representing the new objective function.

        Raises:
            Exception: If obj is an inequality or not a highs_linear_expression.

        Returns:
            A HighsStatus object containing the solve status after maximization.
        """
        if obj is not None:
            # if we have a single variable, wrap it in a linear expression
            expr = highs_linear_expression(obj) if isinstance(obj, highs_var) else obj

            if expr.bounds is not None:
                raise Exception("Objective cannot be an inequality")

            # reset objective
            super().changeColsCost(
                self.numVariables,
                np.arange(self.numVariables, dtype=np.int32),
                np.full(self.numVariables, 0, dtype=np.float64),
            )

            # if we have duplicate variables, add the vals
            idxs, vals = expr.unique_elements()
            super().changeColsCost(len(idxs), idxs, vals)
            super().changeObjectiveOffset(expr.constant or 0.0)

        super().changeObjectiveSense(ObjSense.kMaximize)
        return self.solve()

    @staticmethod
    def internal_get_value(
        array_values: Union[Sequence[float], np.ndarray[Any, np.dtype[np.float64]], readonly_ptr_wrapper_double],
        index_collection: Union[
            Integral, highs_var, highs_cons, highs_linear_expression, Mapping[Any, Any], Sequence[Any], np.ndarray[Any, np.dtype[Any]]
        ],
    ) -> Union[float, bool, Mapping[Any, Any], np.ndarray[Any, np.dtype[np.float64]]]:
        """
        Internal method to get the value of an index from an array of values. Could be value or dual, variable or constraint.
        """
        if isinstance(index_collection, (Integral, highs_var, highs_cons)):
            return array_values[int(index_collection)]

        elif isinstance(index_collection, highs_linear_expression):
            return index_collection.evaluate(array_values)

        elif isinstance(index_collection, Mapping):
            return {k: Highs.internal_get_value(array_values, v) for k, v in index_collection.items()}

        else:
            return np.asarray([Highs.internal_get_value(array_values, v) for v in index_collection])

    @overload
    def val(self, var: Union[Integral, highs_var, highs_cons]) -> float: ...

    @overload
    def val(self, var: highs_linear_expression) -> Union[float, bool]: ...

    @overload
    def val(self, var: Mapping[Any, Any]) -> Mapping[Any, Any]: ...

    @overload
    def val(self, var: Union[Sequence[Any], np.ndarray[Any, np.dtype[Any]]]) -> np.ndarray[Any, np.dtype[np.float64]]: ...

    def val(
        self,
        var: Union[Integral, highs_var, highs_cons, highs_linear_expression, Mapping[Any, Any], Sequence[Any], np.ndarray[Any, np.dtype[Any]]],
    ):
        """
        Gets the value of a variable/index or expression in the solution.

        Args:
            var: A highs_var/index or highs_linear_expression object representing the variable.

        Returns:
            The value of the variable in the solution.
        """
        return Highs.internal_get_value(super().getSolution().col_value, var)

    @overload
    def vals(self, idxs: Union[Integral, highs_var, highs_cons]) -> float: ...

    @overload
    def vals(self, idxs: highs_linear_expression) -> Union[float, bool]: ...

    @overload
    def vals(self, idxs: Mapping[Any, Any]) -> Mapping[Any, Any]: ...

    @overload
    def vals(self, idxs: Union[Sequence[Any], np.ndarray[Any, np.dtype[Any]]]) -> np.ndarray[Any, np.dtype[np.float64]]: ...

    def vals(
        self,
        idxs: Union[Integral, highs_var, highs_cons, highs_linear_expression, Mapping[Any, Any], Sequence[Any], np.ndarray[Any, np.dtype[Any]]],
    ):
        """
        Gets the values of multiple variables in the solution.

        Args:
            idxs: A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.

        Returns:
            If idxs is a Mapping, returns a dict where keys are the same keys from the input idxs and values are the solution values of the corresponding variables. If idxs is an iterable, returns a list of solution values for the variables.
        """
        return Highs.internal_get_value(super().getSolution().col_value, idxs)

    def variableName(self, var: Union[Integral, highs_var]):
        """
        Retrieves the name of a specific variable.

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
            raise Exception("Variable name not found")
        return name

    @overload
    def variableNames(self, idxs: Mapping[Any, Union[highs_var, Integral]]) -> dict[Any, str]: ...

    @overload
    def variableNames(self, idxs: Iterable[Union[highs_var, Integral]]) -> list[str]: ...

    def variableNames(self, idxs: Iterable[Union[highs_var, Integral]]):
        """
        Retrieves the names of multiple variables.

        Args:
            idxs: An iterable of highs_var objects or a mapping where keys are identifiers and values are highs_var objects.

        Raises:
            Exception: If any variable name cannot be found.

        Returns:
            If idxs is a mapping, returns a dict where keys are the same keys from the input idxs and values are the names of the corresponding variables.
            If idxs is an iterable, returns a list of names for the specified variables.
        """
        if isinstance(idxs, Mapping):
            convert: Mapping[Any, Union[highs_var, Integral]] = idxs
            return {key: self.variableName(v) for key, v in convert.items()}
        else:
            return [self.variableName(v) for v in idxs]

    def allVariableNames(self) -> list[str]:
        """
        Retrieves the names of all variables in the model.

        Returns:
            A list of strings representing the names of all variables.
        """
        return super().getLp().col_names_

    def variableValue(
        self,
        var: Union[Integral, highs_var, highs_cons, highs_linear_expression, Mapping[Any, Any], np.ndarray[Any, np.dtype[Any]]],
    ):
        """
        Retrieves the value of a specific variable in the solution.

        Args:
            var: A highs_var object representing the variable.

        Returns:
            The value of the specified variable in the solution.
        """
        return self.val(var)

    def variableValues(
        self,
        idxs: Union[Integral, highs_var, highs_cons, highs_linear_expression, Mapping[Any, Any], np.ndarray[Any, np.dtype[Any]]],
    ):
        """
        Retrieves the values of multiple variables in the solution.

        Args:
            idxs: A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.

        Returns:
            If idxs is a Mapping, returns a dict where keys are the same keys from the input idxs and values are the solution values of the corresponding variables. If idxs is an iterable, returns a list of solution values for the variables.
        """
        return self.vals(idxs)

    def allVariableValues(self):
        """
        Retrieves the values of all variables in the solution.

        Returns:
            A list of values for all variables in the solution.
        """
        return super().getSolution().col_value

    def variableDual(
        self,
        var: Union[Integral, highs_var, highs_cons, highs_linear_expression, Mapping[Any, Any], np.ndarray[Any, np.dtype[Any]]],
    ):
        """
        Retrieves the dual value of a specific variable/index or expression in the solution.

        Args:
            var: A highs_var object representing the variable.

        Returns:
            The dual value of the specified variable in the solution.
        """
        return Highs.internal_get_value(super().getSolution().col_dual, var)

    def variableDuals(
        self,
        idxs: Union[Integral, highs_var, highs_cons, highs_linear_expression, Mapping[Any, Any], np.ndarray[Any, np.dtype[Any]]],
    ):
        """
        Retrieves the dual values of multiple variables in the solution.

        Args:
            idxs: A collection of highs_var objects representing the variables. Can be a Mapping (e.g., dict) where keys are variable names and values are highs_var objects, or an iterable of highs_var objects.

        Returns:
            If idxs is a Mapping, returns a dict where keys are the same keys from the input idxs and values are the dual values of the corresponding variables. If idxs is an iterable, returns a list of dual values for the variables.
        """
        return Highs.internal_get_value(super().getSolution().col_dual, idxs)

    def allVariableDuals(self):
        """
        Retrieves the dual values of all variables in the solution.

        Returns:
            A list of dual values for all variables in the solution.
        """
        return super().getSolution().col_dual

    def constrValue(
        self,
        con: Union[Integral, highs_var, highs_cons, highs_linear_expression, Mapping[Any, Any], np.ndarray[Any, np.dtype[Any]]],
    ):
        """
        Retrieves the value of a specific constraint in the solution.

        Args:
            con: A highs_con object representing the constraint.

        Returns:
            The value of the specified constraint in the solution.
        """
        return Highs.internal_get_value(super().getSolution().row_value, con)

    def constrValues(
        self,
        cons: Union[Integral, highs_var, highs_cons, highs_linear_expression, Mapping[Any, Any], np.ndarray[Any, np.dtype[Any]]],
    ):
        """
        Retrieves the values of multiple constraints in the solution.

        Args:
            cons: A collection of highs_con objects representing the constraints. Can be a Mapping (e.g., dict) where keys are constraint names and values are highs_con objects, or an iterable of highs_con objects.

        Returns:
            If cons is a Mapping, returns a dict where keys are the same keys from the input cons and values are the solution values of the corresponding constraints. If cons is an iterable, returns a list of solution values for the constraints.
        """
        return Highs.internal_get_value(super().getSolution().row_value, cons)

    def allConstrValues(self):
        """
        Retrieves the values of all constraints in the solution.

        Returns:
            A list of values for all constraints in the solution.
        """
        return super().getSolution().row_value

    def constrDual(
        self,
        con: Union[Integral, highs_var, highs_cons, highs_linear_expression, Mapping[Any, Any], np.ndarray[Any, np.dtype[Any]]],
    ):
        """
        Retrieves the dual value of a specific constraint in the solution.

        Args:
            con: A highs_con object representing the constraint.

        Returns:
            The dual value of the specified constraint in the solution.
        """
        return Highs.internal_get_value(super().getSolution().row_dual, con)

    def constrDuals(
        self,
        cons: Union[Integral, highs_var, highs_cons, highs_linear_expression, Mapping[Any, Any], np.ndarray[Any, np.dtype[Any]]],
    ):
        """
        Retrieves the dual values of multiple constraints in the solution.

        Args:
            cons: A collection of highs_con objects representing the constraints. Can be a Mapping (e.g., dict) where keys are constraint names and values are highs_con objects, or an iterable of highs_con objects.

        Returns:
            If cons is a Mapping, returns a dict where keys are the same keys from the input cons and values are the dual values of the corresponding constraints. If cons is an iterable, returns a list of dual values for the constraints.
        """
        return Highs.internal_get_value(super().getSolution().row_dual, cons)

    def allConstrDuals(self):
        """
        Retrieves the dual values of all constraints in the solution.

        Returns:
            A list of dual values for all constraints in the solution.
        """
        return super().getSolution().row_dual

    def addVariable(
        self,
        lb: float = 0,
        ub: float = kHighsInf,
        obj: float = 0.0,
        type: HighsVarType = HighsVarType.kContinuous,
        name: Optional[str] = None,
    ):
        """
        Adds a variable to the model.

        Args:
            lb: Lower bound of the variable (default is 0).
            ub: Upper bound of the variable (default is infinity).
            obj: Objective coefficient of the variable (default is 0).
            type: Type of the variable (continuous, integer; default is continuous).
            name: Optional name for the variable.

        Returns:
            A highs_var object representing the added variable.
        """
        status = super().addCol(obj, lb, ub, 0, np.empty(0, dtype=np.int32), np.empty(0, np.float64))

        if status != HighsStatus.kOk:
            raise Exception("Failed to add variable to the model.")

        var = highs_var(self.numVariables - 1, self)

        if type != HighsVarType.kContinuous:
            super().changeColIntegrality(var.index, type)

        if name is not None:
            super().passColName(var.index, name)

        return var

    @overload
    def addVariables(
        self,
        *nvars: int,
        **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]],
    ) -> HighspyArray: ...

    @overload
    def addVariables(
        self,
        *nvars: Mapping[Any, Any],
        **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]],
    ) -> dict[Any, highs_var]: ...

    @overload
    def addVariables(
        self,
        *nvars: Sequence[Any],
        **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]],
    ) -> Union[dict[Any, highs_var], HighspyArray]: ...

    @overload
    def addVariables(
        self,
        *nvars: Union[int, Mapping[Any, Any], Sequence[Any]],
        **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]],
    ) -> Optional[Union[dict[Any, highs_var], HighspyArray]]: ...

    def addVariables(
        self,
        *nvars: Union[int, Mapping[Any, Any], Sequence[Any]],
        **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]],
    ) -> Optional[Union[dict[Any, highs_var], HighspyArray]]:
        """
        Adds multiple variables to the model.

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

        # if all nvars are scalars, we can assume they are the dimensions
        shape = [n for n in nvars if isinstance(n, int)]

        expanded = [range(n) if isinstance(n, int) else n for n in nvars]
        indices: list[Any] = list(expanded[0] if len(expanded) == 1 else product(*expanded))  # unpack tuple if needed
        N = len(indices)

        if len(shape) != len(nvars):
            shape = [N]

        # parameter can be scalar, array, or mapping lookup (i.e., dictionary, custom class, etc.)
        # scalar: repeat for all N, array: use as is, lookup: convert to array using indices
        def ensure_real(x: Union[Any, Mapping[Any, Any], Sequence[Any]]):
            if isinstance(x, (float, int)):
                return np.full(N, x, dtype=np.float64)

            elif isinstance(x, Mapping):
                mt: Mapping[Any, Any] = x

                if all(isinstance(v, (float, int)) for v in mt.values()):
                    m: Mapping[Any, Union[float, int]] = x
                    return np.fromiter((m[i] for i in indices), np.float64)

            elif isinstance(x, Sequence) and len(x) == N and all(isinstance(v, (float, int)) for v in x):
                return np.asarray(x, dtype=np.float64)

            raise Exception("Invalid parameter.")

        def ensure_HighsVarType(x: Union[Any, Mapping[Any, Any], Sequence[Any]]):
            if x == HighsVarType.kContinuous:
                return None

            elif isinstance(x, HighsVarType):
                return np.full(N, x, dtype=np.uint8)

            elif isinstance(x, Mapping):
                mt: Mapping[Any, Any] = x

                if all(isinstance(v, HighsVarType) for v in mt.values()):
                    m: Mapping[Any, HighsVarType] = x
                    return np.fromiter((m[i] for i in indices), np.uint8)

            elif isinstance(x, Sequence) and len(x) == N and all(isinstance(v, HighsVarType) for v in x):
                return np.asarray(x, dtype=np.uint8)

            raise Exception("Invalid parameter.")

        def ensure_optional_str(x: Optional[Any]):
            if x is None:
                return None
            elif isinstance(x, Sequence):
                mt: Sequence[Any] = x

                if len(mt) == N and all(isinstance(v, str) for v in mt):
                    m: Sequence[str] = mt
                    return m

            raise Exception("Invalid parameter.")

        def ensure_str_or_none(x: Any) -> Optional[str]:
            if x is None or isinstance(x, str):
                return x
            else:
                raise Exception("Invalid parameter.")

        def ensure_bool(x: Any) -> bool:
            if isinstance(x, bool):
                return x
            raise Exception("Invalid parameter.")

        lb = ensure_real(kwargs.get("lb", 0.0))
        ub = ensure_real(kwargs.get("ub", kHighsInf))
        obj = ensure_real(kwargs.get("obj", 0))
        vartype = ensure_HighsVarType(kwargs.get("type", HighsVarType.kContinuous))
        name_prefix = ensure_str_or_none(kwargs.get("name_prefix", None))
        name = ensure_optional_str(kwargs.get("name", None))
        out_array = ensure_bool(kwargs.get("out_array", all(isinstance(n, int) for n in nvars)))

        start_idx = self.numVariables
        idx = np.arange(start_idx, start_idx + N, dtype=np.int32)
        status = super().addCols(
            N,
            obj,
            lb,
            ub,
            0,
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.float64),
        )

        if status != HighsStatus.kOk:
            raise Exception("Failed to add columns to the model.")

        # only set integrality if we have non-continuous variables
        if vartype is not None:
            super().changeColsIntegrality(N, idx, vartype)

        if name or name_prefix:
            names = name or [f"{name_prefix}{i}" for i in indices]

            for i, n in zip(idx, names):
                super().passColName(int(i), str(n))

        return (
            HighspyArray(np.asarray([highs_var(i, self) for i in idx]).reshape(shape), self)
            if out_array
            else {index: highs_var(i, self) for index, i in zip(indices, idx)}
        )

    @overload
    def addIntegrals(
        self,
        *nvars: int,
        **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]],
    ) -> HighspyArray: ...

    @overload
    def addIntegrals(
        self,
        *nvars: Mapping[Any, Any],
        **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]],
    ) -> dict[Any, highs_var]: ...

    @overload
    def addIntegrals(
        self,
        *nvars: Sequence[Any],
        **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]],
    ) -> Union[dict[Any, highs_var], HighspyArray]: ...

    @overload
    def addIntegrals(
        self,
        *nvars: Union[int, Mapping[Any, Any], Sequence[Any]],
        **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]],
    ) -> Optional[Union[dict[Any, highs_var], HighspyArray]]: ...

    def addIntegrals(
        self,
        *nvars: Union[int, Mapping[Any, Any], Sequence[Any]],
        **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]],
    ):
        """
        Alias for the addVariables method, for integer variables.
        """
        kwargs.setdefault("type", HighsVarType.kInteger)
        return self.addVariables(*nvars, **kwargs)

    @overload
    def addBinaries(
        self,
        *nvars: int,
        **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]],
    ) -> HighspyArray: ...

    @overload
    def addBinaries(
        self,
        *nvars: Mapping[Any, Any],
        **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]],
    ) -> dict[Any, highs_var]: ...

    @overload
    def addBinaries(
        self,
        *nvars: Sequence[Any],
        **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]],
    ) -> Union[dict[Any, highs_var], HighspyArray]: ...

    @overload
    def addBinaries(
        self,
        *nvars: Union[int, Mapping[Any, Any], Sequence[Any]],
        **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]],
    ) -> Optional[Union[dict[Any, highs_var], HighspyArray]]: ...

    def addBinaries(
        self,
        *nvars: Union[int, Mapping[Any, Any], Sequence[Any]],
        **kwargs: Union[Any, HighsVarType, Mapping[Any, Any], Sequence[Any]],
    ):
        """
        Alias for the addVariables method, for binary variables.
        """
        kwargs.setdefault("lb", 0)
        kwargs.setdefault("ub", 1)
        kwargs.setdefault("type", HighsVarType.kInteger)

        return self.addVariables(*nvars, **kwargs)

    def addIntegral(self, lb: float = 0.0, ub: float = kHighsInf, obj: float = 0.0, name: Optional[str] = None):
        """
        Alias for the addVariable method, for integer variables.
        """
        return self.addVariable(lb, ub, obj, HighsVarType.kInteger, name)

    def addBinary(self, obj: float = 0.0, name: Optional[str] = None):
        """
        Alias for the addVariable method, for binary variables.
        """
        return self.addVariable(0, 1, obj, HighsVarType.kInteger, name)

    def deleteVariable(
        self,
        var_or_index: Union[Integral, highs_var],
        *args: Union[Mapping[Any, highs_var], Iterable[highs_var], highs_var],
    ):
        """
        Deletes a variable from the model and updates the indices of subsequent variables in provided collections.

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
            if isinstance(collection, Mapping):
                mt: Mapping[Any, highs_var] = collection

                # Update indices in a dictionary of variables
                for _, var in mt.items():
                    if var.index > index:
                        var.index -= 1
                    elif var.index == index:
                        var.index = -1

            elif isinstance(collection, Iterable):
                # Update indices in an iterable of variables
                for var in collection:
                    if var.index > index:
                        var.index -= 1
                    elif var.index == index:
                        var.index = -1

            # If the collection is a single highs_var object, check and update if necessary
            elif collection.index > index:
                collection.index -= 1
            elif collection.index == index:
                collection.index = -1

    def getVariables(self):
        """
        Retrieves all variables in the model.

        Returns:
            A list of highs_var objects, each representing a variable in the model.
        """
        return [highs_var(i, self) for i in range(self.numVariables)]

    @property
    def inf(self):
        """
        Represents infinity in the context of the solver.

        Returns:
            The value used to represent infinity.
        """
        return kHighsInf

    @property
    def numVariables(self):
        """
        Gets the number of variables in the model.

        Returns:
            The number of variables.
        """
        return super().getNumCol()

    @property
    def numConstrs(self):
        """
        Gets the number of constraints in the model.

        Returns:
            The number of constraints.
        """
        return super().getNumRow()

    #
    # add constraints
    #
    def addConstr(self, expr: highs_linear_expression, name: Optional[str] = None):
        """
        Adds a constraint to the model.

        Args:
            expr: A highs_linear_expression to be added.
            name: Optional name of the constraint.

        Returns:
            A highs_con object representing the added constraint.
        """
        con = self.__addRow(expr, self.numConstrs)

        if name is not None:
            super().passRowName(con.index, name)

        return con

    @overload
    def addConstrs(
        self,
        *args: highs_linear_expression,
        **kwargs: Optional[Union[str, Sequence[str]]],
    ) -> list[highs_cons]: ...

    @overload
    def addConstrs(
        self,
        *args: Mapping[Any, highs_linear_expression],
        **kwargs: Optional[Union[str, Sequence[str]]],
    ) -> dict[Any, highs_cons]: ...

    @overload
    def addConstrs(
        self,
        *args: Iterable[highs_linear_expression],
        **kwargs: Optional[Union[str, Sequence[str]]],
    ) -> list[highs_cons]: ...

    def addConstrs(
        self,
        *args: Union[highs_linear_expression, Iterable[highs_linear_expression]],
        **kwargs: Optional[Union[str, Sequence[str]]],
    ):
        """
        Adds multiple constraints to the model.

        Args:
            *args: A sequence of highs_linear_expression to be added.

            **kwargs: Optional keyword arguments.
                name_prefix: Prefix for the constraint names.  Constructed name will be name_prefix + index.
                name: A collection of names for the constraints (list or mapping).

        Returns:
            A highs_con collection array representing the added constraints.
        """
        name_prefix = kwargs.get("name_prefix", None)
        name = kwargs.get("name", None)

        # unpack generator if needed
        generator = args[0] if len(args) == 1 and isinstance(args[0], Iterable) else args
        initial_rows = self.numConstrs

        try:
            if isinstance(generator, Mapping):
                mt: Mapping[Any, highs_linear_expression] = generator
                cons = {key: self.__addRow(expr, initial_rows + count) for count, (key, expr) in enumerate(mt.items())}
            else:
                it: Iterable[Any] = generator
                cons = [self.__addRow(expr, initial_rows + count) for count, expr in enumerate(it)]

            # TODO: Mapping support with constraint names can be improved, e.g., by allowing a name collection to be passed
            if name or name_prefix:
                names = name or [f"{name_prefix}{n}" for n in range(self.numConstrs - initial_rows)]

                for c, n in zip(range(initial_rows, self.numConstrs), names):
                    super().passRowName(int(c), str(n))

        except Exception as e:
            # rollback model if error - remove any constraints that were added
            status = super().deleteRows(
                self.numConstrs - initial_rows,
                np.arange(initial_rows, self.numConstrs, dtype=np.int32),
            )

            if status != HighsStatus.kOk:
                raise Exception("Failed to rollback model after failure.  Model might be in a undeterminate state.")
            else:
                raise e

        return cons

    def __addRow(self, expr: highs_linear_expression, idx: int):
        """
        Internal method to add a constraint to the model.
        """
        if expr.bounds is not None:
            idxs, vals = expr.unique_elements()  # if we have duplicate variables, add the vals
            super().addRow(expr.bounds[0], expr.bounds[1], len(idxs), idxs, vals)
            return highs_cons(idx, self)
        else:
            raise Exception("Constraint bounds must be set via comparison (>=,==,<=).")

    def expr(
        self,
        optional: Optional[Union[float, int, highs_var, highs_linear_expression]] = None,
    ):
        """
        Creates a new highs_linear_expression object.

        Returns:
            A highs_linear_expression object.
        """
        return highs_linear_expression(optional)

    def getExpr(self, cons: Union[Integral, highs_cons]):
        """
        Retrieves the highs_linear_expression of a constraint.

        Args:
            cons: A highs_con object or index representing the constraint.

        Returns:
            A highs_linear_expression object representing the expression of the constraint.
        """
        status, lb, ub, nnz = super().getRow(int(cons))  # type: ignore

        if status != HighsStatus.kOk:
            raise Exception("Error retrieving constraint expression.")

        status, idx, val = super().getRowEntries(int(cons))

        if status != HighsStatus.kOk:
            raise Exception("Error retrieving constraint expression entries.")

        expr = highs_linear_expression()
        expr.bounds = (lb, ub)
        expr.idxs = list(idx)
        expr.vals = list(val)
        return expr

    def chgCoeff(self, cons: Union[highs_cons, Integral], var: Union[highs_var, Integral], val: float):
        """
        Changes the coefficient of a variable in a constraint.

        Args:
            cons: A highs_con object representing the constraint.
            var: A highs_var object representing the variable.
            val: The new coefficient value for the variable in the constraint.
        """
        super().changeCoeff(int(cons), int(var), val)

    def getConstrs(self):
        """
        Retrieves all constraints in the model.

        Returns:
            A list of highs_cons objects, each representing a constraint in the model.
        """
        return [highs_cons(i, self) for i in range(self.numConstrs)]

    def removeConstr(
        self,
        cons_or_index: Union[highs_cons, Integral],
        *args: Union[Mapping[Any, highs_cons], Sequence[highs_cons], highs_cons],
    ):
        """
        Removes a constraint from the model and updates the indices of subsequent constraints in provided collections.

        Args:
            cons_or_index: A highs_cons object or an index representing the constraint to be removed.
            *args: Optional collections (lists, dicts, etc.) of highs_cons objects whose indices need to be updated after the removal.
        """
        # Determine the index of the constraint to delete
        index = int(cons_or_index)

        # Delete the variable from the model if it exists
        if index < self.numConstrs:
            status = super().deleteRows(1, np.asarray([index], dtype=np.int32))

            if status != HighsStatus.kOk:
                raise Exception("Failed to delete constraint from the model.")

            # Update the indices of constraints in the provided collections
            for collection in args:
                if isinstance(collection, Mapping):
                    mt: Mapping[Any, highs_cons] = collection

                    # Update indices in a dictionary of constraints
                    for _, con in mt.items():
                        if con.index > index:
                            con.index -= 1
                        elif con.index == index:
                            con.index = -1

                elif isinstance(collection, Sequence):
                    # Update indices in an iterable of constraints
                    for con in collection:
                        if con.index > index:
                            con.index -= 1
                        elif con.index == index:
                            con.index = -1

                # If the collection is a single highs_cons object, check and update if necessary
                elif collection.index > index:
                    collection.index -= 1

                elif collection.index == index:
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

    def setInteger(self, var_or_collection: Union[highs_var, int, Iterable[Union[highs_var, int]]]):
        """
        Sets a variable/collection to integer.

        Args:
            var_or_collection: A highs_var object/collection representing the variable to be set as integer.
        """
        if isinstance(var_or_collection, Iterable):
            idx = np.fromiter(map(int, var_or_collection), dtype=np.int32)
            super().changeColsIntegrality(len(idx), idx, np.full(len(idx), HighsVarType.kInteger, dtype=np.uint8))
        else:
            super().changeColIntegrality(int(var_or_collection), HighsVarType.kInteger)

    def setContinuous(self, var_or_collection: Union[highs_var, int, Iterable[Union[highs_var, int]]]):
        """
        Sets a variable/collection to continuous.

        Args:
            var_or_collection: A highs_var object/collection representing the variable to be set as continuous.
        """
        if isinstance(var_or_collection, Iterable):
            idx = np.fromiter(map(int, var_or_collection), dtype=np.int32)
            super().changeColsIntegrality(
                len(idx),
                idx,
                np.full(len(idx), HighsVarType.kContinuous, dtype=np.uint8),
            )
        else:
            super().changeColIntegrality(int(var_or_collection), HighsVarType.kContinuous)

    @staticmethod
    def qsum(
        items: Union[Iterable[Union[highs_var, highs_linear_expression]], np.ndarray[Any, np.dtype[np.object_]]],
        initial: Optional[Union[float, int, highs_var, highs_linear_expression]] = None,
    ):
        """
        Performs a faster sum for highs_linear_expressions.

        Args:
            items: A collection of highs_linear_expressions or highs_vars to be summed.
        """
        expr = highs_linear_expression(initial)

        if isinstance(items, np.ndarray):
            X: np.ndarray[Any, np.dtype[np.object_]] = items
            for v in X.flat:
                expr += cast(Union[highs_var, highs_linear_expression], v)
        else:
            for v in items:
                expr += v

        return expr

    ##
    ## Callback support, with optional user interrupt handling
    ##
    @staticmethod
    def __internal_callback(
        callback_type: cb.HighsCallbackType,
        message: str,
        data_out: cb.HighsCallbackDataOut,
        data_in: Optional[cb.HighsCallbackDataIn],
        user_callback_data: Any,
    ):
        user_callback_data.callbacks[int(callback_type)].fire(message, data_out, data_in)

    def enableCallbacks(self):
        """
        Enables callbacks, restarting them if they were previously enabled.
        """
        super().setCallback(Highs.__internal_callback, self)

        # restart callbacks if any exist
        for c in self.callbacks:
            if len(c.callbacks) > 0:
                self.startCallback(c.callback_type)

    def disableCallbacks(self):
        """
        Disables all callbacks.
        """
        status = super().setCallback(None, None)  # this will also stop all callbacks

        if status != HighsStatus.kOk:
            raise Exception("Failed to disable callbacks.")

    def cancelSolve(self):
        """
        If HandleUserInterrupt is enabled, this method will signal the solver to stop.
        """
        self.__solver_should_stop = True

    @property
    def HandleKeyboardInterrupt(self):
        """
        Get/Set whether the solver should handle KeyboardInterrupt (i.e., cancel solve on Ctrl+C). Also enables/disables HandleUserInterrupt.
        """
        return self.__handle_keyboard_interrupt

    @HandleKeyboardInterrupt.setter
    def HandleKeyboardInterrupt(self, value: bool):
        self.__handle_keyboard_interrupt = value
        self.HandleUserInterrupt = value

    @property
    def HandleUserInterrupt(self):
        """
        Get/Set whether the solver should handle user interrupts (i.e., cancel solve on user request)
        """
        return self.__handle_user_interrupt

    @HandleUserInterrupt.setter
    def HandleUserInterrupt(self, value: bool):
        self.__handle_user_interrupt = value

        if value:
            self.cbSimplexInterrupt += self.__user_interrupt_event
            self.cbIpmInterrupt += self.__user_interrupt_event
            self.cbMipInterrupt += self.__user_interrupt_event
        else:
            self.cbSimplexInterrupt -= self.__user_interrupt_event
            self.cbIpmInterrupt -= self.__user_interrupt_event
            self.cbMipInterrupt -= self.__user_interrupt_event

    def __user_interrupt_event(self, e: HighsCallbackEvent):
        if self.__solver_should_stop and e.data_in is not None:
            e.data_in.user_interrupt = True

    @property
    def cbLogging(self):
        return self.callbacks[int(cb.HighsCallbackType.kCallbackLogging)]

    @property
    def cbSimplexInterrupt(self):
        return self.callbacks[int(cb.HighsCallbackType.kCallbackSimplexInterrupt)]

    @property
    def cbIpmInterrupt(self):
        return self.callbacks[int(cb.HighsCallbackType.kCallbackIpmInterrupt)]

    @property
    def cbMipSolution(self):
        return self.callbacks[int(cb.HighsCallbackType.kCallbackMipSolution)]

    @property
    def cbMipImprovingSolution(self):
        return self.callbacks[int(cb.HighsCallbackType.kCallbackMipImprovingSolution)]

    @property
    def cbMipLogging(self):
        return self.callbacks[int(cb.HighsCallbackType.kCallbackMipLogging)]

    @property
    def cbMipInterrupt(self):
        return self.callbacks[int(cb.HighsCallbackType.kCallbackMipInterrupt)]

    @property
    def cbMipGetCutPool(self):
        return self.callbacks[int(cb.HighsCallbackType.kCallbackMipGetCutPool)]

    @property
    def cbMipDefineLazyConstraints(self):
        return self.callbacks[int(cb.HighsCallbackType.kCallbackMipDefineLazyConstraints)]

    # callback setters are required for +=/-= syntax
    # e.g., h.cbLogging += my_callback
    @cbLogging.setter
    def cbLogging(self, value: HighsCallback):
        if self.cbLogging is not value:
            raise Exception("Cannot set callback directly.  Use .subscribe(callback) instead.")

    @cbSimplexInterrupt.setter
    def cbSimplexInterrupt(self, value: HighsCallback):
        if self.cbSimplexInterrupt is not value:
            raise Exception("Cannot set callback directly.  Use .subscribe(callback) instead.")

    @cbIpmInterrupt.setter
    def cbIpmInterrupt(self, value: HighsCallback):
        if self.cbIpmInterrupt is not value:
            raise Exception("Cannot set callback directly.  Use .subscribe(callback) instead.")

    @cbMipSolution.setter
    def cbMipSolution(self, value: HighsCallback):
        if self.cbMipSolution is not value:
            raise Exception("Cannot set callback directly.  Use .subscribe(callback) instead.")

    @cbMipImprovingSolution.setter
    def cbMipImprovingSolution(self, value: HighsCallback):
        if self.cbMipImprovingSolution is not value:
            raise Exception("Cannot set callback directly.  Use .subscribe(callback) instead.")

    @cbMipLogging.setter
    def cbMipLogging(self, value: HighsCallback):
        if self.cbMipLogging is not value:
            raise Exception("Cannot set callback directly.  Use .subscribe(callback) instead.")

    @cbMipInterrupt.setter
    def cbMipInterrupt(self, value: HighsCallback):
        if self.cbMipInterrupt is not value:
            raise Exception("Cannot set callback directly.  Use .subscribe(callback) instead.")

    @cbMipGetCutPool.setter
    def cbMipGetCutPool(self, value: HighsCallback):
        if self.cbMipGetCutPool is not value:
            raise Exception("Cannot set callback directly.  Use .subscribe(callback) instead.")

    @cbMipDefineLazyConstraints.setter
    def cbMipDefineLazyConstraints(self, value: HighsCallback):
        if self.cbMipDefineLazyConstraints is not value:
            raise Exception("Cannot set callback directly.  Use .subscribe(callback) instead.")


##
## Callback support
##
class HighsCallbackEvent(object):
    __slots__ = ["message", "data_out", "data_in", "user_data"]

    def __init__(
        self,
        message: str,
        data_out: cb.HighsCallbackDataOut,
        data_in: Optional[cb.HighsCallbackDataIn],
        user_data: Optional[Any],
    ):
        self.message = message
        self.data_out = data_out
        self.data_in = data_in
        self.user_data = user_data

    def val(
        self,
        var_expr: Union[Integral, highs_var, highs_cons, highs_linear_expression, Mapping[Any, Any], np.ndarray[Any, np.dtype[Any]]],
    ):
        """
        Gets the value(s) of a variable/index or expression in the callback solution.
        """
        return Highs.internal_get_value(self.data_out.mip_solution, var_expr)


class HighsCallback(object):
    __slots__ = ["callbacks", "user_callback_data", "highs", "callback_type"]

    def __init__(self, callback_type: cb.HighsCallbackType, highs: Highs):
        self.callbacks: list[Callable[[HighsCallbackEvent], None]] = []
        self.user_callback_data: list[Any] = []
        self.callback_type = callback_type
        self.highs = highs

    def subscribe(
        self,
        callback: Callable[[HighsCallbackEvent], None],
        user_data: Optional[Any] = None,
    ):
        """
        Subscribes a callback to the event.

        Args:
            callback: The callback function to be executed.
            user_data: Optional user data to be passed to the callback.
        """
        if len(self.callbacks) == 0:
            status = self.highs.startCallback(self.callback_type)

            if status != HighsStatus.kOk:
                raise Exception("Failed to start callback.")

        self.callbacks.append(callback)
        self.user_callback_data.append(user_data)
        return self

    def unsubscribe(self, callback: Callable[[HighsCallbackEvent], None]):
        """
        Unsubscribes a callback from the event.

        Args:
            callback: The callback function to be removed.
        """
        try:
            idx = self.callbacks.index(callback)
            del self.callbacks[idx]
            del self.user_callback_data[idx]

            if len(self.callbacks) == 0:
                self.highs.stopCallback(self.callback_type)

        except ValueError:
            pass

        return self

    def unsubscribe_by_data(self, user_data: Optional[Any]):
        """
        Unsubscribes a callback by user data.

        Args:
            user_data: The user data corresponding to the callback(s) to be removed.
        """
        idx = reversed([i for i, ud in enumerate(self.user_callback_data) if ud == user_data])

        for i in idx:
            del self.callbacks[i]
            del self.user_callback_data[i]

        if len(self.callbacks) == 0:
            self.highs.stopCallback(self.callback_type)

        return self

    def __iadd__(self, callback: Callable[[HighsCallbackEvent], None]):
        return self.subscribe(callback)

    def __isub__(self, callback: Callable[[HighsCallbackEvent], None]):
        return self.unsubscribe(callback)

    def clear(self):
        """
        Unsubscribes all callbacks from the event.
        """
        self.callbacks = []
        self.user_callback_data = []
        self.highs.stopCallback(self.callback_type)

    def fire(
        self,
        message: str,
        data_out: cb.HighsCallbackDataOut,
        data_in: cb.HighsCallbackDataIn,
    ):
        """
        Fires the event, executing all subscribed callbacks.
        """
        e = HighsCallbackEvent(message, data_out, data_in, None)

        for fn, user_data in zip(self.callbacks, self.user_callback_data):
            e.user_data = user_data
            fn(e)


class HighspyArray(ndarray_object_type):
    """
    A numpy array wrapper for highs_var/highs_linear_expression objects.

    This provides additional type information for static analysis, and also allows faster sum operations.
    """

    def __new__(cls, input_array: np.ndarray[Any, np.dtype[np.object_]], highs: Optional[Highs]):
        obj = np.asarray(input_array).view(cls)
        obj.highs = highs
        return obj

    def __array_finalize__(self, obj: Optional[Any]):
        self.highs = getattr(obj, "highs", None)

    @overload
    def __getitem__(self, key: Union[SupportsIndex, tuple[SupportsIndex, ...]]) -> highs_linear_expression: ...  # type: ignore

    @overload
    def __getitem__(
        self,
        key: Union[
            np.ndarray[Any, np.dtype[np.integer[Any]]],
            np.ndarray[Any, np.dtype[np.bool_]],
            tuple[Union[np.ndarray[Any, np.dtype[np.integer[Any]]], np.ndarray[Any, np.dtype[np.bool_]]], ...],
        ],
    ) -> HighspyArray: ...

    @overload
    def __getitem__(self, key: Union[None, slice, SupportsIndex, tuple[Union[None, slice, SupportsIndex], ...]]) -> HighspyArray: ...

    @overload
    def __getitem__(self, key: Any) -> HighspyArray: ...  # type: ignore

    def __getitem__(self, key: Any) -> Union[HighspyArray, highs_linear_expression]:  # type: ignore
        return super(HighspyArray, self).__getitem__(key)  # type: ignore

    def __ge__(self, other: Any) -> HighspyArray:  # type: ignore
        return cast(HighspyArray, np.greater_equal(self, other, dtype=np.object_))

    def __le__(self, other: Any) -> HighspyArray:  # type: ignore
        return cast(HighspyArray, np.less_equal(self, other, dtype=np.object_))

    def __eq__(self, other: Any) -> HighspyArray:
        return cast(HighspyArray, np.equal(self, other, dtype=np.object_))

    @overload
    def sum(self, axis: None = None, dtype: Optional[Any] = None, out: None = ...) -> highs_linear_expression: ...

    @overload
    def sum(self, axis: Any, dtype: Optional[Any] = None, out: HighspyArray = ...) -> HighspyArray: ...

    def sum(  # type: ignore
        self,
        axis: Optional[int] = None,
        dtype: Optional[Any] = None,
        out: Optional[np.ndarray[Any, np.dtype[np.object_]]] = None,
        **unused_kwargs: Any,
    ) -> Union[HighspyArray, highs_linear_expression]:
        if self.highs is not None:
            if axis is not None:
                return HighspyArray(np.apply_along_axis(self.highs.qsum, axis, self, initial=unused_kwargs.get("initial", None)), self.highs)
            else:
                return self.highs.qsum(self, unused_kwargs.get("initial", None))
        else:
            raise Exception("Cannot sum without a Highs object.")


# highs variable
class highs_var(object):
    """
    Variable index wrapper for HiGHS
    """

    __slots__ = ["index", "highs"]

    def __init__(self, i: int, highs: Highs):
        self.index = i
        self.highs = highs

    def __repr__(self):
        return f"highs_var({self.index})"

    @property
    def name(self):
        return self.highs.variableName(self)

    @name.setter
    def name(self, value: str):
        self.highs.passColName(self.index, value)

    def __int__(self):
        return int(self.index)

    def __hash__(self):
        return int(self.index)

    def __neg__(self):
        expr = highs_linear_expression()
        expr.idxs = [self.index]
        expr.vals = [-1.0]
        return expr

    def __add__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        expr = highs_linear_expression(self)
        expr += other
        return expr

    def __radd__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        expr = highs_linear_expression(self)
        expr += other
        return expr

    def __mul__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        expr = highs_linear_expression(self)
        expr *= other
        return expr

    def __rmul__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        expr = highs_linear_expression(self)
        expr *= other
        return expr

    def __rsub__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        expr = highs_linear_expression(other)
        expr -= self
        return expr

    def __sub__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        expr = highs_linear_expression(self)
        expr -= other
        return expr

    # self <= other
    def __le__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        if isinstance(other, highs_linear_expression):
            return other.__ge__(self)
        else:
            return highs_linear_expression(self).__le__(other)

    # self == other
    def __eq__(self, other: Any) -> highs_linear_expression:  # type: ignore
        if isinstance(other, highs_linear_expression):
            return other.__eq__(self)
        else:
            return highs_linear_expression(self).__eq__(other)

    # self != other
    def __ne__(self, other: Optional[Any]):
        if other is None:
            return True
        else:
            raise Exception("Invalid comparison.")

    # self >= other
    def __ge__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        if isinstance(other, highs_linear_expression):
            return other.__le__(self)
        else:
            return highs_linear_expression(self).__ge__(other)


# highs constraint
class highs_cons(object):
    """
    Constraint index wrapper for HiGHS
    """

    __slots__ = ["index", "highs"]

    def __init__(self, i: int, highs: Highs):
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
    def name(self, value: str):
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
    Linear constraint builder for HiGHS
    """

    __slots__ = ["idxs", "vals", "constant", "bounds"]

    @overload
    def __init__(self, other: None = None) -> None: ...

    @overload
    def __init__(self, other: Union[float, int]) -> None: ...

    @overload
    def __init__(self, other: highs_var) -> None: ...

    @overload
    def __init__(self, other: highs_linear_expression) -> None: ...

    def __init__(
        self,
        other: Optional[Union[float, int, highs_var, highs_linear_expression]] = None,
    ):
        self.constant: Optional[float] = None  # constant is only valid when bounds are None
        self.bounds: Optional[tuple[float, float]] = None  # bounds are only valid when constant is None

        if other is None:
            self.idxs: list[int] = []
            self.vals: list[float] = []

        elif isinstance(other, highs_linear_expression):
            self.idxs = list(other.idxs)
            self.vals = list(other.vals)
            self.constant = other.constant
            self.bounds = other.bounds if other.bounds is not None else None

        elif isinstance(other, highs_var):
            self.idxs = [other.index]
            self.vals = [1.0]

        else:
            self.idxs = []
            self.vals = []
            self.constant = float(other)

    def simplify(self):
        """
        Simplifies the linear expression by combining duplicate variables.
        """
        copy = highs_linear_expression()
        copy.idxs, copy.vals = (v.tolist() for v in self.unique_elements())
        copy.bounds = self.bounds if self.bounds is not None else None
        copy.constant = self.constant
        return copy

    def copy(self):
        """
        Creates a copy of the linear expression.
        """
        return highs_linear_expression(self)

    def evaluate(
        self,
        values: Union[Sequence[float], np.ndarray[Any, np.dtype[np.float64]], readonly_ptr_wrapper_double],
    ) -> Union[float, bool]:
        """
        Evaluates the linear expression given a solution array (values).
        """
        result = sum(v * values[c] for c, v in zip(self.idxs, self.vals)) + (self.constant or 0.0)
        return result if self.bounds is None else (self.bounds[0] <= result <= self.bounds[1])

    def __repr__(self):
        # display duplicate variables
        v = str.join("  ", [f"{c}_v{x}" for x, c in zip(self.idxs, self.vals)])

        if self.bounds is None:
            return f"{v}" + (f"  {self.constant}" if self.constant is not None else "")
        elif self.bounds[0] == self.bounds[1]:
            return f"{v} == {self.bounds[0] - (self.constant or 0.0)}"
        else:
            return f"{self.bounds[0]} <= {v} <= {self.bounds[1]}"

    def __str__(self):
        # display unique variables (values are totaled)
        idxs, vals = self.unique_elements()
        v = str.join("  ", [f"{c}_v{x}" for x, c in zip(idxs, vals)])

        if self.bounds is None:
            return f"{v}" + (f"  {self.constant}" if self.constant is not None else "")
        elif self.bounds[0] == self.bounds[1]:
            return f"{v} == {self.bounds[0] - (self.constant or 0.0)}"
        else:
            return f"{self.bounds[0]} <= {v} <= {self.bounds[1]}"

    # self != other
    def __ne__(self, other: Optional[Any]):
        if other is None:
            return True
        else:
            raise Exception("Invalid comparison.")

    # self == other
    def __eq__(self, other: Any) -> highs_linear_expression:  # type: ignore
        if self.bounds is not None:
            raise Exception("Bounds have already been set.")

        # self == c
        elif isinstance(other, (float, int)):
            copy = highs_linear_expression(self)
            copy.bounds = (
                float(other) - (self.constant or 0.0),
                float(other) - (self.constant or 0.0),
            )
            copy.constant = None
            return copy

        # self == other
        elif isinstance(other, highs_linear_expression):
            if other.bounds is not None:
                raise Exception("Bounds have already been set.")

            copy = highs_linear_expression()

            # prefer most idxs, constant, 'left'
            if len(other.idxs) > len(self.idxs) or (len(other.idxs) == len(self.idxs) and other.constant is None and self.constant is not None):
                copy.idxs = other.idxs + self.idxs
                copy.vals = other.vals + [-v for v in self.vals]
                copy.bounds = (
                    (self.constant or 0.0) - (other.constant or 0.0),
                    (self.constant or 0.0) - (other.constant or 0.0),
                )
            else:
                copy.idxs = self.idxs + other.idxs
                copy.vals = self.vals + [-v for v in other.vals]
                copy.bounds = (
                    (other.constant or 0.0) - (self.constant or 0.0),
                    (other.constant or 0.0) - (self.constant or 0.0),
                )

            return copy

        # self == x
        elif isinstance(other, highs_var):
            copy = highs_linear_expression()

            if len(self.idxs) == 0 or len(self.idxs) == 1 and self.constant is not None:
                copy.idxs = [other.index] + self.idxs
                copy.vals = [1.0] + [-v for v in self.vals]
                copy.bounds = ((self.constant or 0.0), (self.constant or 0.0))
            else:
                copy.idxs = self.idxs + [other.index]
                copy.vals = self.vals + [-1.0]
                copy.bounds = (-(self.constant or 0.0), -(self.constant or 0.0))

            return copy

        # support expr == [lb, ub] --> lb <= expr <= ub
        elif hasattr(other, "__getitem__") and hasattr(other, "__len__") and len(other) == 2:
            if not (isinstance(other[0], (float, int)) and isinstance(other[1], (float, int))):
                raise Exception("Provided bounds were not valid numbers.")

            copy = highs_linear_expression(self)
            copy.bounds = (
                float(other[0]) - (copy.constant or 0.0),
                float(other[1]) - (copy.constant or 0.0),
            )
            copy.constant = None
            return copy

        else:
            raise Exception("Unknown comparison.")

    # self <= other
    def __le__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        if self.bounds is not None:
            raise Exception("Bounds have already been set.")

        elif self.__is_active_chain():
            other = other if isinstance(other, highs_linear_expression) else highs_linear_expression(other)
            order = self.__get_chain(other, False)  # [self, LHS, inner, RHS, other], ignores None values

            # inner <= (self == RHS) <= other
            # LHS <= (self == inner) <= other
            if self.__is_equal_except_bounds(order[2]):
                return self.__compose_chain(*order[1:])

            # self <= (other == inner) <= RHS
            # self <= (other == LHS) <= inner
            elif other.__is_equal_except_bounds(order[1]):
                return self.__compose_chain(*order[:-1])

        self.__reset_chain(self, other, None)

        # self <= c
        if isinstance(other, (float, int)):
            copy = highs_linear_expression(self)
            copy.bounds = (-kHighsInf, float(other) - (copy.constant or 0.0))
            copy.constant = None
            return copy

        # self <= other
        elif isinstance(other, highs_linear_expression):
            if other.bounds is None:
                copy = highs_linear_expression()

                # prefer most idxs, constant, 'left'
                if len(other.idxs) > len(self.idxs) or (len(other.idxs) == len(self.idxs) and other.constant is None and self.constant is not None):
                    copy.idxs = other.idxs + self.idxs
                    copy.vals = other.vals + [-v for v in self.vals]
                    copy.bounds = (
                        (self.constant or 0.0) - (other.constant or 0.0),
                        kHighsInf,
                    )
                else:
                    copy.idxs = self.idxs + other.idxs
                    copy.vals = self.vals + [-v for v in other.vals]
                    copy.bounds = (
                        -kHighsInf,
                        (other.constant or 0.0) - (self.constant or 0.0),
                    )

                return copy
            else:
                raise Exception("Bounds have already been set.")

        # self <= x
        else:  # other is highs_var
            copy = highs_linear_expression()

            if len(self.idxs) == 0 or len(self.idxs) == 1 and self.constant is not None:
                copy.idxs = [other.index] + self.idxs
                copy.vals = [1.0] + [-v for v in self.vals]
                copy.bounds = ((self.constant or 0.0), kHighsInf)
            else:
                copy.idxs = self.idxs + [other.index]
                copy.vals = self.vals + [-1.0]
                copy.bounds = (-kHighsInf, -(self.constant or 0.0))

            return copy

    # other <= self
    def __ge__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        if self.bounds is not None:
            raise Exception("Bounds have already been set.")

        elif self.__is_active_chain():
            other = other if isinstance(other, highs_linear_expression) else highs_linear_expression(other)
            order = self.__get_chain(other, True)  # [other, LHS, inner, RHS, self], ignores None values

            # other <= (self == LHS) <= inner
            # other <= (self == inner) <= RHS
            if self.__is_equal_except_bounds(order[1]):
                return self.__compose_chain(*order[:-1])

            # LHS <= (other == inner) <= self
            # inner <= (other == RHS) <= self
            elif other.__is_equal_except_bounds(order[2]):
                return self.__compose_chain(*order[1:])

        self.__reset_chain(None, other, self)

        # c <= self
        if isinstance(other, (float, int)):
            copy = highs_linear_expression(self)
            copy.bounds = (float(other) - (self.constant or 0.0), kHighsInf)
            copy.constant = None
            return copy

        # other <= self
        elif isinstance(other, highs_linear_expression):
            if other.bounds is None:
                copy = highs_linear_expression()

                # prefer most idxs, constant, 'left'
                if len(self.idxs) > len(other.idxs) or (len(self.idxs) == len(other.idxs) and self.constant is None and other.constant is not None):
                    copy.idxs = self.idxs + other.idxs
                    copy.vals = self.vals + [-v for v in other.vals]
                    copy.bounds = (
                        (other.constant or 0.0) - (self.constant or 0.0),
                        kHighsInf,
                    )
                else:
                    copy.idxs = other.idxs + self.idxs
                    copy.vals = other.vals + [-v for v in self.vals]
                    copy.bounds = (
                        -kHighsInf,
                        (self.constant or 0.0) - (other.constant or 0.0),
                    )

                return copy
            else:
                raise Exception("Bounds have already been set.")

        # x <= self
        else:  # other is highs_var
            copy = highs_linear_expression()

            if len(self.idxs) > 1:
                copy.idxs = self.idxs + [other.index]
                copy.vals = self.vals + [-1.0]
                copy.bounds = (-(self.constant or 0.0), kHighsInf)
            else:
                copy.idxs = [other.index] + self.idxs
                copy.vals = [1.0] + [-v for v in self.vals]
                copy.bounds = (-kHighsInf, (self.constant or 0.0))

            return copy

    def __radd__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        return self + other

    # (LHS <= self <= RHS) + (LHS <= other <= RHS)
    def __add__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        copy = highs_linear_expression(self)
        copy += other
        return copy

    def __neg__(self):
        copy = highs_linear_expression()
        copy.idxs = list(self.idxs)
        copy.vals = [-v for v in self.vals]
        copy.constant = -self.constant if self.constant is not None else None

        if self.bounds is not None:
            copy.bounds = (-self.bounds[1], -self.bounds[0])

        return copy

    def __rmul__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        return self * other

    def __mul__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        copy = highs_linear_expression(self)
        copy *= other
        return copy

    # other - self
    def __rsub__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        copy = highs_linear_expression(other)
        copy -= self
        return copy

    def __sub__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        copy = highs_linear_expression(self)
        copy -= other
        return copy

    def unique_elements(self):
        """
        Collects unique variables and sums their corresponding values.  Keeps all values (including zeros).
        """
        # sort by groups for fast unique
        groups = np.asarray(self.idxs, dtype=np.int32)
        order = np.argsort(groups, kind="stable")
        groups = groups[order]

        # get unique groups
        index = np.ones(len(groups), dtype=bool)
        index[:-1] = groups[1:] != groups[:-1]

        if index.all():
            values = np.asarray(self.vals, dtype=np.float64)
            values = values[order]
            return groups, values
        else:
            values = np.asarray(self.vals, dtype=np.float64)
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
    def __iadd__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        if isinstance(other, highs_var):
            self.idxs.append(other.index)
            self.vals.append(1.0)
            return self

        elif isinstance(other, highs_linear_expression):
            if self.constant is not None and other.bounds is not None or self.bounds is not None and other.constant is not None:
                raise Exception("""Cannot add a bounded constraint to a constraint with a constant, i.e., (lb <= expr1 <= ub) + (expr2 + c). 
                    Unsure of your intent. Did you want: lb + c <= expr1 + expr2 <= ub + c?  Try: (lb <= expr1 <= ub) + (expr2 == c) instead.""")

            self.idxs.extend(other.idxs)
            self.vals.extend(other.vals)

            if self.constant is not None or other.constant is not None:
                self.constant = (self.constant or 0.0) + (other.constant or 0.0)

            # (l1 <= expr1 <= u1) + (l2 <= expr2 <= u2)  -->  l1 + l2 <= expr1 + expr2 <= u1 + u2
            if self.bounds is not None and other.bounds is not None:
                self.bounds = (self.bounds[0] + other.bounds[0], self.bounds[1] + other.bounds[1])

            # (expr1) + (lb <= expr2 <= ub)  -->  lb <= expr1 + expr2 <= ub
            elif self.bounds is None and other.bounds is not None:
                self.bounds = other.bounds

            return self

        else:
            if self.bounds is not None:
                raise Exception("""Cannot add a constant to a bounded constraint, i.e., (lb <= expr <= ub) + c. 
                    Unsure of your intent. Did you want: lb + c <= expr <= ub + c?  Try: (lb <= expr <= ub) + (highs_linear_expression() == c) instead.""")

            self.constant = float(other) + (self.constant or 0.0)
            return self

    def __isub__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        if isinstance(other, highs_var):
            self.idxs.append(other.index)
            self.vals.append(-1.0)
            return self

        elif isinstance(other, highs_linear_expression):
            if self.constant is not None and other.bounds is not None or self.bounds is not None and other.constant is not None:
                raise Exception("""Cannot subtract a bounded constraint to a constraint with a constant, i.e., (lb <= expr1 <= ub) - (expr2 + c). 
                    Unsure of your intent. Did you want: lb - c <= expr1 - expr2 <= ub - c?  Try: (lb <= expr1 <= ub) - (expr2 == c) instead.""")

            self.idxs.extend(other.idxs)
            self.vals.extend([-v for v in other.vals])

            if self.constant is not None or other.constant is not None:
                self.constant = (self.constant or 0.0) - (other.constant or 0.0)

            # (l1 <= expr1 <= u1) - (l2 <= expr2 <= u2)  -->  l1 - u2 <= expr1 - expr2 <= u1 - l2
            if self.bounds is not None and other.bounds is not None:
                self.bounds = (self.bounds[0] - other.bounds[1], self.bounds[1] - other.bounds[0])

            # expr1 - (lb <= expr2 <= ub)  -->  -ub <= expr1 - expr2 <= -lb
            elif self.bounds is None and other.bounds is not None:
                self.bounds = (-other.bounds[1], -other.bounds[0])

            return self

        else:
            if self.bounds is not None:
                raise Exception("""Cannot subtract a constant to a bounded constraint, i.e., (lb <= expr <= ub) - c. 
                    Unsure of your intent. Did you want: lb - c <= expr <= ub - c?  Try: (lb <= expr <= ub) - (highs_linear_expression() == c) instead.""")

            self.constant = (self.constant or 0.0) - float(other)
            return self

    def __imul__(self, other: Union[float, int, highs_var, highs_linear_expression]):
        if isinstance(other, (float, int)):
            scale = float(other)

        # other is a constant expression, so treat as a scalar
        elif isinstance(other, highs_linear_expression) and other.idxs == [] and other.constant is not None:
            scale = float(other.constant)
        else:
            scale = None

        if scale is not None:
            self.vals = [scale * v for v in self.vals]

            if self.constant is not None:
                self.constant *= scale

            if self.bounds is not None:
                # negative scale reverses bounds
                if scale >= 0:
                    self.bounds = (scale * self.bounds[0], scale * self.bounds[1])
                else:
                    self.bounds = (scale * self.bounds[1], scale * self.bounds[0])

            return self

        # self only has a constant, so treat as a scalar
        elif self.idxs == [] and self.constant is not None:
            scale = self.constant

            if isinstance(other, highs_linear_expression):
                self.idxs = other.idxs
                self.vals = [scale * v for v in other.vals]

                if other.constant is not None:
                    self.constant = other.constant * scale
                else:
                    self.constant = None

                if other.bounds is not None:
                    # negative scale reverses bounds
                    if scale >= 0:
                        self.bounds = (scale * other.bounds[0], scale * other.bounds[1])
                    else:
                        self.bounds = (scale * other.bounds[1], scale * other.bounds[0])

            elif isinstance(other, highs_var):
                self.idxs = [other.index]
                self.vals = [scale]
                self.constant = None

            else:
                raise Exception("Unexpected parameters.")

            return self

        elif isinstance(other, highs_var):
            raise Exception("Only linear expressions are allowed.")

        else:
            raise Exception("Unexpected parameters.")

    # The following is needed to support chained comparison, i.e., lb <= expr <= ub. This is interpreted
    # as '__bool__(lb <= expr) and (expr <= ub)'; returning (expr <= ub), since __bool__(lb <= expr) == True.
    #
    # We essentially want to "rewrite" this as '(lb <= expr) <= ub', while keeping the expr instance immutable.
    # As a slight hack, we can use a shared (thread local) object to keep track of the chain.
    #
    # Whenever we perform an inequality, we first check if the current expression is part of a chain.
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
        LHS = getattr(highs_linear_expression.__chain, "left", None)
        EXR = getattr(highs_linear_expression.__chain, "inner", None)
        RHS = getattr(highs_linear_expression.__chain, "right", None)

        highs_linear_expression.__chain.left = highs_linear_expression(LHS) if LHS is not None else None
        highs_linear_expression.__chain.inner = highs_linear_expression(EXR) if EXR is not None else None
        highs_linear_expression.__chain.right = highs_linear_expression(RHS) if RHS is not None else None
        return True

    def __is_equal_except_bounds(self, other: highs_linear_expression) -> bool:
        return self.idxs == other.idxs and self.vals == other.vals and self.constant == other.constant

    def __is_active_chain(self):
        return getattr(highs_linear_expression.__chain, "check", None) is not None

    def __reset_chain(
        self,
        LHS: Optional[Any] = None,
        inner: Optional[Any] = None,
        RHS: Optional[Any] = None,
    ) -> None:
        highs_linear_expression.__chain.check = None
        highs_linear_expression.__chain.left = LHS
        highs_linear_expression.__chain.inner = inner
        highs_linear_expression.__chain.right = RHS

    def __get_chain(self, other: highs_linear_expression, is_ge_than: bool):
        LHS = getattr(highs_linear_expression.__chain, "left", None)
        RHS = getattr(highs_linear_expression.__chain, "right", None)
        inner = getattr(highs_linear_expression.__chain, "inner", None)

        # assume (LHS is None) ^ (RHS is None) == 1, i.e., only LHS or RHS is set
        assert (LHS is None) ^ (RHS is None) == 1
        order = np.asarray([expr for expr in [other, LHS, inner, RHS, self] if expr is not None])

        if not is_ge_than:  # swap order for: self <= other
            order[0], order[-1] = order[-1], order[0]

        return order

    # left <= self <= right
    def __compose_chain(
        self,
        left: highs_linear_expression,
        inner: highs_linear_expression,
        right: highs_linear_expression,
    ):
        self.__reset_chain()
        assert not isinstance(inner, highs_linear_expression) or inner.bounds is None, "Bounds already set in chain comparison."

        # check to see if we have a valid chain, i.e., left.idxs "==" right.idxs
        # we can assume that left and right are both highs_linear_expression
        LHS_vars, LHS_vals = left.reduced_elements()
        RHS_vars, RHS_vals = right.reduced_elements()

        if not np.array_equal(LHS_vars, RHS_vars) or not np.array_equal(LHS_vals, RHS_vals):
            raise Exception("Mismatched variables in chain comparison.")

        if len(LHS_vars) > len(inner.idxs):
            copy = highs_linear_expression()
            copy.idxs = LHS_vars.tolist() + inner.idxs
            copy.vals = LHS_vals.tolist() + [-v for v in inner.vals]
            copy.bounds = (
                (inner.constant or 0.0) - (right.constant or 0.0),
                (inner.constant or 0.0) - (left.constant or 0.0),
            )
        else:
            copy = highs_linear_expression(inner)
            copy.idxs.extend(LHS_vars)
            copy.vals.extend([-v for v in LHS_vals])
            copy.bounds = (
                (left.constant or 0.0) - (copy.constant or 0.0),
                (right.constant or 0.0) - (copy.constant or 0.0),
            )
            copy.constant = None

        return copy


def qsum(
    items: Iterable[Union[highs_var, highs_linear_expression]],
    initial: Optional[Union[float, int, highs_var, highs_linear_expression]] = None,
):
    """
    Performs a faster sum for highs_linear_expressions.

    Args:
        items: A collection of highs_linear_expressions or highs_vars to be summed.
    """
    return Highs.qsum(items, initial)
