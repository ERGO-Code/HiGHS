from .highs_bindings import (
    # enum classes
    ObjSense,
    MatrixFormat,
    HessianFormat,
    SolutionStatus,
    BasisValidity,
    HighsModelStatus,
    HighsBasisStatus,
    HighsVarType,
    HighsOptionType,
    HighsInfoType,
    HighsStatus,
    HighsLogType,
    CallbackTuple,
    # classes
    HighsSparseMatrix,
    HighsLp,
    HighsHessian,
    HighsModel,
    HighsSolution,
    HighsBasis,
    HighsInfo,
    HighsOptions,
    _Highs,
    # Constants
    kHighsInf,
    kHighsIInf,
)


class Highs(_Highs):
    def __init__(self):
        super().__init__()
        self._log_callback_tuple = CallbackTuple()

    def setLogCallback(self, func, callback_data):
        self._log_callback_tuple.callback = func
        self._log_callback_tuple.callback_data = callback_data
        super().setLogCallback(self._log_callback_tuple)
