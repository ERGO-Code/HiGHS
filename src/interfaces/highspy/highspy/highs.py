from .highs_bindings import (
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
    HighsSparseMatrix,
    HighsLp,
    HighsHessian,
    HighsModel,
    HighsSolution,
    HighsBasis,
    HighsInfo,
    HighsOptions,
    _Highs,
    kHighsInf,
    HIGHS_VERSION_MAJOR,
    HIGHS_VERSION_MINOR,
    HIGHS_VERSION_PATCH,
)


class Highs(_Highs):
    def __init__(self):
        super().__init__()
        self._log_user_callback_tuple = CallbackTuple()

    def setLogCallback(self, func
#                       , deprecated// V2.0 remove
                       ):
        self._log_user_callback_tuple.callback = func
        self._log_user_callback_tuple.log_deprecated = log_deprecated
        super().setLogCallback(self._log_user_callback_tuple)
