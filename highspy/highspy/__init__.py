import os as _os
import sys as _sys

if _sys.platform == "win32" and _os.path.exists(_os.path.join(_os.path.dirname(__file__), "cudalin.dll")):
    # Python 3.8+ no longer searches PATH when loading extension DLLs.
    # cudalin.dll is present only in GPU builds; register the CUDA bin directory
    # so cuda*.dll etc. are found when _core is loaded.
    _cuda_path = _os.environ.get("CUDA_PATH") or _os.environ.get("CUDA_HOME")
    if _cuda_path:
        _cuda_bin = _os.path.join(_cuda_path, "bin")
        if _os.path.isdir(_cuda_bin):
            _os.add_dll_directory(_cuda_bin)
    del _cuda_path
del _os, _sys

from ._core import (
    HIGHS_VERSION_MAJOR,
    HIGHS_VERSION_MINOR,
    HIGHS_VERSION_PATCH,
    BasisValidity,
    HessianFormat,
    HighsBasis,
    HighsBasisStatus,
    HighsHessian,
    HighsIis,
    HighsInfo,
    HighsInfoType,
    HighsLinearObjective,
    HighsLogType,
    HighsLp,
    HighsModel,
    HighsModelStatus,
    HighsObjectiveSolution,
    HighsOptions,
    HighsOptionType,
    HighsPresolveStatus,
    HighsRanging,
    HighsRangingRecord,
    HighsSolution,
    HighsSparseMatrix,
    HighsStatus,
    HighsVarType,
    IisBoundStatus,
    IisStatus,
    IisStrategy,
    MatrixFormat,
    ObjSense,
    SolutionStatus,
    _Highs,  # type: ignore
    cb,  # type: ignore
    kBasisValidityInvalid,
    kBasisValidityValid,
    kHighsIInf,
    kHighsInf,
    kHighsUndefined,
    kSolutionStatusFeasible,
    kSolutionStatusInfeasible,
    kSolutionStatusNone,
    simplex_constants,  # type: ignore
)
from .highs import (
    Highs,
    HighsCallback,
    HighsCallbackEvent,
    HighsError,
    HighspyArray,
    HighsStatusError,
    highs_cons,
    highs_linear_expression,
    highs_var,
    qsum,
)

__all__ = [
    "HIGHS_VERSION_MAJOR",
    "HIGHS_VERSION_MINOR",
    "HIGHS_VERSION_PATCH",
    "BasisValidity",
    "HessianFormat",
    "Highs",
    "HighsBasis",
    "HighsBasisStatus",
    "HighsCallback",
    "HighsCallbackEvent",
    "HighsError",
    "HighsHessian",
    "HighsIis",
    "HighsInfo",
    "HighsInfoType",
    "HighsLinearObjective",
    "HighsLogType",
    "HighsLp",
    "HighsModel",
    "HighsModelStatus",
    "HighsObjectiveSolution",
    "HighsOptionType",
    "HighsOptions",
    "HighsPresolveStatus",
    "HighsRanging",
    "HighsRangingRecord",
    "HighsSolution",
    "HighsSparseMatrix",
    "HighsStatus",
    "HighsStatusError",
    "HighsVarType",
    "HighspyArray",
    "IisBoundStatus",
    "IisStatus",
    "IisStrategy",
    "MatrixFormat",
    "ObjSense",
    "SolutionStatus",
    "_Highs",
    "__doc__",
    "cb",
    "highs_cons",
    "highs_linear_expression",
    "highs_var",
    "kBasisValidityInvalid",
    "kBasisValidityValid",
    "kHighsIInf",
    "kHighsInf",
    "kHighsUndefined",
    "kSolutionStatusFeasible",
    "kSolutionStatusInfeasible",
    "kSolutionStatusNone",
    "qsum",
    "simplex_constants",
]
