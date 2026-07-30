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
    IisStrategy,
    IisBoundStatus,
    IisStatus,
    HighsSparseMatrix,
    HighsLp,
    HighsHessian,
    HighsModel,
    HighsInfo,
    HighsOptions,
    _Highs,  # type: ignore
    HighsIis,
    HighsSolution,
    HighsObjectiveSolution,
    HighsBasis,
    HighsRangingRecord,
    HighsRanging,
    kHighsInf,
    kHighsIInf,
    kHighsUndefined,
    HighsLinearObjective,
    HIGHS_VERSION_MAJOR,
    HIGHS_VERSION_MINOR,
    HIGHS_VERSION_PATCH,
    simplex_constants,  # type: ignore
    cb,  # type: ignore
    kSolutionStatusNone,
    kSolutionStatusInfeasible,
    kSolutionStatusFeasible,
    kBasisValidityInvalid,
    kBasisValidityValid,
)

from .highs import Highs, HighsCallbackEvent, HighsCallback, HighspyArray, highs_var, highs_cons, highs_linear_expression

__all__ = [
    "__doc__",
    "ObjSense",
    "MatrixFormat",
    "HessianFormat",
    "SolutionStatus",
    "BasisValidity",
    "HighsModelStatus",
    "HighsPresolveStatus",
    "HighsBasisStatus",
    "HighsVarType",
    "HighsOptionType",
    "HighsInfoType",
    "HighsStatus",
    "HighsLogType",
    "IisStrategy",
    "IisBoundStatus",
    "IisStatus",
    "HighsSparseMatrix",
    "HighsLp",
    "HighsHessian",
    "HighsModel",
    "HighsInfo",
    "HighsOptions",
    "_Highs",
    "Highs",
    "HighsIis",
    "HighsSolution",
    "HighsObjectiveSolution",
    "HighsBasis",
    "HighsRangingRecord",
    "HighsRanging",
    "kHighsInf",
    "kHighsIInf",
    "kHighsUndefined",
    "HighsLinearObjective",
    "HIGHS_VERSION_MAJOR",
    "HIGHS_VERSION_MINOR",
    "HIGHS_VERSION_PATCH",
    "simplex_constants",
    "cb",
    "kSolutionStatusNone",
    "kSolutionStatusInfeasible",
    "kSolutionStatusFeasible",
    "kBasisValidityInvalid",
    "kBasisValidityValid",
    "HighsCallbackEvent",
    "HighsCallback",
    "HighspyArray",
    "highs_var",
    "highs_cons",
    "highs_linear_expression",
]
