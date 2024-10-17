from __future__ import annotations
import numpy
import typing
from . import cb
from . import simplex_constants

__all__ = [
    "BasisValidity",
    "HIGHS_VERSION_MAJOR",
    "HIGHS_VERSION_MINOR",
    "HIGHS_VERSION_PATCH",
    "HessianFormat",
    "HighsBasis",
    "HighsBasisStatus",
    "HighsHessian",
    "HighsIis",
    "HighsIisInfo",
    "HighsInfo",
    "HighsInfoType",
    "HighsLogType",
    "HighsLp",
    "HighsLpMods",
    "HighsModel",
    "HighsModelStatus",
    "HighsObjectiveSolution",
    "HighsOptionType",
    "HighsOptions",
    "HighsPresolveStatus",
    "HighsRanging",
    "HighsRangingRecord",
    "HighsScale",
    "HighsSolution",
    "HighsSparseMatrix",
    "HighsStatus",
    "HighsVarType",
    "IisBoundStatus",
    "IisStrategy",
    "MatrixFormat",
    "ObjSense",
    "SolutionStatus",
    "cb",
    "kBasisValidityInvalid",
    "kBasisValidityValid",
    "kHighsIInf",
    "kHighsInf",
    "kSolutionStatusFeasible",
    "kSolutionStatusInfeasible",
    "kSolutionStatusNone",
    "readonly_ptr_wrapper_double",
    "simplex_constants",
]

class BasisValidity:
    """
    Members:
      kBasisValidityInvalid
      kBasisValidityValid
    """

    __members__: typing.ClassVar[dict[str, BasisValidity]]
    kBasisValidityInvalid: typing.ClassVar[BasisValidity]
    kBasisValidityValid: typing.ClassVar[BasisValidity]

    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class HessianFormat:
    """
    Members:
      kTriangular
      kSquare
    """

    __members__: typing.ClassVar[dict[str, HessianFormat]]
    kSquare: typing.ClassVar[HessianFormat]
    kTriangular: typing.ClassVar[HessianFormat]

    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class HighsBasis:
    alien: bool
    col_status: list[HighsBasisStatus]
    debug_id: int
    debug_origin_name: str
    debug_update_count: int
    row_status: list[HighsBasisStatus]
    valid: bool
    was_alien: bool

    def __init__(self) -> None: ...

class HighsBasisStatus:
    """
    Members:
      kLower
      kBasic
      kUpper
      kZero
      kNonbasic
    """

    __members__: typing.ClassVar[dict[str, HighsBasisStatus]]
    kBasic: typing.ClassVar[HighsBasisStatus]
    kLower: typing.ClassVar[HighsBasisStatus]
    kNonbasic: typing.ClassVar[HighsBasisStatus]
    kUpper: typing.ClassVar[HighsBasisStatus]
    kZero: typing.ClassVar[HighsBasisStatus]

    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class HighsHessian:
    dim_: int
    format_: HessianFormat
    index_: list[int]
    start_: list[int]
    value_: list[float]

    def __init__(self) -> None: ...

class HighsIis:
    col_bound: list[int]
    col_index: list[int]
    info: list[HighsIisInfo]
    row_bound: list[int]
    row_index: list[int]
    strategy: int
    valid: bool

    def __init__(self) -> None: ...
    def invalidate(self) -> None: ...

class HighsIisInfo:
    simplex_iterations: int
    simplex_time: float

    def __init__(self) -> None: ...

class HighsInfo:
    basis_validity: int
    crossover_iteration_count: int
    dual_solution_status: int
    ipm_iteration_count: int
    max_complementarity_violation: float
    max_dual_infeasibility: float
    max_integrality_violation: float
    max_primal_infeasibility: float
    mip_dual_bound: float
    mip_gap: float
    mip_node_count: int
    num_dual_infeasibilities: int
    num_primal_infeasibilities: int
    objective_function_value: float
    pdlp_iteration_count: int
    primal_solution_status: int
    qp_iteration_count: int
    simplex_iteration_count: int
    sum_complementarity_violations: float
    sum_dual_infeasibilities: float
    sum_primal_infeasibilities: float
    valid: bool

    def __init__(self) -> None: ...

class HighsInfoType:
    """
    Members:
      kInt64
      kInt
      kDouble
    """

    __members__: typing.ClassVar[dict[str, HighsInfoType]]
    kDouble: typing.ClassVar[HighsInfoType]
    kInt: typing.ClassVar[HighsInfoType]
    kInt64: typing.ClassVar[HighsInfoType]

    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class HighsLogType:
    """
    Members:
      kInfo
      kDetailed
      kVerbose
      kWarning
      kError
    """

    __members__: typing.ClassVar[dict[str, HighsLogType]]
    kDetailed: typing.ClassVar[HighsLogType]
    kError: typing.ClassVar[HighsLogType]
    kInfo: typing.ClassVar[HighsLogType]
    kVerbose: typing.ClassVar[HighsLogType]
    kWarning: typing.ClassVar[HighsLogType]

    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class HighsLp:
    a_matrix_: HighsSparseMatrix
    col_cost_: list[float]
    col_lower_: list[float]
    col_names_: list[str]
    col_upper_: list[float]
    integrality_: list[HighsVarType]
    is_moved_: bool
    is_scaled_: bool
    model_name_: str
    mods_: HighsLpMods
    num_col_: int
    num_row_: int
    offset_: float
    row_lower_: list[float]
    row_names_: list[str]
    row_upper_: list[float]
    scale_: HighsScale
    sense_: ObjSense

    def __init__(self) -> None: ...

class HighsLpMods:
    pass

class HighsModel:
    hessian_: HighsHessian
    lp_: HighsLp

    def __init__(self) -> None: ...

class HighsModelStatus:
    """
    Members:
      kNotset
      kLoadError
      kModelError
      kPresolveError
      kSolveError
      kPostsolveError
      kModelEmpty
      kOptimal
      kInfeasible
      kUnboundedOrInfeasible
      kUnbounded
      kObjectiveBound
      kObjectiveTarget
      kTimeLimit
      kIterationLimit
      kUnknown
      kSolutionLimit
      kInterrupt
      kMemoryLimit
    """

    __members__: typing.ClassVar[dict[str, HighsModelStatus]]
    kInfeasible: typing.ClassVar[HighsModelStatus]
    kInterrupt: typing.ClassVar[HighsModelStatus]
    kIterationLimit: typing.ClassVar[HighsModelStatus]
    kLoadError: typing.ClassVar[HighsModelStatus]
    kMemoryLimit: typing.ClassVar[HighsModelStatus]
    kModelEmpty: typing.ClassVar[HighsModelStatus]
    kModelError: typing.ClassVar[HighsModelStatus]
    kNotset: typing.ClassVar[HighsModelStatus]
    kObjectiveBound: typing.ClassVar[HighsModelStatus]
    kObjectiveTarget: typing.ClassVar[HighsModelStatus]
    kOptimal: typing.ClassVar[HighsModelStatus]
    kPostsolveError: typing.ClassVar[HighsModelStatus]
    kPresolveError: typing.ClassVar[HighsModelStatus]
    kSolutionLimit: typing.ClassVar[HighsModelStatus]
    kSolveError: typing.ClassVar[HighsModelStatus]
    kTimeLimit: typing.ClassVar[HighsModelStatus]
    kUnbounded: typing.ClassVar[HighsModelStatus]
    kUnboundedOrInfeasible: typing.ClassVar[HighsModelStatus]
    kUnknown: typing.ClassVar[HighsModelStatus]

    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class HighsObjectiveSolution:
    col_value: list[float]
    objective: float

    def __init__(self) -> None: ...

class HighsOptionType:
    """
    Members:
      kBool
      kInt
      kDouble
      kString
    """

    __members__: typing.ClassVar[dict[str, HighsOptionType]]
    kBool: typing.ClassVar[HighsOptionType]
    kDouble: typing.ClassVar[HighsOptionType]
    kInt: typing.ClassVar[HighsOptionType]
    kString: typing.ClassVar[HighsOptionType]

    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class HighsOptions:
    allow_unbounded_or_infeasible: bool
    allowed_matrix_scale_factor: int
    dual_feasibility_tolerance: float
    highs_analysis_level: int
    highs_debug_level: int
    infinite_bound: float
    infinite_cost: float
    ipm_iteration_limit: int
    ipm_optimality_tolerance: float
    ipx_dualize_strategy: int
    large_matrix_value: float
    log_dev_level: int
    log_file: str
    log_to_console: bool
    mip_abs_gap: float
    mip_detect_symmetry: bool
    mip_feasibility_tolerance: float
    mip_heuristic_effort: float
    mip_lp_age_limit: int
    mip_max_improving_sols: int
    mip_max_leaves: int
    mip_max_nodes: int
    mip_max_stall_nodes: int
    mip_min_cliquetable_entries_for_parallelism: int
    mip_min_logging_interval: float
    mip_pool_age_limit: int
    mip_pool_soft_limit: int
    mip_pscost_minreliable: int
    mip_rel_gap: float
    mip_report_level: int
    objective_bound: float
    objective_target: float
    output_flag: bool
    parallel: str
    presolve: str
    primal_feasibility_tolerance: float
    random_seed: int
    ranging: str
    run_crossover: str
    simplex_crash_strategy: int
    simplex_dual_edge_weight_strategy: int
    simplex_dualize_strategy: int
    simplex_iteration_limit: int
    simplex_max_concurrency: int
    simplex_min_concurrency: int
    simplex_permute_strategy: int
    simplex_price_strategy: int
    simplex_primal_edge_weight_strategy: int
    simplex_scale_strategy: int
    simplex_strategy: int
    simplex_update_limit: int
    small_matrix_value: float
    solution_file: str
    solver: str
    threads: int
    time_limit: float
    write_model_file: str
    write_model_to_file: bool
    write_solution_style: int
    write_solution_to_file: bool

    def __init__(self) -> None: ...

class HighsPresolveStatus:
    """
    Members:
      kNotPresolved
      kNotReduced
      kInfeasible
      kUnboundedOrInfeasible
      kReduced
      kReducedToEmpty
      kTimeout
      kNullError
      kOptionsError
    """

    __members__: typing.ClassVar[dict[str, HighsPresolveStatus]]
    kInfeasible: typing.ClassVar[HighsPresolveStatus]
    kNotPresolved: typing.ClassVar[HighsPresolveStatus]
    kNotReduced: typing.ClassVar[HighsPresolveStatus]
    kNullError: typing.ClassVar[HighsPresolveStatus]
    kOptionsError: typing.ClassVar[HighsPresolveStatus]
    kReduced: typing.ClassVar[HighsPresolveStatus]
    kReducedToEmpty: typing.ClassVar[HighsPresolveStatus]
    kTimeout: typing.ClassVar[HighsPresolveStatus]
    kUnboundedOrInfeasible: typing.ClassVar[HighsPresolveStatus]

    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class HighsRanging:
    col_bound_dn: HighsRangingRecord
    col_bound_up: HighsRangingRecord
    col_cost_dn: HighsRangingRecord
    col_cost_up: HighsRangingRecord
    row_bound_dn: HighsRangingRecord
    row_bound_up: HighsRangingRecord
    valid: bool

    def __init__(self) -> None: ...

class HighsRangingRecord:
    in_var_: list[int]
    objective_: list[float]
    ou_var_: list[int]
    value_: list[float]

    def __init__(self) -> None: ...

class HighsScale:
    pass

class HighsSolution:
    col_dual: list[float]
    col_value: list[float]
    dual_valid: bool
    row_dual: list[float]
    row_value: list[float]
    value_valid: bool

    def __init__(self) -> None: ...

class HighsSparseMatrix:
    format_: MatrixFormat
    index_: list[int]
    num_col_: int
    num_row_: int
    p_end_: list[int]
    start_: list[int]
    value_: list[float]

    def __init__(self) -> None: ...

class HighsStatus:
    """
    Members:
      kError
      kOk
      kWarning
    """

    __members__: typing.ClassVar[dict[str, HighsStatus]]
    kError: typing.ClassVar[HighsStatus]
    kOk: typing.ClassVar[HighsStatus]
    kWarning: typing.ClassVar[HighsStatus]

    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class HighsVarType:
    """
    Members:
      kContinuous
      kInteger
      kSemiContinuous
      kSemiInteger
    """

    __members__: typing.ClassVar[dict[str, HighsVarType]]
    kContinuous: typing.ClassVar[HighsVarType]
    kInteger: typing.ClassVar[HighsVarType]
    kSemiContinuous: typing.ClassVar[HighsVarType]
    kSemiInteger: typing.ClassVar[HighsVarType]

    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class IisBoundStatus:
    """
    Members:
      kIisBoundStatusDropped
      kIisBoundStatusNull
      kIisBoundStatusFree
      kIisBoundStatusLower
      kIisBoundStatusUpper
      kIisBoundStatusBoxed
    """

    __members__: typing.ClassVar[dict[str, IisBoundStatus]]
    kIisBoundStatusBoxed: typing.ClassVar[IisBoundStatus]
    kIisBoundStatusDropped: typing.ClassVar[IisBoundStatus]
    kIisBoundStatusFree: typing.ClassVar[IisBoundStatus]
    kIisBoundStatusLower: typing.ClassVar[IisBoundStatus]
    kIisBoundStatusNull: typing.ClassVar[IisBoundStatus]
    kIisBoundStatusUpper: typing.ClassVar[IisBoundStatus]

    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class IisStrategy:
    """
    Members:
      kIisStrategyMin
      kIisStrategyFromLpRowPriority
      kIisStrategyFromLpColPriority
      kIisStrategyMax
    """

    __members__: typing.ClassVar[dict[str, IisStrategy]]
    kIisStrategyFromLpColPriority: typing.ClassVar[IisStrategy]
    kIisStrategyFromLpRowPriority: typing.ClassVar[IisStrategy]
    kIisStrategyMax: typing.ClassVar[IisStrategy]
    kIisStrategyMin: typing.ClassVar[IisStrategy]

    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class MatrixFormat:
    """
    Members:
      kColwise
      kRowwise
      kRowwisePartitioned
    """

    __members__: typing.ClassVar[dict[str, MatrixFormat]]
    kColwise: typing.ClassVar[MatrixFormat]
    kRowwise: typing.ClassVar[MatrixFormat]
    kRowwisePartitioned: typing.ClassVar[MatrixFormat]

    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class ObjSense:
    """
    Members:
      kMinimize
      kMaximize
    """

    __members__: typing.ClassVar[dict[str, ObjSense]]
    kMaximize: typing.ClassVar[ObjSense]
    kMinimize: typing.ClassVar[ObjSense]

    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class SolutionStatus:
    """
    Members:
      kSolutionStatusNone
      kSolutionStatusInfeasible
      kSolutionStatusFeasible
    """

    __members__: typing.ClassVar[dict[str, SolutionStatus]]
    kSolutionStatusFeasible: typing.ClassVar[SolutionStatus]
    kSolutionStatusInfeasible: typing.ClassVar[SolutionStatus]
    kSolutionStatusNone: typing.ClassVar[SolutionStatus]

    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class _Highs:
    def __init__(self) -> None: ...
    def addCol(
        self,
        cost: float,
        lower_bound: float,
        upper_bound: float,
        num_elements: int,
        indices: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        values: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
    ) -> HighsStatus: ...
    def addCols(
        self,
        num_cols: int,
        costs: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        lower_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        upper_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        num_elements: int,
        starts: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        indices: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        values: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
    ) -> HighsStatus: ...
    def addRow(
        self,
        lower_bound: float,
        upper_bound: float,
        num_elements: int,
        indices: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        values: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
    ) -> HighsStatus: ...
    def addRows(
        self,
        num_rows: int,
        lower_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        upper_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        num_elements: int,
        starts: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        indices: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        values: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
    ) -> HighsStatus: ...
    def addVar(self, lower_bound: float, upper_bound: float) -> HighsStatus: ...
    def addVars(
        self,
        num_vars: int,
        lower_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        upper_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
    ) -> HighsStatus: ...
    def basisStatusToString(self, basis_status: HighsBasisStatus) -> str: ...
    def basisValidityToString(self, validity: int) -> str: ...
    def changeCoeff(self, row: int, col: int, value: float) -> HighsStatus: ...
    def changeColBounds(self, col: int, lower_bound: float, upper_bound: float) -> HighsStatus: ...
    def changeColCost(self, col: int, cost: float) -> HighsStatus: ...
    def changeColIntegrality(self, col: int, var_type: HighsVarType) -> HighsStatus: ...
    def changeColsBounds(
        self,
        num_cols: int,
        cols: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        lower_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        upper_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
    ) -> HighsStatus: ...
    def changeColsCost(
        self,
        num_cols: int,
        cols: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        costs: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
    ) -> HighsStatus: ...
    def changeColsIntegrality(
        self,
        num_cols: int,
        cols: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        var_types: numpy.ndarray[typing.Any, numpy.dtype[numpy.uint8]],
    ) -> HighsStatus: ...
    def changeObjectiveOffset(self, offset: float) -> HighsStatus: ...
    def changeObjectiveSense(self, sense: ObjSense) -> HighsStatus: ...
    def changeRowBounds(self, row: int, lower_bound: float, upper_bound: float) -> HighsStatus: ...
    def clear(self) -> HighsStatus: ...
    def clearModel(self) -> HighsStatus: ...
    def clearSolver(self) -> HighsStatus: ...
    def crossover(self, solution: HighsSolution) -> HighsStatus: ...
    def deleteCols(self, num_cols: int, cols: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]]) -> HighsStatus: ...
    def deleteRows(self, num_rows: int, rows: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]]) -> HighsStatus: ...
    def deleteVars(self, num_vars: int, vars: list[int]) -> HighsStatus: ...
    def feasibilityRelaxation(
        self,
        global_lower_penalty: float,
        global_upper_penalty: float,
        global_rhs_penalty: float,
        local_lower_penalty: typing.Any = ...,
        local_upper_penalty: typing.Any = ...,
        local_rhs_penalty: typing.Any = ...,
    ) -> HighsStatus: ...
    def getBasis(self) -> HighsBasis: ...
    def getCol(self, col: int) -> tuple[HighsStatus, float, float, float, int]: ...
    def getColByName(self, name: str) -> tuple[HighsStatus, int]: ...
    def getColEntries(
        self, col: int
    ) -> tuple[HighsStatus, numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]], numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]]: ...
    def getColName(self, col: int) -> tuple[HighsStatus, str]: ...
    def getCols(
        self, num_cols: int, cols: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]]
    ) -> tuple[
        HighsStatus,
        int,
        numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        int,
    ]: ...
    def getColsEntries(
        self, num_cols: int, cols: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]]
    ) -> tuple[
        HighsStatus,
        numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
    ]: ...
    def getHessianNumNz(self) -> int: ...
    def getIis(self, iis: HighsIis) -> HighsStatus: ...
    def getInfinity(self) -> float: ...
    def getInfo(self) -> HighsInfo: ...
    def getInfoType(self, info: str) -> tuple[HighsStatus, HighsInfoType]: ...
    def getInfoValue(self, info: str) -> tuple[HighsStatus, object]: ...
    def getLp(self) -> HighsLp: ...
    def getModel(self) -> HighsModel: ...
    def getModelPresolveStatus(self) -> HighsPresolveStatus: ...
    def getModelStatus(self) -> HighsModelStatus: ...
    def getNumCol(self) -> int: ...
    def getNumNz(self) -> int: ...
    def getNumRow(self) -> int: ...
    def getObjectiveOffset(self) -> tuple[HighsStatus, float]: ...
    def getObjectiveSense(self) -> tuple[HighsStatus, ObjSense]: ...
    def getObjectiveValue(self) -> float: ...
    def getOptionType(self, name: str) -> tuple[HighsStatus, HighsOptionType]: ...
    def getOptionValue(self, name: str) -> tuple[HighsStatus, object]: ...
    def getOptions(self) -> HighsOptions: ...
    def getPresolvedLp(self) -> HighsLp: ...
    def getRanging(self) -> tuple[HighsStatus, HighsRanging]: ...
    def getRow(self, row_index: int) -> tuple[HighsStatus, float, float, int]: ...
    def getRowByName(self, row_name: str) -> tuple[HighsStatus, int]: ...
    def getRowEntries(
        self, row_index: int
    ) -> tuple[HighsStatus, numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]], numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]]: ...
    def getRowName(self, row_index: int) -> tuple[HighsStatus, str]: ...
    def getRows(
        self, num_rows: int, row_indices: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]]
    ) -> tuple[
        HighsStatus,
        int,
        numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        int,
    ]: ...
    def getRowsEntries(
        self, num_rows: int, row_indices: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]]
    ) -> tuple[
        HighsStatus,
        numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
    ]: ...
    def getRunTime(self) -> float: ...
    def getSavedMipSolutions(self) -> list[HighsObjectiveSolution]: ...
    def getSolution(self) -> HighsSolution: ...
    def githash(self) -> str: ...
    def modelStatusToString(self, model_status: HighsModelStatus) -> str: ...
    def passColName(self, col_index: int, col_name: str) -> HighsStatus: ...
    @typing.overload
    def passHessian(self, hessian: HighsHessian) -> HighsStatus: ...
    @typing.overload
    def passHessian(
        self,
        dim: int,
        nnz: int,
        format: int,
        q_start: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        q_index: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        q_value: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
    ) -> HighsStatus: ...
    @typing.overload
    def passModel(self, model: HighsModel) -> HighsStatus: ...
    @typing.overload
    def passModel(
        self,
        num_cols: int,
        num_rows: int,
        nnz: int,
        qnnz: int,
        a_format: int,
        q_format: int,
        sense: int,
        offset: float,
        col_costs: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        col_lower_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        col_upper_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        row_lower_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        row_upper_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        a_start: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        a_index: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        a_value: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        q_start: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        q_index: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        q_value: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        integrality: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
    ) -> HighsStatus: ...
    @typing.overload
    def passModel(self, lp: HighsLp) -> HighsStatus: ...
    @typing.overload
    def passModel(
        self,
        num_cols: int,
        num_rows: int,
        num_entries: int,
        a_format: int,
        sense: int,
        offset: float,
        col_costs: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        col_lower_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        col_upper_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        row_lower_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        row_upper_bounds: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        a_start: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        a_index: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        a_value: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
        integrality: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
    ) -> HighsStatus: ...
    def passOptions(self, options: HighsOptions) -> HighsStatus: ...
    def passRowName(self, row_index: int, row_name: str) -> HighsStatus: ...
    @typing.overload
    def postsolve(self, solution: HighsSolution, basis: HighsBasis) -> HighsStatus: ...
    @typing.overload
    def postsolve(self, solution: HighsSolution) -> HighsStatus: ...
    def presolve(self) -> HighsStatus: ...
    def readBasis(self, file_path: str) -> HighsStatus: ...
    def readModel(self, file_path: str) -> HighsStatus: ...
    def readOptions(self, file_path: str) -> HighsStatus: ...
    def readSolution(self, file_path: str, style: int) -> HighsStatus: ...
    @staticmethod
    def resetGlobalScheduler(blocking: bool) -> None: ...
    def resetOptions(self) -> HighsStatus: ...
    def run(self) -> HighsStatus: ...
    @typing.overload
    def setBasis(self, basis: HighsBasis) -> HighsStatus: ...
    @typing.overload
    def setBasis(self) -> HighsStatus: ...
    def setCallback(
        self,
        callback: typing.Callable[
            [cb.HighsCallbackType, str, cb.HighsCallbackDataOut, cb.HighsCallbackDataIn, typing.Any],
            None,
        ]
        | None,
        callback_data: typing.Any | None,
    ) -> HighsStatus: ...
    @typing.overload
    def setOptionValue(self, option_name: str, option_value: bool) -> HighsStatus: ...
    @typing.overload
    def setOptionValue(self, option_name: str, option_value: int) -> HighsStatus: ...
    @typing.overload
    def setOptionValue(self, option_name: str, option_value: float) -> HighsStatus: ...
    @typing.overload
    def setOptionValue(self, option_name: str, option_value: str) -> HighsStatus: ...
    @typing.overload
    def setSolution(self, solution: HighsSolution) -> HighsStatus: ...
    @typing.overload
    def setSolution(
        self,
        num_entries: int,
        indices: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]],
        values: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]],
    ) -> HighsStatus: ...
    def solutionStatusToString(self, solution_status: int) -> str: ...
    def startCallback(self, callback_type: cb.HighsCallbackType) -> HighsStatus: ...
    def startCallbackInt(self, callback_type: int) -> HighsStatus: ...
    def stopCallback(self, callback_type: cb.HighsCallbackType) -> HighsStatus: ...
    def stopCallbackInt(self, callback_type: int) -> HighsStatus: ...
    def version(self) -> str: ...
    def versionMajor(self) -> int: ...
    def versionMinor(self) -> int: ...
    def versionPatch(self) -> int: ...
    def writeBasis(self, file_path: str) -> HighsStatus: ...
    def writeInfo(self, file_path: str) -> HighsStatus: ...
    def writeModel(self, file_path: str) -> HighsStatus: ...
    def writeOptions(self, file_path: str) -> HighsStatus: ...
    def writePresolvedModel(self, file_path: str) -> HighsStatus: ...
    def writeSolution(self, file_path: str, style: int) -> HighsStatus: ...

class readonly_ptr_wrapper_double:
    def __bool__(self) -> bool: ...
    def __getitem__(self, idx: int) -> float: ...
    def __init__(self, ptr) -> None: ...
    def to_array(self, size: int) -> numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]: ...

HIGHS_VERSION_MAJOR: int
HIGHS_VERSION_MINOR: int
HIGHS_VERSION_PATCH: int
kBasisValidityInvalid: BasisValidity
kBasisValidityValid: BasisValidity
kHighsIInf: int
kHighsInf: float
kSolutionStatusFeasible: SolutionStatus
kSolutionStatusInfeasible: SolutionStatus
kSolutionStatusNone: SolutionStatus
