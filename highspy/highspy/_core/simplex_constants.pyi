"""
Submodule for simplex constants
"""

from __future__ import annotations
import typing

__all__ = [
    "EdgeWeightMode",
    "SimplexEdgeWeightStrategy",
    "SimplexNlaOperation",
    "SimplexPivotalRowRefinementStrategy",
    "SimplexPriceStrategy",
    "SimplexPrimalCorrectionStrategy",
    "SimplexSolvePhase",
    "SimplexStrategy",
    "SimplexUnscaledSolutionStrategy",
    "kNumSimplexNlaOperation",
    "kSimplexEdgeWeightStrategyChoose",
    "kSimplexEdgeWeightStrategyDantzig",
    "kSimplexEdgeWeightStrategyDevex",
    "kSimplexEdgeWeightStrategyMax",
    "kSimplexEdgeWeightStrategyMin",
    "kSimplexEdgeWeightStrategySteepestEdge",
    "kSimplexInfeasibilityProofRefinementAlsoScaledLp",
    "kSimplexInfeasibilityProofRefinementMax",
    "kSimplexInfeasibilityProofRefinementMin",
    "kSimplexInfeasibilityProofRefinementNo",
    "kSimplexInfeasibilityProofRefinementUnscaledLp",
    "kSimplexNlaBtranBasicFeasibilityChange",
    "kSimplexNlaBtranEp",
    "kSimplexNlaBtranFull",
    "kSimplexNlaBtranPse",
    "kSimplexNlaFtran",
    "kSimplexNlaFtranBfrt",
    "kSimplexNlaFtranDse",
    "kSimplexNlaNull",
    "kSimplexNlaPriceAp",
    "kSimplexNlaPriceFull",
    "kSimplexPriceStrategyCol",
    "kSimplexPriceStrategyMax",
    "kSimplexPriceStrategyMin",
    "kSimplexPriceStrategyRow",
    "kSimplexPriceStrategyRowSwitch",
    "kSimplexPriceStrategyRowSwitchColSwitch",
    "kSimplexPrimalCorrectionStrategyAlways",
    "kSimplexPrimalCorrectionStrategyInRebuild",
    "kSimplexPrimalCorrectionStrategyNone",
    "kSimplexStrategyChoose",
    "kSimplexStrategyDual",
    "kSimplexStrategyDualMulti",
    "kSimplexStrategyDualPlain",
    "kSimplexStrategyDualTasks",
    "kSimplexStrategyMax",
    "kSimplexStrategyMin",
    "kSimplexStrategyNum",
    "kSimplexStrategyPrimal",
    "kSimplexUnscaledSolutionStrategyDirect",
    "kSimplexUnscaledSolutionStrategyMax",
    "kSimplexUnscaledSolutionStrategyMin",
    "kSimplexUnscaledSolutionStrategyNone",
    "kSimplexUnscaledSolutionStrategyNum",
    "kSimplexUnscaledSolutionStrategyRefine",
    "kSolvePhase1",
    "kSolvePhase2",
    "kSolvePhaseError",
    "kSolvePhaseExit",
    "kSolvePhaseMax",
    "kSolvePhaseMin",
    "kSolvePhaseOptimal",
    "kSolvePhaseOptimalCleanup",
    "kSolvePhasePrimalInfeasibleCleanup",
    "kSolvePhaseTabooBasis",
    "kSolvePhaseUnknown",
]

class EdgeWeightMode:
    """
    Members:
      kDantzig
      kDevex
      kSteepestEdge
      kCount
    """

    __members__: typing.ClassVar[dict[str, EdgeWeightMode]]
    kCount: typing.ClassVar[EdgeWeightMode]
    kDantzig: typing.ClassVar[EdgeWeightMode]
    kDevex: typing.ClassVar[EdgeWeightMode]
    kSteepestEdge: typing.ClassVar[EdgeWeightMode]
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

class SimplexEdgeWeightStrategy:
    """
    Members:
      kSimplexEdgeWeightStrategyMin
      kSimplexEdgeWeightStrategyChoose
      kSimplexEdgeWeightStrategyDantzig
      kSimplexEdgeWeightStrategyDevex
      kSimplexEdgeWeightStrategySteepestEdge
      kSimplexEdgeWeightStrategyMax
    """

    __members__: typing.ClassVar[dict[str, SimplexEdgeWeightStrategy]]
    kSimplexEdgeWeightStrategyChoose: typing.ClassVar[SimplexEdgeWeightStrategy]
    kSimplexEdgeWeightStrategyDantzig: typing.ClassVar[SimplexEdgeWeightStrategy]
    kSimplexEdgeWeightStrategyDevex: typing.ClassVar[SimplexEdgeWeightStrategy]
    kSimplexEdgeWeightStrategyMax: typing.ClassVar[SimplexEdgeWeightStrategy]
    kSimplexEdgeWeightStrategyMin: typing.ClassVar[SimplexEdgeWeightStrategy]
    kSimplexEdgeWeightStrategySteepestEdge: typing.ClassVar[SimplexEdgeWeightStrategy]
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

class SimplexNlaOperation:
    """
    Members:
      kSimplexNlaNull
      kSimplexNlaBtranFull
      kSimplexNlaPriceFull
      kSimplexNlaBtranBasicFeasibilityChange
      kSimplexNlaBtranEp
      kSimplexNlaPriceAp
      kSimplexNlaFtran
      kSimplexNlaFtranBfrt
      kSimplexNlaFtranDse
      kSimplexNlaBtranPse
      kNumSimplexNlaOperation
    """

    __members__: typing.ClassVar[dict[str, SimplexNlaOperation]]
    kNumSimplexNlaOperation: typing.ClassVar[SimplexNlaOperation]
    kSimplexNlaBtranBasicFeasibilityChange: typing.ClassVar[SimplexNlaOperation]
    kSimplexNlaBtranEp: typing.ClassVar[SimplexNlaOperation]
    kSimplexNlaBtranFull: typing.ClassVar[SimplexNlaOperation]
    kSimplexNlaBtranPse: typing.ClassVar[SimplexNlaOperation]
    kSimplexNlaFtran: typing.ClassVar[SimplexNlaOperation]
    kSimplexNlaFtranBfrt: typing.ClassVar[SimplexNlaOperation]
    kSimplexNlaFtranDse: typing.ClassVar[SimplexNlaOperation]
    kSimplexNlaNull: typing.ClassVar[SimplexNlaOperation]
    kSimplexNlaPriceAp: typing.ClassVar[SimplexNlaOperation]
    kSimplexNlaPriceFull: typing.ClassVar[SimplexNlaOperation]
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

class SimplexPivotalRowRefinementStrategy:
    """
    Members:
      kSimplexInfeasibilityProofRefinementMin
      kSimplexInfeasibilityProofRefinementNo
      kSimplexInfeasibilityProofRefinementUnscaledLp
      kSimplexInfeasibilityProofRefinementAlsoScaledLp
      kSimplexInfeasibilityProofRefinementMax
    """

    __members__: typing.ClassVar[dict[str, SimplexPivotalRowRefinementStrategy]]
    kSimplexInfeasibilityProofRefinementAlsoScaledLp: typing.ClassVar[
        SimplexPivotalRowRefinementStrategy
    ]
    kSimplexInfeasibilityProofRefinementMax: typing.ClassVar[
        SimplexPivotalRowRefinementStrategy
    ]
    kSimplexInfeasibilityProofRefinementMin: typing.ClassVar[
        SimplexPivotalRowRefinementStrategy
    ]
    kSimplexInfeasibilityProofRefinementNo: typing.ClassVar[
        SimplexPivotalRowRefinementStrategy
    ]
    kSimplexInfeasibilityProofRefinementUnscaledLp: typing.ClassVar[
        SimplexPivotalRowRefinementStrategy
    ]
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

class SimplexPriceStrategy:
    """
    Members:
      kSimplexPriceStrategyMin
      kSimplexPriceStrategyCol
      kSimplexPriceStrategyRow
      kSimplexPriceStrategyRowSwitch
      kSimplexPriceStrategyRowSwitchColSwitch
      kSimplexPriceStrategyMax
    """

    __members__: typing.ClassVar[dict[str, SimplexPriceStrategy]]
    kSimplexPriceStrategyCol: typing.ClassVar[SimplexPriceStrategy]
    kSimplexPriceStrategyMax: typing.ClassVar[SimplexPriceStrategy]
    kSimplexPriceStrategyMin: typing.ClassVar[SimplexPriceStrategy]
    kSimplexPriceStrategyRow: typing.ClassVar[SimplexPriceStrategy]
    kSimplexPriceStrategyRowSwitch: typing.ClassVar[SimplexPriceStrategy]
    kSimplexPriceStrategyRowSwitchColSwitch: typing.ClassVar[SimplexPriceStrategy]
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

class SimplexPrimalCorrectionStrategy:
    """
    Members:
      kSimplexPrimalCorrectionStrategyNone
      kSimplexPrimalCorrectionStrategyInRebuild
      kSimplexPrimalCorrectionStrategyAlways
    """

    __members__: typing.ClassVar[dict[str, SimplexPrimalCorrectionStrategy]]
    kSimplexPrimalCorrectionStrategyAlways: typing.ClassVar[
        SimplexPrimalCorrectionStrategy
    ]
    kSimplexPrimalCorrectionStrategyInRebuild: typing.ClassVar[
        SimplexPrimalCorrectionStrategy
    ]
    kSimplexPrimalCorrectionStrategyNone: typing.ClassVar[
        SimplexPrimalCorrectionStrategy
    ]
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

class SimplexSolvePhase:
    """
    Members:
      kSolvePhaseMin
      kSolvePhaseError
      kSolvePhaseExit
      kSolvePhaseUnknown
      kSolvePhaseOptimal
      kSolvePhase1
      kSolvePhase2
      kSolvePhasePrimalInfeasibleCleanup
      kSolvePhaseOptimalCleanup
      kSolvePhaseTabooBasis
      kSolvePhaseMax
    """

    __members__: typing.ClassVar[dict[str, SimplexSolvePhase]]
    kSolvePhase1: typing.ClassVar[SimplexSolvePhase]
    kSolvePhase2: typing.ClassVar[SimplexSolvePhase]
    kSolvePhaseError: typing.ClassVar[SimplexSolvePhase]
    kSolvePhaseExit: typing.ClassVar[SimplexSolvePhase]
    kSolvePhaseMax: typing.ClassVar[SimplexSolvePhase]
    kSolvePhaseMin: typing.ClassVar[SimplexSolvePhase]
    kSolvePhaseOptimal: typing.ClassVar[SimplexSolvePhase]
    kSolvePhaseOptimalCleanup: typing.ClassVar[SimplexSolvePhase]
    kSolvePhasePrimalInfeasibleCleanup: typing.ClassVar[SimplexSolvePhase]
    kSolvePhaseTabooBasis: typing.ClassVar[SimplexSolvePhase]
    kSolvePhaseUnknown: typing.ClassVar[SimplexSolvePhase]
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

class SimplexStrategy:
    """
    Members:
      kSimplexStrategyMin
      kSimplexStrategyChoose
      kSimplexStrategyDual
      kSimplexStrategyDualPlain
      kSimplexStrategyDualTasks
      kSimplexStrategyDualMulti
      kSimplexStrategyPrimal
      kSimplexStrategyMax
      kSimplexStrategyNum
    """

    __members__: typing.ClassVar[dict[str, SimplexStrategy]]
    kSimplexStrategyChoose: typing.ClassVar[SimplexStrategy]
    kSimplexStrategyDual: typing.ClassVar[SimplexStrategy]
    kSimplexStrategyDualMulti: typing.ClassVar[SimplexStrategy]
    kSimplexStrategyDualPlain: typing.ClassVar[SimplexStrategy]
    kSimplexStrategyDualTasks: typing.ClassVar[SimplexStrategy]
    kSimplexStrategyMax: typing.ClassVar[SimplexStrategy]
    kSimplexStrategyMin: typing.ClassVar[SimplexStrategy]
    kSimplexStrategyNum: typing.ClassVar[SimplexStrategy]
    kSimplexStrategyPrimal: typing.ClassVar[SimplexStrategy]
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

class SimplexUnscaledSolutionStrategy:
    """
    Members:
      kSimplexUnscaledSolutionStrategyMin
      kSimplexUnscaledSolutionStrategyNone
      kSimplexUnscaledSolutionStrategyRefine
      kSimplexUnscaledSolutionStrategyDirect
      kSimplexUnscaledSolutionStrategyMax
      kSimplexUnscaledSolutionStrategyNum
    """

    __members__: typing.ClassVar[dict[str, SimplexUnscaledSolutionStrategy]]
    kSimplexUnscaledSolutionStrategyDirect: typing.ClassVar[
        SimplexUnscaledSolutionStrategy
    ]
    kSimplexUnscaledSolutionStrategyMax: typing.ClassVar[
        SimplexUnscaledSolutionStrategy
    ]
    kSimplexUnscaledSolutionStrategyMin: typing.ClassVar[
        SimplexUnscaledSolutionStrategy
    ]
    kSimplexUnscaledSolutionStrategyNone: typing.ClassVar[
        SimplexUnscaledSolutionStrategy
    ]
    kSimplexUnscaledSolutionStrategyNum: typing.ClassVar[
        SimplexUnscaledSolutionStrategy
    ]
    kSimplexUnscaledSolutionStrategyRefine: typing.ClassVar[
        SimplexUnscaledSolutionStrategy
    ]
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

kNumSimplexNlaOperation: SimplexNlaOperation
kSimplexEdgeWeightStrategyChoose: SimplexEdgeWeightStrategy
kSimplexEdgeWeightStrategyDantzig: SimplexEdgeWeightStrategy
kSimplexEdgeWeightStrategyDevex: SimplexEdgeWeightStrategy
kSimplexEdgeWeightStrategyMax: SimplexEdgeWeightStrategy
kSimplexEdgeWeightStrategyMin: SimplexEdgeWeightStrategy
kSimplexEdgeWeightStrategySteepestEdge: SimplexEdgeWeightStrategy
kSimplexInfeasibilityProofRefinementAlsoScaledLp: SimplexPivotalRowRefinementStrategy
kSimplexInfeasibilityProofRefinementMax: SimplexPivotalRowRefinementStrategy
kSimplexInfeasibilityProofRefinementMin: SimplexPivotalRowRefinementStrategy
kSimplexInfeasibilityProofRefinementNo: SimplexPivotalRowRefinementStrategy
kSimplexInfeasibilityProofRefinementUnscaledLp: SimplexPivotalRowRefinementStrategy
kSimplexNlaBtranBasicFeasibilityChange: SimplexNlaOperation
kSimplexNlaBtranEp: SimplexNlaOperation
kSimplexNlaBtranFull: SimplexNlaOperation
kSimplexNlaBtranPse: SimplexNlaOperation
kSimplexNlaFtran: SimplexNlaOperation
kSimplexNlaFtranBfrt: SimplexNlaOperation
kSimplexNlaFtranDse: SimplexNlaOperation
kSimplexNlaNull: SimplexNlaOperation
kSimplexNlaPriceAp: SimplexNlaOperation
kSimplexNlaPriceFull: SimplexNlaOperation
kSimplexPriceStrategyCol: SimplexPriceStrategy
kSimplexPriceStrategyMax: SimplexPriceStrategy
kSimplexPriceStrategyMin: SimplexPriceStrategy
kSimplexPriceStrategyRow: SimplexPriceStrategy
kSimplexPriceStrategyRowSwitch: SimplexPriceStrategy
kSimplexPriceStrategyRowSwitchColSwitch: SimplexPriceStrategy
kSimplexPrimalCorrectionStrategyAlways: SimplexPrimalCorrectionStrategy
kSimplexPrimalCorrectionStrategyInRebuild: SimplexPrimalCorrectionStrategy
kSimplexPrimalCorrectionStrategyNone: SimplexPrimalCorrectionStrategy
kSimplexStrategyChoose: SimplexStrategy
kSimplexStrategyDual: SimplexStrategy
kSimplexStrategyDualMulti: SimplexStrategy
kSimplexStrategyDualPlain: SimplexStrategy
kSimplexStrategyDualTasks: SimplexStrategy
kSimplexStrategyMax: SimplexStrategy
kSimplexStrategyMin: SimplexStrategy
kSimplexStrategyNum: SimplexStrategy
kSimplexStrategyPrimal: SimplexStrategy
kSimplexUnscaledSolutionStrategyDirect: SimplexUnscaledSolutionStrategy
kSimplexUnscaledSolutionStrategyMax: SimplexUnscaledSolutionStrategy
kSimplexUnscaledSolutionStrategyMin: SimplexUnscaledSolutionStrategy
kSimplexUnscaledSolutionStrategyNone: SimplexUnscaledSolutionStrategy
kSimplexUnscaledSolutionStrategyNum: SimplexUnscaledSolutionStrategy
kSimplexUnscaledSolutionStrategyRefine: SimplexUnscaledSolutionStrategy
kSolvePhase1: SimplexSolvePhase
kSolvePhase2: SimplexSolvePhase
kSolvePhaseError: SimplexSolvePhase
kSolvePhaseExit: SimplexSolvePhase
kSolvePhaseMax: SimplexSolvePhase
kSolvePhaseMin: SimplexSolvePhase
kSolvePhaseOptimal: SimplexSolvePhase
kSolvePhaseOptimalCleanup: SimplexSolvePhase
kSolvePhasePrimalInfeasibleCleanup: SimplexSolvePhase
kSolvePhaseTabooBasis: SimplexSolvePhase
kSolvePhaseUnknown: SimplexSolvePhase
