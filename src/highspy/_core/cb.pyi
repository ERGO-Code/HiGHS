"""
Callback interface submodule
"""

from __future__ import annotations
import highspy._core
import typing

__all__ = [
    "HighsCallbackDataIn",
    "HighsCallbackDataOut",
    "HighsCallbackType",
    "kCallbackIpmInterrupt",
    "kCallbackLogging",
    "kCallbackMax",
    "kCallbackMin",
    "kCallbackMipDefineLazyConstraints",
    "kCallbackMipGetCutPool",
    "kCallbackMipImprovingSolution",
    "kCallbackMipInterrupt",
    "kCallbackMipLogging",
    "kCallbackMipSolution",
    "kCallbackSimplexInterrupt",
    "kNumCallbackType",
]

class HighsCallbackDataIn:
    user_interrupt: int
    def __init__(self) -> None: ...

class HighsCallbackDataOut:
    ipm_iteration_count: int
    log_type: int
    mip_dual_bound: float
    mip_gap: float
    mip_node_count: int
    mip_primal_bound: float
    objective_function_value: float
    pdlp_iteration_count: int
    running_time: float
    simplex_iteration_count: int
    def __init__(self) -> None: ...
    @property
    def mip_solution(self) -> highspy._core.readonly_ptr_wrapper_double: ...

class HighsCallbackType:
    """
    Members:
      kCallbackMin
      kCallbackLogging
      kCallbackSimplexInterrupt
      kCallbackIpmInterrupt
      kCallbackMipSolution
      kCallbackMipImprovingSolution
      kCallbackMipLogging
      kCallbackMipInterrupt
      kCallbackMipGetCutPool
      kCallbackMipDefineLazyConstraints
      kCallbackMax
      kNumCallbackType
    """

    __members__: typing.ClassVar[dict[str, HighsCallbackType]]
    kCallbackIpmInterrupt: typing.ClassVar[HighsCallbackType]
    kCallbackLogging: typing.ClassVar[HighsCallbackType]
    kCallbackMax: typing.ClassVar[HighsCallbackType]
    kCallbackMin: typing.ClassVar[HighsCallbackType]
    kCallbackMipDefineLazyConstraints: typing.ClassVar[HighsCallbackType]
    kCallbackMipGetCutPool: typing.ClassVar[HighsCallbackType]
    kCallbackMipImprovingSolution: typing.ClassVar[HighsCallbackType]
    kCallbackMipInterrupt: typing.ClassVar[HighsCallbackType]
    kCallbackMipLogging: typing.ClassVar[HighsCallbackType]
    kCallbackMipSolution: typing.ClassVar[HighsCallbackType]
    kCallbackSimplexInterrupt: typing.ClassVar[HighsCallbackType]
    kNumCallbackType: typing.ClassVar[HighsCallbackType]
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

kCallbackIpmInterrupt: HighsCallbackType
kCallbackLogging: HighsCallbackType
kCallbackMax: HighsCallbackType
kCallbackMin: HighsCallbackType
kCallbackMipDefineLazyConstraints: HighsCallbackType
kCallbackMipGetCutPool: HighsCallbackType
kCallbackMipImprovingSolution: HighsCallbackType
kCallbackMipInterrupt: HighsCallbackType
kCallbackMipLogging: HighsCallbackType
kCallbackMipSolution: HighsCallbackType
kCallbackSimplexInterrupt: HighsCallbackType
kNumCallbackType: HighsCallbackType
