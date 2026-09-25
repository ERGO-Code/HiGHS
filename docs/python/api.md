# API Reference

The public `highspy` package re-exports the main solver class, model data
structures, enums, constants, and callback helpers.

## Top-Level Package

```{eval-rst}
.. automodule:: highspy
   :inherited-members: object, Exception, BaseException, ndarray
```

## High-Level Python Interface

`highspy.Highs` inherits most of its solver API (`readModel`, `run`,
`getSolution`, `setOptionValue`, and so on) from the pybind11-bound `_Highs`
base class, so those inherited members are included below too.

```{eval-rst}
.. automodule:: highspy.highs
   :inherited-members: object, Exception, BaseException, ndarray
```

## Core Solver Bindings

The {py:mod}`highspy._core` module contains the pybind11 bindings to the
HiGHS C++ API. In an installed build these entries include the signatures
exposed by the compiled extension.

```{eval-rst}
.. automodule:: highspy._core
```

## Callback Bindings

```{eval-rst}
.. automodule:: highspy._core.cb
```

## Simplex Constants

```{eval-rst}
.. automodule:: highspy._core.simplex_constants
```
