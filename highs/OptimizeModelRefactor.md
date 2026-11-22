# Refactoring of `Highs::calledOptimizeModel()`

## Summary

- `calledOptimizeModel()` was a single function of 926 lines, with 34
  `return` statements (29 of them through `returnFromOptimizeModel()`). It is now a short driver that calls one
  function per phase. The LP part is a pipeline of `lp*` functions.
- The phases share their state through a `Context` struct, instead of
  through locals.
- Option changes go through a new `HighsUndoStack`, which undoes them
  on every return path, including exceptions. Each function that
  changes options has its own local one, so changes don't outlive the
  function that makes them.
- All modifications ("mods") of the incumbent model are applied in one
  place and undone in one place. Bounds with small inconsistencies
  are still corrected permanently, as before.
- `runPresolve()` and `runPostsolve()` take the `PresolveComponent`
  and the LP as arguments, so that later LP phases (for example LP
  folding) can run them on an LP other than `model_.lp_`.
- There is now a single exit, `returnFromOptimizeModel(context)`, once
  the solve is marked as running.

Most of the diff is code moved without changes. The behaviour changes
are listed in [Behaviour changes](#behaviour-changes), divided into
fixes and changes in behaviour.

## Files changed

| File | Change |
|---|---|
| `highs/lp_data/Highs.cpp` | Split of `calledOptimizeModel()`, mods, single exit, `runPresolve`/`runPostsolve` overloads, `callRunPostsolve` error path, `run()` unscaling |
| `highs/Highs.h` | `Context`, `LpState`, declarations of the new private methods, flags `optimize_model_running_` and `model_mods_applied_`, overloads |
| `highs/util/HighsUndoStack.h` | New |
| `cmake/sources.cmake` | Lists `util/HighsUndoStack.h` |
| `highs/lp_data/HighsSolve.cpp` | The `SOLVE_CATCH_CALL` macros replaced by the function templates `solveCatchCall()`, `solveCatchCallLpQp()` and `solveCatchCallMip()`, the last replacing the inline `try`/`catch` around the MIP solver. They don't log: a caught exception is stored in a `SolveExceptionRecord` and the caller logs it with `reportSolveException()`. `HighsUndoStack` for option changes |
| `highs/lp_data/HighsInterface.cpp` | `HighsUndoStack` for option changes, `has_infinite_cost_` kept up to date in multi-objective solves |
| `highs/simplex/HApp.h` | `HighsUndoStack` for option changes |
| `highs/HighsRun.md` | Describes where the mods are now applied, and drops the ToDo for moving them |

## Before

```mermaid
flowchart TD
    subgraph called["calledOptimizeModel() - one function, 926 lines"]
        direction TB
        a1["Debug level, name hashes, warm start,<br/>re-entry check"] --> a2["handleInfCost()<br/>(mods)"]
        a2 --> a3["exactResizeModel()"]
        a3 --> a4["completeSolutionFromDiscreteAssignment()<br/>(MIP with a solution)"]
        a4 --> a5["Reset info and timer;<br/>empty model, bounds, assessLp"]
        a5 --> a6["assessSemiVariables()<br/>(mods)"]
        a6 --> a7{"Problem type"}
        a7 -->|QP| a8["relaxSemiVariables() or MIQP error,<br/>convexity check,<br/>try { callSolveQp() } catch"]
        a7 -->|MIP| a9["relaxSemiVariables() and fall through to LP,<br/>or callSolveMip()"]
        a7 -->|LP| a10["iCrash (inline), basis set-up,<br/>timedSolveLp lambda"]
        a9 -.->|relaxation| a10
        a10 --> a11["No presolve: solve"]
        a10 --> a12["runPresolve(), then switch over<br/>11 presolve statuses (solve or exit)"]
        a12 --> a13["Postsolve and simplex hot start<br/>in one block"]
        a13 --> a14["PDLP clean-up"]
        a11 --> a15["KKT check, report times"]
        a14 --> a15
    end
    called --> r["returnFromOptimizeModel(status, undo_mods)<br/>reached from 29 places"]
```

## After: from `optimizeHighs()` to the problem type

```mermaid
flowchart TD
    optimizeHighs["optimizeHighs()"] -->|multiple objectives| multi["multiobjectiveSolve()"]
    optimizeHighs -->|single objective| optimizeModel["optimizeModel()<br/>initializeMultiThreading(), profiling"]
    multi -->|once per objective| optimizeModel

    optimizeModel --> tryCatch["optimizeModelTryCatch()<br/>creates Context"]
    optimizeLp["optimizeLp()<br/>(LP-only entry point)"] --> tryCatch
    tryCatch --> called["calledOptimizeModel(context)"]
    tryCatch -.->|exception| catch["handleCatch:<br/>clearSolver() keeping model status,<br/>optimizeModelUndoMods()"]

    called --> preChecks["optimizeModelPreChecks()<br/>re-entry check, debug level,<br/>name hashes, warm start, exactResizeModel(),<br/>model_.assess() assert"]
    preChecks -->|"called while running"| refuse["return kError,<br/>nothing touched"]
    preChecks --> applyMods["optimizeModelApplyMods()<br/>infinite costs, semi-variables"]
    applyMods --> complete["optimizeModelCompleteSolution()<br/>MIP with a (partial) solution:<br/>completeSolutionFromDiscreteAssignment(),<br/>may call optimizeModel() again"]
    applyMods -->|error| early["returnFromOptimizeModelBeforeRunning()<br/>invalidate data, undo mods"]
    complete -->|error| early
    complete --> running["optimize_model_running_ = true"]
    running --> prepare["optimizeModelPrepare()<br/>reset info and timer; empty model,<br/>bounds, debug assessLp"]

    prepare -->|QP| qp["optimizeModelQp()<br/>MIQP and convexity checks, callSolveQp()"]
    prepare -->|MIP| mip["optimizeModelMip()"]
    prepare -->|LP| lp["optimizeModelLp()<br/>(see below)"]
    mip -->|solve_relaxation| relax["optimizeModelLp()<br/>on the LP relaxation"]
    mip -->|otherwise| callMip["callSolveMip()"]

    prepare -->|early exit| returnFrom
    qp --> returnFrom["returnFromOptimizeModel(context)<br/>checks, optimize_model_running_ = false,<br/>optimizeModelUndoMods(), returnFromHighs()"]
    callMip --> returnFrom
    relax --> returnFrom
    lp --> returnFrom
```

Any change to `options_` during the solve is made through a
`HighsUndoStack`: a local one in the function that makes the change,
or `context.options_undo` for a change that must last until the end
of the solve. Anything still recorded in `context.options_undo` is
undone when the `Context` is destroyed at the end of
`optimizeModelTryCatch()`.

## After: `optimizeModelLp()`

```mermaid
flowchart TD
    setup["lpSetup(context)<br/>iCrash (lpICrash), basis from solution"] --> usePresolve{"lpUsePresolve(context, incumbent_lp)"}

    usePresolve -->|"no: useful basis,<br/>presolve off,<br/>or no matrix"| unreduced["lpSolveUnreduced(context, incumbent_lp)"]
    unreduced --> finish

    usePresolve -->|yes| presolve["lpPresolve(context, presolve_, lp)<br/>runPresolve(presolve_, incumbent_lp)"]
    presolve --> handle{"lpHandlePresolveStatus()"}
    handle -->|"infeasible, timeout,<br/>out of memory, error"| stop["return<br/>(skips lpFinish)"]
    handle -->|unbounded or infeasible| uoi["lpResolveUnboundedOrInfeasible(context, model_.lp_)<br/>primal simplex on the original LP"]
    uoi --> stop
    handle -->|reduced| prepareReduced["lpPrepareReducedLp(context, reduced LP)<br/>lp now points to the reduced LP"]
    handle -->|reduced to empty| empty["lpSetReducedToEmptySolution()"]
    handle -->|not presolved / not reduced| solveReduced

    prepareReduced --> solveReduced["lpSolveReduced(context, status, *lp)"]
    empty --> solveReduced
    solveReduced --> postsolve["lpPostsolve(context, presolve_, status, incumbent_lp)<br/>runPostsolve(presolve_, incumbent_lp)<br/>sets lp_state.postsolve_recovered_basis"]
    postsolve --> cleanup["lpCleanupOriginal(context, incumbent_lp)<br/>simplex hot start if a basis was recovered,<br/>PDLP clean-up (disabled, as before)"]
    cleanup --> finish["lpFinish(context, incumbent_lp)<br/>KKT check, lpReportTimes(), final status"]
```

Every function that can stop the solve returns `true` to do so, with
`context.status` set (it stays `kOk` when the stop isn't an error, for
example infeasibility detected by presolve), and the caller returns
straight away. Every LP
solve goes through `lpTimedSolve()`, which replaces the
`timedSolveLp` lambda.

## New mechanisms

### `Context` and `LpState`

`Highs::Context` is created in `optimizeModelTryCatch()` and passed
down by reference. It holds the status, the timings and the postsolve
iteration count that used to be locals, `min_highs_debug_level`,
`undo_mods` (see [Mods](#mods)),
`options_undo` (only for option changes that last the whole solve) and
`lp_state`.

`Highs::LpState` holds the values shared by the `lp*` functions: the
debug logging set-up, `solver_will_use_basis`,
`may_require_basis_postsolve`, `no_incumbent_lp_solution_or_basis`,
the pivot threshold from the presolved solve, the PDLP iteration count
of the presolved solve and `postsolve_recovered_basis`.

Postsolve was split from the simplex hot start and the PDLP clean-up.
The copy of `info_` that `latest` kept for the PDLP clean-up
(`presolved_lp_info`) is replaced by
`presolved_lp_pdlp_iteration_count`, the only field that
`tryPdlpCleanup()` read. It is now set in `lpPostsolve()`, where the
copy was taken, rather than after solving the reduced LP.
`postsolve_recovered_basis` is new: it replaces the
nesting of the hot start inside the branch of postsolve that recovered
a basis. Nothing ran between setting the recovered basis and the hot
start, and `LpState` is created fresh for each call, so the flag can't
leak from one solve to the next.

The other place where existing code was cut up is the switch over the
presolve status. Its cases are now spread over
`lpHandlePresolveStatus()` (the statuses that stop the solve),
`lpPrepareReducedLp()` and `lpSetReducedToEmptySolution()` (called
from `lpPresolve()`) and `lpSolveReduced()` (the solves). For
`kReduced` and `kReducedToEmpty`, the preparation and the solve are
consecutive, as before. These two splits deserve the closest look.

### `HighsUndoStack`

[util/HighsUndoStack.h](util/HighsUndoStack.h) records the old value
of a variable before it is changed, and restores it later. It replaces
saving and restoring values by hand, including copies of the whole
`HighsOptions`:

```cpp
HighsUndoStack options_undo;
options_undo.set(options_.solver, kSimplexString);
...
options_undo.reset();
```

The destructor undoes anything still recorded, so a change is also
undone on early returns and exceptions.

- **Local stacks** are used wherever options are changed and restored
  in the same function: `lpTimedSolve()`, `lpSolveReduced()`,
  `lpCleanupOriginal()` (one each for the debug check, the hot start
  and the PDLP clean-up), `lpResolveUnboundedOrInfeasible()`,
  `completeSolutionFromDiscreteAssignment()`, and in
  `HighsInterface.cpp`, `HighsSolve.cpp` and `simplex/HApp.h`. The
  explicit `reset()` is where the old code restored the values, so the
  order of restores on the normal paths is unchanged.
- **`context.options_undo`** is only for a change that must last until
  the end of the solve. The only one is raising `highs_debug_level` in
  `optimizeModelPreChecks()` (see B10).

The old code sometimes restored a whole copy of `HighsOptions`
(`options_ = save_options`). The new code restores only the fields it
changed. I found no other option written on those paths: the writes in
`HApp.h`, `HEkkDebug.cpp` and `solveLpWithSelection()` (through its
`RestoreOptions` struct, in `HighsSolve.cpp`) restore their own
values.

### Single exit and `optimize_model_running_`

- The new flag `optimize_model_running_` replaces
  `called_return_from_optimize_model`.
  - It is set after `optimizeModelCompleteSolution()`, which may call
    `optimizeModel()` itself.
  - It is cleared in `returnFromOptimizeModel()` and in `handleCatch`.
  - It is checked first in `optimizeModelPreChecks()`, before
    anything is changed.
- Once the flag is set, the only exit is
  `returnFromOptimizeModel(context)`, which reads the status from the
  context.
- Errors before the flag is set go through
  `returnFromOptimizeModelBeforeRunning()`. It invalidates info, run
  data, solution and basis, sets the model status and undoes the mods.
  The refusal of a call made while a solve is running doesn't use it,
  since that would destroy the running solve's data.
- The assert at the start of `returnFromOptimizeModel()` is removed,
  because the single exit makes it redundant.
- The check of `called_return_from_optimize_model` in
  `returnFromHighs()` is also removed. It caught API calls that end in
  `returnFromHighs()` while a solve was running, which in practice
  means calls from a callback. A future PR will refuse calls to
  methods that modify the `Highs` instance from within a callback,
  which covers this case in release builds too.

### Mods

- `optimizeModelApplyMods()` makes all changes to the incumbent model,
  in this order: infinite costs, then assessing or (with
  `solve_relaxation`) relaxing semi-variables. Bounds with small
  inconsistencies are corrected permanently by
  `processedBoundsFindsInfeasibility()` in `optimizeModelPrepare()`,
  as before, so they are not a mod.
- Everything needed to undo a mod is recorded in `model_.lp_.mods_`
  (`HighsLpMods`). The flags that say whether mods are applied live
  elsewhere: `HighsLp::has_infinite_cost_`, `model_mods_applied_` and
  `context.undo_mods`.
- `optimizeModelUndoMods()` undoes them. It is called from
  `returnFromOptimizeModel()`, `returnFromOptimizeModelBeforeRunning()`
  and `handleCatch`. It calls `restoreInfCost()`, then
  `HighsLp::unapplyMods()`, which restores the semi-variables.
- This undo order is not the reverse of the apply order, since
  infinite costs are applied first but restored first. See
  [Known limitations](#known-limitations-not-introduced-here).
- The new flag `model_mods_applied_` stops the nested `optimizeModel()`
  call made by `completeSolutionFromDiscreteAssignment()` from applying
  or undoing them a second time.

### Presolve and postsolve on any LP

These overloads are new:

```cpp
runPresolve(PresolveComponent&, HighsLp&, force_lp_presolve, force_presolve)
runPostsolve(PresolveComponent&, const HighsLp&)
```

The old signatures are wrappers that pass `presolve_` and
`model_.lp_`, so `presolve()`, `postsolve()` and the presolve getters
are unchanged. `presolve_log_` and `model_presolve_status_` are only
set when the component is `presolve_`. `model_.isEmpty()` and
`model_.isMip()` became `lp.isEmpty()` and `lp.isMip()`, which are the
same for `model_.lp_`. (`HighsLp::isEmpty()` tests `<= 0` rather than
`== 0`, which is the same for valid dimensions.)

## Behaviour changes

### Fixes

| # | Change | Before | After |
|---|---|---|---|
| F1 | MIP primal feasibility tolerance (`solveMip`) | The "restore" assigned the MIP value again, so after any MIP solve the user's `primal_feasibility_tolerance` was left at `mip_feasibility_tolerance` (1e-7 became 1e-6 with defaults), and later LP solves used it | Restored |
| F2 | Infinite costs and a user-supplied MIP solution | If the nested call from `completeSolutionFromDiscreteAssignment()` made semi-variable mods, its undo also ran `restoreInfCost()` on the outer call's record, so the outer MIP was solved with the infinite costs back in place | The nested call leaves mods alone |
| F3 | Re-entry check | Ran after `clearSolver()` (without warm start) and other changes, so a refused nested call still cleared the running solve's data | Runs first, changes nothing |
| F4 | Early returns that skipped `returnFromOptimizeModel()`: unknown iCrash strategy, `cleanBounds` error | Left the flag stuck and mods applied | Go through the single exit |
| F5 | Failure of `completeSolutionFromDiscreteAssignment()` | Returned without undoing infinite-cost mods | Undone |
| F6 | Errors before the solve runs: infinite costs that can't be handled, failure to complete the solution | Infinite costs: `info_.valid` set to `true` over stale info from an earlier solve, run data not invalidated. Failure to complete: returned with whatever the nested call left. (An invalid semi-variable already invalidated everything, through `returnFromOptimizeModel()` with `kSolveError`.) | Info, run data, solution and basis invalidated |
| F7 | `callRunPostsolve()` errors | Called `returnFromOptimizeModel()` outside a solve, which hit its assert in debug | Own clean-up: invalidates info, run data, solution and basis |
| F8 | `Highs::run()` with an error from `optimizeHighs()` or from `model_.assess()` | Returned without `userUnscale()`, leaving the model scaled | Unscaled, then the error is returned (but see B11) |
| F9 | Multiple objectives (`multiobjectiveSolve()`, `returnFromLexicographicOptimization()`) | `col_cost_` replaced without updating `has_infinite_cost_`: an assert in debug, or infinite costs reaching the solvers | Flag recomputed whenever the costs change |
| F10 | Option changes in `HApp.h`, `HighsInterface.cpp`, `HighsSolve.cpp` and the LP pipeline | Not restored if an exception was thrown in between | Restored |

The handling of exceptions in `optimizeModelTryCatch()` was already
fixed in `latest` (3dde8fdc93): `handleCatch` keeps the model status
across `clearSolver()`, undoes the mods and clears the re-entry flag.
This branch keeps that behaviour, with `context.undo_mods` and
`model_mods_applied_` in place of `mods_were_clear`.

### Changes in behaviour

These change results or output in cases that weren't obviously wrong
before.

| # | Change | Effect |
|---|---|---|
| B1 | **Order of mods.** Semi-variables are assessed (or relaxed) together with the infinite-cost handling, before `completeSolutionFromDiscreteAssignment()`. Before, only infinite costs came first. This is what makes F2 possible. | The completion step sees the assessed bounds and types: (a) a semi-continuous variable with lower bound 0 is now continuous, so the user's value for it is no longer fixed; (b) a semi-variable with inconsistent bounds is fixed at 0, so a non-zero user value makes the user's solution infeasible; (c) an upper bound above `kMaxSemiVariableUpper` has been tightened (or the lower bound relaxed with `solve_relaxation`) before the feasibility check. The completed solution, which is the MIP solver's starting point, can differ from before. The model solved by the MIP solver is unchanged. |
| B2 | **Order of checks.** Semi-variable assessment now runs before the empty-model, bound-infeasibility and `assessLp` checks. | A model with both an invalid semi-variable and one of those problems now gets the semi-variable error (`kSolveError`). `processedBoundsFindsInfeasibility()` now sees assessed bounds. I found no case where its result changes. |
| B3 | **Invalid semi-variable error** returns before the solve runs (follows from B2). | It goes through `returnFromOptimizeModelBeforeRunning()`, not `returnFromOptimizeModel()`, so there are no debug checks, no `reportSolvedLpQpStats()` and no `returnFromHighs()`. The status is still `kError` with `kSolveError`. |
| B4 | **Failure to complete a MIP solution** now invalidates the solution (F6). | The user's partial solution was already cleared by `completeSolutionFromDiscreteAssignment()` before its nested call. What's new is that the solution, basis and info left by the nested call are invalidated. The model status left by the nested call is kept. |
| B5 | **QP exceptions** are handled by the single catch in `optimizeModelTryCatch()`, not by a catch around `callSolveQp()`. | `clearSolver()` instead of `returnFromOptimizeModel()`: no debug checks, no stats report, no `returnFromHighs()`. The status is the same (`kMemoryLimit` or `kSolveError`). Exceptions thrown inside the solvers themselves are still caught by `solveCatchCall()` in `HighsSolve.cpp`. |
| B6 | **`callRunPostsolve()` errors** (F7). | No `reportSolvedLpQpStats()` output. Info, solution and basis are always invalidated. Before, that only happened if the model status was an error status. |
| B7 | **Re-entry refusal** logs a user error instead of a developer error and an assert. | Visible in the log. |
| B8 | **`run_data_.valid`** in `Highs::run()` is set after unscaling. | Stays `false` if `userUnscale()` returns an error. Before, it was `true`. |
| B9 | **`assessLp` warning** in `optimizeModelPrepare()` no longer reaches the QP/MIP status. | None in practice: it is asserted never to happen, and the assert in `returnFromOptimizeModel()` rejected it anyway. |
| B10 | **Raised `highs_debug_level`** is undone at the end of the solve. | None in practice: it only applies when the option is below its minimum, which option checks prevent. |
| B11 | **Unscaling after an error in `Highs::run()`** (F8). | `userUnscale()` calls `userScaleSolution(data, true)`. When nothing was solved (a `model_.assess()` failure, or an error from `optimizeHighs()` before `info_` is invalidated, such as a refused re-entry), it multiplies the stale `info_.objective_function_value` by the unscaling factor, although `userScale()` never scaled it. It also reruns the KKT check on the stale solution and logs "After solving the user-scaled model...". On these paths, only the model scaling should be undone, or `userScaleSolution()` called without `update_kkt`. |
| B12 | **MIP debug checks with the restored tolerance** (F1). | With `highs_debug_level` at least cheap, `debugHighsSolution()` in `returnFromOptimizeModel()` now recomputes the MIP solution's infeasibilities with `primal_feasibility_tolerance`, not `mip_feasibility_tolerance`. An optimal MIP with violations between the two (1e-7 and 1e-6 with defaults) now gives "primal infeasibilities but model status is Optimal" and `kError`. The clobbered tolerance used to hide this. |
| B13 | **Order of exception messages in the log.** | In the IPM and PDLP paths of `solveLp()`, the exception message is now printed after the error from `interpretCallStatus()`. In `solveMip()`, it is printed after the profiling clock is stopped. |

## Other points for the review

- **Public API.** `HighsStatus calledOptimizeModel();` was in the
  public section of `Highs` (under "advanced methods: only for internal
  use"). It is now private and takes a `Context&`. Nothing in the
  repository calls it.
- **Installed header.** `Highs.h` now includes `util/HighsUndoStack.h`,
  so the new header has to be installed with it. It is listed in
  `cmake/sources.cmake`. `highs/meson.build` lists no headers, so
  nothing is needed there.
- **Unchanged dead code.** The PDLP clean-up in `lpCleanupOriginal()`
  is still disabled by `consider_pdlp_cleanup = false`, and
  `possibly_use_log_dev_level_2` is still `false`. Both were moved
  without change.

### Known limitations (not introduced here)

These are in `latest` too, but they contradict the aim that the user's
model is returned unchanged.

- **Undo order.** `optimizeModelUndoMods()` calls `restoreInfCost()`
  before `unapplyMods()`, although infinite costs are applied first.
  Example: a semi-variable with bounds [l, u] and a `-inf` cost,
  minimising, with `solve_relaxation`. `handleInfCost()` sets the
  bounds to [u, u], then relaxing saves the lower bound u and sets it
  to 0. Undoing restores [l, u], then sets the lower bound back to u,
  so the user's bounds end up as [u, u]. Swapping the two calls isn't
  enough, since `unapplyMods()` ends by clearing all of `mods_`,
  including the infinite-cost records.
- **Rounded bounds for infinite costs.** For an integer column,
  `handleInfCost()` records the rounded bounds (`ceil` of the lower,
  `floor` of the upper), so `restoreInfCost()` writes the rounded
  values back into the user's model.
- **Corrected bounds are permanent.** Bounds with inconsistencies less
  than the primal feasibility tolerance are changed in the user's
  model by `processedBoundsFindsInfeasibility()`, and not restored.
  Restoring them is left for a separate PR, since the kept simplex
  basis, `presolve()` after `run()` and the MIP sub-LPs would all have
  to handle the restored bounds.

## Where each part went

These are moves. Apart from the items listed under
[Behaviour changes](#behaviour-changes), the only edits are locals
becoming `Context`/`LpState` fields, `return
returnFromOptimizeModel(...)` becoming `return true`, and manual
option save/restore becoming `HighsUndoStack`.

| Old part of `calledOptimizeModel()` | Now |
|---|---|
| Re-entry check, debug level, name hashes, warm start, `exactResizeModel()`, `model_.assess()` assert (added in `latest`, moved here so that it checks the model before any mods) | `optimizeModelPreChecks()` |
| `handleInfCost()`, `assessSemiVariables()`, `relaxSemiVariables()` (in the QP and MIP branches) | `optimizeModelApplyMods()`, undone by `optimizeModelUndoMods()` |
| Call of `completeSolutionFromDiscreteAssignment()` | `optimizeModelCompleteSolution()` |
| Reset info and timer, empty model, infeasible bounds, `assessLp` (the large-value check was removed in `latest`, which checks the model with `model_.assess()` in `run()`) | `optimizeModelPrepare()` |
| QP branch | `optimizeModelQp()` |
| MIP branch | `optimizeModelMip()` |
| Debug log level set-up, basis set-up | `lpSetup()` |
| iCrash (inline) | `lpICrash()`, called from `lpSetup()` |
| `timedSolveLp` lambda | `lpTimedSolve()` |
| Choice of whether to presolve | `lpUsePresolve()` |
| Solve without presolve | `lpSolveUnreduced()` |
| `runPresolve()` and reduction reporting | `lpPresolve()` |
| Terminal presolve statuses | `lpHandlePresolveStatus()`, with `lpResolveUnboundedOrInfeasible()` |
| Assess and `cleanBounds` on the reduced LP | `lpPrepareReducedLp()` |
| Trivial solution when reduced to empty | `lpSetReducedToEmptySolution()` |
| Solves in the presolve switch | `lpSolveReduced()` |
| Postsolve | `lpPostsolve()` |
| Simplex hot start after postsolve, PDLP clean-up | `lpCleanupOriginal()` |
| KKT check, info, final status | `lpFinish()` |
| Report times | `lpReportTimes()` |

Dead code:

- The `use_highs_crossover` branch in iCrash (constant `false`) is
  removed.
- In `solveMip()`, the save and restore of `log_dev_level` around a
  commented-out change is replaced by a `HighsUndoStack` with the
  change still commented out, so it records nothing.
