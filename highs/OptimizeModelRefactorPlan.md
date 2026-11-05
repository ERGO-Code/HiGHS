# Plan for re-applying the `optimizeModel` refactoring

Steps in order. Steps 1–6 should not change results, so a test run
that matches the previous one is enough after each. Steps 7 and 8
change behaviour.

1. **Options undoer.**
   - Add `util/HighsUndoStack.h`, and list it in `cmake/sources.cmake`
     and `highs/meson.build`.
   - Replace the manual save and restore of options with `set`,
     `mark` and `resetTo` in `Highs.cpp`, `HighsInterface.cpp`,
     `HighsSolve.cpp` and `simplex/HApp.h`.
   - Each block takes its own local mark. A lambda that writes to a
     shared outer mark leaves options unrestored.

2. **`Context` struct.**
   - Create it in `optimizeModelTryCatch()` and pass it to
     `calledOptimizeModel()`.
   - It holds the status, the timings, the debug level, `undo_mods`
     and `options_undo`.
   - This only moves locals into the struct.

3. **Extract the top-level phases, without reordering anything.**
   - `optimizeModelPreChecks`, `optimizeModelPrepare`,
     `optimizeModelQp`, `optimizeModelMip`, `optimizeModelLp`,
     `lpICrash`, `lpReportTimes`.
   - `optimizeModelCompleteSolution`, holding the call to
     `completeSolutionFromDiscreteAssignment()` for a MIP with a user
     solution, so that `calledOptimizeModel()` only calls phases. It
     must stay before the point where the solve is marked as running,
     since it may call `optimizeModel()`.
   - Convention: a `bool` return of `true` means stop, with
     `context.status` set.
   - Use the final names (`lp*`, `Context`, `LpState`) from the
     start, to avoid a rename pass later.

4. **Make presolve and postsolve parameterised.**
   - Add overloads `runPresolve(PresolveComponent&, HighsLp&, ...)`
     and `runPostsolve(PresolveComponent&, const HighsLp&)`.
   - Keep the old signatures as wrappers that pass `presolve_` and
     `model_.lp_`, so `presolve()` and `postsolve()` are unchanged.
   - Set `presolve_log_` only when the component is `presolve_`.

5. **Split up the LP pipeline as pure moves.**
   - Add `LpState` and turn the solve lambda into `lpTimedSolve()`.
   - Extract `lpSetup`, `lpUsePresolve`, `lpSolveUnreduced` and
     `lpFinish`.
   - Extract `lpPresolve`, which takes `HighsLp*&` and moves it to
     the reduced LP, together with `lpHandlePresolveStatus`,
     `lpResolveUnboundedOrInfeasible` (always on the original LP),
     `lpPrepareReducedLp` and `lpSetReducedToEmptySolution`.
   - Extract `lpSolveReduced`, `lpPostsolve` and `lpCleanupOriginal`.
     This splits postsolve from the simplex hot start, so add
     `postsolve_recovered_basis` and `presolved_lp_info` to `LpState`.
     It is the only step that splits existing code apart, so review
     it most carefully.

6. **Single exit.**
   - `calledOptimizeModel()` calls `returnFromOptimizeModel(context)`
     from one place, and that function reads the status from the
     context.
   - `callRunPostsolve()` uses its own error clean-up (invalidate
     info, run data, solution and basis, then return `kError`)
     instead of `returnFromOptimizeModel()`. This also fixes the
     `assert(!called_return_from_optimize_model)` that fails on those
     paths in debug builds.
   - Replace `called_return_from_optimize_model` with
     `optimize_model_running_`, used only to stop `optimizeModel()`
     being called from within a solve. Set it where the old flag was
     set false, clear it in `returnFromOptimizeModel()` and
     `handleCatch`, and check it at the start of
     `optimizeModelPreChecks()`, before anything is changed, with a
     user error and no assert. Remove the assert at the start of
     `returnFromOptimizeModel()` and the check in `returnFromHighs()`,
     which the single exit makes redundant.

7. **Reorder the mods.** This is the main change in behaviour.
   - `optimizeModelApplyMods()` and `optimizeModelUndoMods()` handle
     infinite costs and semi-variables before
     `optimizeModelCompleteSolution()`, which undoes them itself if
     completing the solution fails.
   - `model_mods_applied_` stops the nested `optimizeModel()` call
     from applying or undoing them.
   - The undo runs in `returnFromOptimizeModel()`, on the early exits
     and in `handleCatch`.
   - `multiobjectiveSolve()` recomputes `has_infinite_cost_`
     whenever it changes the costs.
   - This fixes the nested call restoring the outer call's infinite
     costs. It also changes what the completion step sees: a
     semi-continuous variable with a zero lower bound is no longer
     fixed at the user's value. Check this with the team before
     applying it, and consider a separate PR.

8. **Exceptions.**
   - Remove the QP-only `try/catch` around `callSolveQp()`.
   - In `handleCatch`, keep the model status across `clearSolver()`,
     so that LP, MIP and QP all report `kMemoryLimit` or
     `kSolveError` after an exception.

9. **Docs.**
   - Update `HighsRun.md`, and add the flow diagram and the
     refactoring note, adjusted to the steps that are kept.
