# Logic flow of `Highs::optimizeHighs()`

This describes the call sequence from `Highs::optimizeHighs()`, which
`Highs::run()` calls after any file reading and user scaling, down to
the solvers.

## From `Highs::optimizeHighs()` to the problem type

```mermaid
flowchart TD
    optimizeHighs["optimizeHighs()"] -->|multiple objectives| multi["multiobjectiveSolve()"]
    optimizeHighs -->|single objective| optimizeModel["optimizeModel()<br/>initializeMultiThreading(), profiling"]
    multi -->|once per objective| optimizeModel

    optimizeModel --> tryCatch["optimizeModelTryCatch()<br/>creates Context"]
    optimizeLp["optimizeLp()<br/>(LP-only entry point)"] --> tryCatch
    tryCatch --> called["calledOptimizeModel(context)"]
    tryCatch -.->|exception| catch["clearSolver()<br/>optimizeModelUndoMods()"]

    called --> preChecks["optimizeModelPreChecks()"]
    preChecks --> applyMods["optimizeModelApplyMods()<br/>infinite costs, semi-variables"]
    applyMods -->|"MIP with a<br/>(partial) solution"| complete["completeSolutionFromDiscreteAssignment()<br/>may call optimizeModel() again"]
    applyMods --> prepare
    complete --> prepare["optimizeModelPrepare()<br/>empty model, bound infeasibility, large values"]

    prepare -->|QP| qp["optimizeModelQp()<br/>callSolveQp()"]
    prepare -->|MIP| mip["optimizeModelMip()"]
    prepare -->|LP| lp["optimizeModelLp()<br/>(see below)"]
    mip -->|solve_relaxation| relax["optimizeModelLp()<br/>on the LP relaxation"]
    mip -->|otherwise| callMip["callSolveMip()"]

    prepare -->|early exit| returnFrom
    qp --> returnFrom["returnFromOptimizeModel(status, &context)<br/>optimizeModelUndoMods(), returnFromHighs()"]
    callMip --> returnFrom
    relax --> returnFrom
    lp --> returnFrom
```

Early exits from `optimizeModelPreChecks()` and
`optimizeModelApplyMods()` return straight to `optimizeModelTryCatch()`
without calling `returnFromOptimizeModel()`. `optimizeModelApplyMods()`
undoes any mods itself before returning.

Any change to `options_` during the optimization is made through
`context.options_undo`, and anything still recorded there is undone
when the `Context` is destroyed at the end of
`optimizeModelTryCatch()`.

## `optimizeModelLp()`

```mermaid
flowchart TD
    setup["lpSetup(context)<br/>iCrash (lpICrash), basis from solution"] --> usePresolve{"lpUsePresolve(context, incumbent_lp)"}

    usePresolve -->|"no: useful basis,<br/>presolve off,<br/>or no matrix"| unreduced["lpSolveUnreduced(context, incumbent_lp)"]
    unreduced --> finish

    usePresolve -->|yes| presolve["lpPresolve(context, presolve_, incumbent_lp)<br/>runPresolve(presolve_, incumbent_lp)"]
    presolve --> handle{"lpHandlePresolveStatus()"}
    handle -->|"infeasible, timeout,<br/>out of memory, error"| stop["return<br/>(skips lpFinish)"]
    handle -->|unbounded or infeasible| uoi["lpResolveUnboundedOrInfeasible(context, model_.lp_)<br/>primal simplex on the original LP"]
    uoi --> stop
    handle -->|reduced| prepareReduced["lpPrepareReducedLp(context, reduced LP)"]
    handle -->|reduced to empty| empty["lpSetReducedToEmptySolution()"]
    handle -->|not presolved / not reduced| solveReduced

    prepareReduced --> solveReduced["lpSolveReduced(context, status, presolved_lp)<br/>presolved_lp = reduced LP if reduced, else incumbent_lp"]
    empty --> solveReduced
    solveReduced --> postsolve["lpPostsolve(context, presolve_, status, incumbent_lp)<br/>runPostsolve(presolve_, incumbent_lp)<br/>sets lp_state.postsolve_recovered_basis"]
    postsolve --> cleanup["lpCleanupOriginal(context, incumbent_lp)<br/>simplex hot start if a basis was recovered,<br/>PDLP clean-up (currently disabled)"]
    cleanup --> finish["lpFinish(context, incumbent_lp)<br/>KKT check, lpReportTimes(), final status"]
```

Each `lp*` function that can stop the solve returns `true` to do so,
with `context.status` set. The driver then returns straight away, so
`lpFinish()` is skipped.

Every LP solve in this pipeline goes through
`lpTimedSolve(context, lp, description, time)`, which calls
`callSolveLp()` and then `solveLp()` to dispatch to simplex, IPX, PDLP
or HiPO.

State shared between the `lp*` functions is held in `context.lp_state`
(`Highs::LpState`).
