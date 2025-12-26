namespace Highs.Enums;

/// <summary>
/// This defines the status of the model after a call to run
/// </summary>
public enum ModelStatus
{
    /// <summary>
    /// The model status has not been set
    /// </summary>
    Notset = 0,
    LoadError,
    /// <summary>
    /// There is an error in the model
    /// </summary>
    ModelError,
    PresolveError,
    /// <summary>
    /// There has been an error when solving the model
    /// </summary>
    SolveError,
    PostsolveError,
    /// <summary>
    /// The model is empty
    /// </summary>
    ModelEmpty,
    /// <summary>
    /// The model has been solved to optimality
    /// </summary>
    Optimal,
    /// <summary>
    /// The model is infeasible
    /// </summary>
    Infeasible,
    /// <summary>
    /// The model is unbounded or infeasible
    /// </summary>
    UnboundedOrInfeasible,
    /// <summary>
    /// The model is unbounded
    /// </summary>
    Unbounded,
    /// <summary>
    /// The bound on the model objective value has been reached
    /// </summary>
    ObjectiveBound,
    /// <summary>
    /// The target value for the model objective has been reached
    /// </summary>
    ObjectiveTarget,
    /// <summary>
    /// The run time limit has been reached
    /// </summary>
    TimeLimit,
    /// <summary>
    /// The iteration limit has been reached
    /// </summary>
    IterationLimit,
    /// <summary>
    /// The model status is unknown
    /// </summary>
    Unknown,
    /// <summary>
    /// The MIP solver has reached the limit on the number of LPs solved
    /// </summary>
    SolutionLimit,
    /// <summary>
    /// The solver has been interrupted by the user
    /// </summary>
    Interrupt,
    /// <summary>
    /// The solver has been unable to allocate sufficient memory
    /// </summary>
    MemoryLimit
}
