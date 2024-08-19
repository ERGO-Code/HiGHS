# [Enums](@id structures-enums)

The members of the fundamental HiGHS enums are defined below. If `Enum` refers
to a particular enum, and `Member` to a particular member, the members are
available as follows.

 * Python: `highspy.Enum.Member`
 * C++: `Enum::Member`

Members for internal use only are not documented.

## HighsStatus

This is (part of) the return value of most HiGHS methods:

 * `kError`: The method has exposed an error
 * `kOk`: The method has completed successfully
 * `kWarning`: The method has recovered from an unusual event, or has terminated
  due to reaching a time or iteration limit

## MatrixFormat

This defines the format of a [HighsSparseMatrix](@ref):

 * `kColwise`: The matrix is stored column-wise
 * `kRowwise`: The matrix is stored row-wise

## ObjSense

This defines optimization sense of a [HighsLp](@ref):

 * `kMinimize`: The objective is to be minimized
 * `kMaximize`: The objective is to be maximized

## HighsVarType

This defines the feasible values of a variable within a model:

 * `kContinuous`: The variable can take continuous values between its bounds
 * `kInteger`: The variable must take integer values between its bounds
 * `kSemiContinuous`: The variable must be zero or take continuous values between its bounds
 * `kSemiInteger`: The variable must be zero or take integer values between its bounds

## SolutionStatus

This defines the nature of any primal or dual solution information:

 * `kSolutionStatusNone`: There is no solution information
 * `kSolutionStatusInfeasible`: The solution is not feasible
 * `kSolutionStatusFeasible`: The solution is feasible

## BasisValidity

This defines the nature of any basis information:

 * `kBasisValidityInvalid`: There is no basisn information
 * `kBasisValidityValid`: The basis information is valid

## HighsModelStatus

This defines the status of the model after a call to `run`

 * `kNotset`: The model status has not been set
 * `kModelError`: There is an error in the model
 * `kSolveError`: There has been an error when solving the model
 * `kModelEmpty`: The model is empty
 * `kOptimal`: The model has been solved to optimality
 * `kInfeasible`: The model is infeasible
 * `kUnboundedOrInfeasible`: The model is unbounded or infeasible
 * `kUnbounded`: The model is unbounded
 * `kObjectiveBound`: The bound on the model objective value has been reached
 * `kObjectiveTarget`: The target value for the model objective has been reached
 * `kTimeLimit`: The run time limit has been reached
 * `kIterationLimit`: The iteration limit has been reached
 * `kSolutionLimit`: The MIP solver has reached the limit on the number of LPs solved
 * `kInterrupt`: The solver has been interrupted by the user
 * `kMemoryLimit`: The solver has been unable to allocate sufficient memory
 * `kUnknown`: The model status is unknown

## HighsBasisStatus

This defines the status of a variable (or slack variable for a constraint) in a
basis:

 * `kLower`: The variable is nonbasic at its lower bound (or fixed value)
 * `kBasic`: The variable is basic
 * `kUpper`: The variable is at its upper bound
 * `kZero`: A free variable is nonbasic and set to zero
 * `kNonbasic`: The variable is nonbasic

## HighsOptionType

This defines the types of option values that control HiGHS:

 * `kBool`: The option type is boolean
 * `kInt`: The option type is integer
 * `kDouble`: The option type is double
 * `kString`: The option type is string

## HighsInfoType

This defines the types of (scalar) information available after a call to `run`:

 * `kInt64`: The information type is 64-bit integer
 * `kInt`: The information type is integer
 * `kDouble`: The information type is double

