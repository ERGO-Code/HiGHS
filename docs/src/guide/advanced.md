# [Advanced features](@id guide-advanced)


## Simplex tableau data

HiGHS has a suite of methods for operations with the invertible
representation of the current basis matrix ``B``. To use
these requires knowledge of the corresponding (ordered) basic
variables. This is obtained using the
method
`getBasicVariables`
, with non-negative values being
columns and negative values corresponding to row indices plus one [so
-1 indicates row 0]. Methods
`getBasisInverseRow`
and
`getBasisInverseCol`
yield a specific row or column
of ``B^{-1}``. Methods
`getBasisSolve`
and
`getBasisTransposeSolve`
yield the solution
of ``Bx=b`` and ``B^{T}x=b`` respectively. Finally, the
methods
`getReducedRow`
and
`getReducedColumn`
yield a specific row or column of ``B^{-1}A``. In all cases,
HiGHS can return the number and indices of the nonzeros in the result.

## [Irreducible infeasibility system (IIS) detection](@id highs-iis)

An Irreducible infeasibility system (IIS) consists of a set of
variables and a set of constraints in a model, together with
variable/constraint bound information, that cannot be satisfied (so is
infeasible). It is irreducible in that if any constraint or variable
bound is removed, then the system can be satisfied (so is feasible).

HiGHS has an IIS facility that is under development. Currently it can
only be used for LPs. The full IIS calculation is expensive, since it
requires the solution of multiple LPs. Although there is a prototype
implementation, it is not as robust or efficient as it will
be. Otherwise, there is a simple, cheap test that looks for
infeasibility due to incompatible variable or constraint bounds, or
constraint bounds that cannot be satisfied given the range of values
on the constraint activity implied by bounds on variables.

The choice of IIS strategy is defined by the [iis_strategy](@ref option-iis-strategy) option. This a bit map

- 0 => "light strategy", which is always performed when Highs::getIis is called
- 1 => From dual ray, which is currently unavailable
- 2 => From the whole LP (solving an elasticity LP repeatedly (fixing positive elastic variables at zero) until no more elastic variables are positive, and using the fixed elastic variables to determine a set of infeasible rows, for which there is a corresponding set of columns with nonzeros in those rows that form an infeasibility set (IS)
- 4 => Attempt to reduce the IS to an ISS
- 8 => Prioritize low numbers of columns (rather than low numbers of rows) when reducing the IS

Hence, by just setting the 2-bit, an IS is formed reliably, and at no great expense (for an LP)

### IIS-related methods in the `Highs` class

- `const HighsLp& getIisLp()`: Return a const reference to the internal IIS LP instance
- `HighsStatus getIis(HighsIis& iis)`: Try to find an IIS for the incumbent model. Gets the internal [HighsIis](@ref) instance, returning `HighsStatus::kError` if the calculation failed. Note that if the incumbent model is found to be feasible, this is a "success", and `HighsStatus::kOk` is returned.
- ` HighsStatus writeIisModel(const std::string& filename = "")`: Write out the internal IIS LP instance to a file.



