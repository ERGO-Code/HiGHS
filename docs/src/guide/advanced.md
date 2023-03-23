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
of ``Bx=b`` and ``Bx=b`` respectively. Finally, the
methods
`getReducedRow`
and
`getReducedColumn`
yield a specific row or column of ``B^{-1}A``. In all cases,
HiGHS can return the number and indices of the nonzeros in the result.

