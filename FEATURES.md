## Code changes

Prompted by [#2821](https://github.com/ERGO-Code/HiGHS/issues/2821),
the treatment of Hessian matrix anomalies has been changed. Firstly,
any duplicate entries in the Hessian are now summed.

- When a square Hessian is read from the `QMATRIX` section of an MPS
  file, or passed by a user, any asymmetry results in
  `Highs::readModel` or `Highs::passHessian` returning
  `HighsStatus::kError`. Previously HiGHS would use $$(Q+Q^T)/2$$ as
  the Hessian.

- A triangular Hessian, whether read from the `QUADOBJ` section of an
  MPS file, or passed by a user, was previsouly assumed to be given by
  only lower triangular entries, with any entries in the upper
  triangle being ignored. Now, any entries in the upper triangle of
  the Hessian are accepted, being added to any corresponding entries
  in the lower triangle. If there are entries in the upper triangle,
  their number is logged in a warning message, which also states the
  number of any summations, and `Highs::readModel` or
  `Highs::passHessian` will return `HighsStatus::kWarning`.

Following PR [#2854](https://github.com/ERGO-Code/HiGHS/pull/2854),
HiPO is now capable of solving convex QP problems. Option
solver="qpasm" selects the previous active-set QP solver, while
solver="hipo" or solver="ipm" selects the HiPO solver.

## Build changes

