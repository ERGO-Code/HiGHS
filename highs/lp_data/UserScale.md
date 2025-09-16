# HiGHS - User bound and cost scaling

User bound and cost scaling have been introduced to assess solver
performance if a user's model were better scaled. It is achieved using
the integer `HighsOptions::user_bound_scale` and
`HighsOptions::user_cost_scale` options, where the scaling factor
itself is two to the power of the option value. For simplicity, the
scaling is applied to the incumbent model - as if the model had been
well scaled by the user.

Scaling of costs is simple, since all costs are scaled. Scaling of LPs
is also simple since all finite column and row bounds are scaled, so
the constraint matrix is unchanged. However, nonzero (finite) bounds
on non-continuous variables cannot be scaled, since it would change
the number of feasible values. For a non-continuous variable, its
contributions to the objective and constraint activites are scaled by
scaling the cost coefficient and matrix values.

In summary, the incumbent data affected by user bound and cost scaling
are as follows

- User bound scaling affects
  - Costs of non-continuous variables &#8212; using `userScaleNonContinuousCosts`
  - Matrix coefficients of non-continuous variables &#8212; using `userScaleNonContinuousMatrix`
  - Bounds of continuous variables &#8212; using `userScaleContinuousBounds`
  - All row bounds &#8212; using `userScaleRowBounds`
- User cost scaling affects
  - All column costs &#8212; using `userScaleCosts`

The only exception is the case of changes driven by a
`HighsIndexCollection`, where scaling is performed on individual
entries

## Where user scaling takes place

User scaling typically takes place when a model has been passed to
HiGHS and the values of `user_bound_scale` or `user_cost_scale` are
nonzero, or when the values of `user_bound_scale` or `user_cost_scale`
are changed after a model has been passed to HiGHS. However, for
consistency, when the incumbent model has been scaled that must be
maintained when other operations take place.

Although scaling is more likely to be done to reduce the magnitude of
bounds and costs, extreme scaling up can make increase the magnitude
of bounds and costs beyond the value at which HiGHS considers them to
be infinite. This (arguably) changes the model fundamentally. Hence,
in the typical cases of user scaling, an initial pass is made to
determine whether it generates infinite bounds and costs, returning an
error if this occurs.

### `Highs::passModel`



### `Highs::addColsInterface`

When new columns are added to a model, they have no integrality, so
the scaling of their bounds and costs is simple. Hence,
`user_bound_scale` is applied uniformly to the (local) lower and upper
bounds of the new columns, and `user_cost_scale` is applied uniformly
to the (local) costs of the new columns.

### `Highs::addRowsInterface`

When new rows are added to a model, any column integrality must be
respected. Hence `user_bound_scale` must be applied to the columns of
non-continuous variables in the (local) matrix, and must be applied
uniformly to the lower and upper bounds of the new rows

### `Highs::changeIntegralityInterface`

When integrality changes, the simplest thing to do is first remove any
scaling specific to non-continuous variables &#8212; by using the negation
of `user_bound_scale` to scale the costs and matrix coefficient for
non-continuous variables. Then, after updating the integrality, use
`user_bound_scale` to scale the costs and matrix coefficient for
non-continuous variables.

### `Highs::changeCostsInterface`

When costs change, `user_bound_scale` must be applied to the changed costs of non-continuous variables, and `user_cost_scale` must be applied uniformly.

### `Highs::changeColBoundsInterface`

When column bounds change, `user_bound_scale` must be applied to the changed bounds of continuous variables

### `Highs::changeRowBoundsInterface`

When row bounds change, `user_bound_scale` must be applied uniformly to the changed bounds

### `Highs::changeCoefficientInterface

The coefficient is scaled using `user_bound_scale` if is corresponds to a non-continuous variable

### `Highs::optionChangeAction`

When `user_bound_scale` or `user_cost_scale` change, the difference in
value must be applied as for a whole model.

### `Highs::passModel`


