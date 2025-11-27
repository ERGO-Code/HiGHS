# HighsIis

Irreducible infeasibility system (IIS) data are communicated via an instance of the `HighsIis` class.

- `valid_`: The data in the `HighsIis` instance is valid
- `strategy_`: The IIS strategy used
- `col_index_`: The indices of model columns in the IIS
- `row_index_`: The indices of model rows in the IIS
- `col_bound_`: The bounds on each column that define the IIS
- `row_bound_`: The bounds on each row that define the IIS
- `col_status_`: Indicates whether a column in the model is in an IIS, may be in an IIS, or is not in an IIS
- `row_status_`: Indicates whether a row in the model is in an IIS, may be in an IIS, or is not in an IIS
- `info_`: Data on the time and number of simplex iterations required to form the IIS
- `model_`: A [HighsModel](@ref) consisting of the variables, constraints and bounds in the IIS. Currently only its [HighsLp](@ref) instance is relevant


