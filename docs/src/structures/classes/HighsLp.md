# HighsLp

An LP model is communicated via an instance of the `HighsLp` class

- `num_col_`: Scalar of type integer - Number of columns in the model
- `num_row_`: Scalar of type integer - Number of rows in the model
- `col_cost_`: Vector of type double - Coefficients of the linear term in the objective function
- `col_lower_`: Vector of type double - Lower bounds on the variables
- `col_upper_`: Vector of type double - Upper bounds on the variables
- `row_lower_`: Vector of type double - Lower bounds on the constraints
- `row_upper_`: Vector of type double - Upper bounds on the constraints
- `a_matrix_`: Instance of [HighsSparseMatrix](@ref) class - Constraint matrix
- `sense_`: Scalar of type [ObjSense](@ref) - Optimization sense of the model
- `offset_`: Scalar of type double - Constant term in the objective function
- `model_name_`: Scalar of type string - Name of the model
- `objective_name_`: Scalar of type string - Name of the objective function
- `col_names_`: Vector of type string - Names of the variables
- `row_names_`: Vector of type string - Names of the constraints
- `integrality_`: Vector of type [HighsVarType](@ref) - Type of each variable
