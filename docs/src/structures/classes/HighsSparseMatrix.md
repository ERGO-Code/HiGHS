# HighsSparseMatrix

The constraint matrix of an LP model is communicated via an instance of the HighsSparseMatrix class

- `format_`: Scalar of [MatrixFormat](@ref) type - Format of the matrix
- `num_col_ `: Scalar of integer type - Number of columns in the matrix
- `num_row_`: Scalar of integer type - Number of rows in the matrix
- `start_`: Vector of integer type - Start of each compressed vector in the matrix
- `index_`: Vector of integer type - Indices of the nonzeros in the matrix
- `value_`: Vector of double type - Values of the nonzeros in the matrix
