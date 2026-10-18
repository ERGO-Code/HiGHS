# HighsHessian

A Hessian matrix is communicated via an instance of the HighsHessian class.

- `dim_`: Scalar of type integer - Dimension of the Hessian
- `format_`: Scalar of [HessianFormat](@ref) type - Format of the Hessian
- `start_`: Vector of integer type - Start of each compressed column in the Hessian
- `index_`: Vector of integer type - Indices of the nonzeros in the Hessian
- `value_`: Vector of double type - Values of the nonzeros in the Hessian
