The constraint matrix of an LP model is communicated via an instance of the HighsSparseMatrix class 

- format\_: Scalar of MatrixFormat type - Format of the matrix
- num\_col\_ : Scalar of integer type - Number of columns in the matrix
- num\_row\_: Scalar of integer type - Number of rows in the matrix
- start\_: Vector of integer type - Start of each compressed vector in the matrix
- index\_: Vector of integer type - Indices of the nonzeros in the matrix
- value\_: Vector of double type - Values of the nonzeros in the matrix
