The solution of a model is communicated via an instance of the HighsSolution class

- value\_valid: Scalar of type bool - Indicates whether the column and row values are valid
- dual\_valid: Scalar of type bool - Indicates whether the column and row [duals](https://ergo-code.github.io/HiGHS/terminology.html#Dual-values) are valid
- col\_value: Vector of type double - Values of the columns (variables)
- col\_dual: Vector of type double - Duals of the columns (variables)
- row\_value: Vector of type double - Values of the rows (constraints)
- row\_dual: Vector of type double - Duals of the rows (constraints)
