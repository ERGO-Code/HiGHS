# HighsLp

An LP model is communicated via an instance of the HighsLp class

- num\_col\_: Scalar of type integer - Number of columns in the model
- num\_row\_: Scalar of type integer - Number of rows in the model
- col\_cost\_: Vector of type double - Coefficients of the linear term in the objective function
- col\_lower\_: Vector of type double - Lower bounds on the variables
- col\_upper\_: Vector of type double - Upper bounds on the variables
- row\_lower\_: Vector of type double - Lower bounds on the constraints
- row\_upper\_: Vector of type double - Upper bounds on the constraints
- a\_matrix\_: Instance of [HighsSparseMatrix](@ref) class - Constraint matrix
- sense\_: Scalar of type [ObjSense](@ref) - Optimization sense of the model
- offset\_: Scalar of type double - Constant term in the objective function
- model\_name\_: Scalar of type string - Name of the model
- objective\_name\_: Scalar of type string - Name of the objective function
- col\_names\_: Vector of type string - Names of the variables
- row\_names\_: Vector of type string - Names of the constraints
- integrality\_: Vector of type [HighsVarType](@ref) - Type of each variable
