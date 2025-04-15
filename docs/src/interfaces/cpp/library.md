The HiGHS library is defined in the [`highs/Highs.h`](https://github.com/ERGO-Code/HiGHS/blob/master/src/Highs.h) header file. It contains the definition of the methods and members of the class. 

## Define model

Models in HiGHS are defined as an instance of the `HighsModel` class. This consists of one instance of the `HighsLp` class, and one instance of the `HighsHessian` class. Communication of models to and from HiGHS is possible via instances of the `HighsLp` or `HighsModel` class. In the C and other interfaces, communication of models is via scalar values and addresses of arrays.

In C++, the neatest way of passing a model to HiGHS is to create an instance of the `HighsModel` class, populate its data, and call
``` cpp
Highs::passModel(const HighsModel& model)
```

or create and populate an instance of the `HighsLp` class, and call
``` cpp
Highs::passModel(const HighsLp& lp)
```

For reading models from a file, use
``` cpp
Highs::readModel(const std::string& filename)
```

Below is an example of building a `HighsModel`
``` cpp
  // Create and populate a HighsModel instance for the LP
  
  // Min    f  =  x_0 +  x_1 + 3
  // s.t.                x_1 <= 7
  //        5 <=  x_0 + 2x_1 <= 15
  //        6 <= 3x_0 + 2x_1
  // 0 <= x_0 <= 4; 1 <= x_1
  
  // Although the first constraint could be expressed as an upper
  // bound on x_1, it serves to illustrate a non-trivial packed
  // column-wise matrix.
  
  HighsModel model;
  model.lp_.num_col_ = 2;
  model.lp_.num_row_ = 3;
  model.lp_.sense_ = ObjSense::kMinimize;
  model.lp_.offset_ = 3;
  model.lp_.col_cost_ = {1.0, 1.0};
  model.lp_.col_lower_ = {0.0, 1.0};
  model.lp_.col_upper_ = {4.0, 1.0e30};
  model.lp_.row_lower_ = {-1.0e30, 5.0, 6.0};
  model.lp_.row_upper_ = {7.0, 15.0, 1.0e30};
  
  // Here the orientation of the matrix is column-wise
  model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
  // a_start_ has num_col_+1 entries, and the last entry is the number
  // of nonzeros in A, allowing the number of nonzeros in the last
  // column to be defined
  model.lp_.a_matrix_.start_ = {0, 2, 5};
  model.lp_.a_matrix_.index_ = {1, 2, 0, 1, 2};
  model.lp_.a_matrix_.value_ = {1.0, 3.0, 1.0, 2.0, 2.0};
```

## Solve model

``` cpp
  // Create a Highs instance
  Highs highs;
  HighsStatus return_status;
  
  // Pass the model to HiGHS
  return_status = highs.passModel(model);
  assert(return_status==HighsStatus::kOk);
  
  // Get a const reference to the LP data in HiGHS
  const HighsLp& lp = highs.getLp();
  
  // Solve the model
  return_status = highs.run();
  assert(return_status==HighsStatus::kOk);
  
  // Get the model status
  const HighsModelStatus& model_status = highs.getModelStatus();
  assert(model_status==HighsModelStatus::kOptimal);
```

Solution information:

``` cpp
  const HighsInfo& info = highs.getInfo();
  cout << "Simplex iteration count: " << info.simplex_iteration_count << endl;
  cout << "Objective function value: " << info.objective_function_value << endl;
  cout << "Primal  solution status: " << highs.solutionStatusToString(info.primal_solution_status) << endl;
  cout << "Dual    solution status: " << highs.solutionStatusToString(info.dual_solution_status) << endl;
  cout << "Basis: " << highs.basisValidityToString(info.basis_validity) << endl;
```

## Integrality variables

To indicate that variables must take integer values use the `HighsLp::integrality` vector.
``` cpp
  model.lp_.integrality_.resize(lp.num_col_);
  for (int col=0; col < lp.num_col_; col++)
    model.lp_.integrality_[col] = HighsVarType::kInteger;

  highs.passModel(model);
```
