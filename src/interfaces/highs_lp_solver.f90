module highs_lp_solver
  interface
    function Highs_call (n, m, nz, cc, cl, cu, rl, ru, as, ai, av, cv, cd, rv, rd, cbs, rbs) result(s) bind (c, name='Highs_call')
      use iso_c_binding
      integer ( c_int ), VALUE :: n
      integer ( c_int ), VALUE :: m
      integer ( c_int ), VALUE :: nz
      real ( c_double ) :: cc(*)
      real ( c_double ) :: cl(*)
      real ( c_double ) :: cu(*)
      real ( c_double ) :: rl(*)
      real ( c_double ) :: ru(*)
      integer ( c_int ) :: as(*)
      integer ( c_int ) :: ai(*)
      real ( c_double ) :: av(*)
      real ( c_double ) :: cv(*)
      real ( c_double ) :: cd(*)
      real ( c_double ) :: rv(*)
      real ( c_double ) :: rd(*)
      integer ( c_int ) :: cbs(*)
      integer ( c_int ) :: rbs(*)
      integer ( c_int ) :: s
    end function Highs_call

    function Highs_create () result ( h ) bind( c, name='Highs_create' )
      use iso_c_binding
      type(c_ptr) :: h
    end function Highs_create
      
    subroutine Highs_destroy ( h ) bind( c, name='Highs_destroy' )
      use iso_c_binding
      type(c_ptr), VALUE :: h 
    end subroutine Highs_destroy

    function Highs_run ( h ) result ( s ) bind( c, name='Highs_run' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ) :: s
    end function Highs_run

    function Highs_readFromFile ( h, f ) result ( s ) bind ( c, name='Highs_readFromFile' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: f(*)
      integer ( c_int ) :: s
    end function Highs_readFromFile

    function Highs_writeToFile ( h, f ) result ( s ) bind ( c, name='Highs_writeToFile' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: f(*)
      integer ( c_int ) :: s
    end function Highs_writeToFile

    function Highs_loadModel ( h, n, m, nz, cc, cl, cu, rl, ru, as, ai, av) result ( s ) bind ( c, name='Highs_loadModel' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: n
      integer ( c_int ), VALUE :: m
      integer ( c_int ), VALUE :: nz
      real ( c_double ) :: cc(*)
      real ( c_double ) :: cl(*)
      real ( c_double ) :: cu(*)
      real ( c_double ) :: rl(*)
      real ( c_double ) :: ru(*)
      integer ( c_int ) :: as(*)
      integer ( c_int ) :: ai(*)
      real ( c_double ) :: av(*)
      integer ( c_int ) :: s
    end function Highs_loadModel

    function Highs_setOptionValue ( h, o, v ) result( s ) bind ( c, name='Highs_setOptionValue' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      character( c_char ) :: o(*)
      character( c_char ) :: v(*)
      integer ( c_int ) :: s
    end function Highs_setOptionValue

    function Highs_getSolution (h, cv, cd, rv, rd) result ( s ) bind ( c, name='Highs_getSolution' )
      use iso_c_binding
      type(c_ptr), VALUE :: h
      real ( c_double ) :: cv(*)
      real ( c_double ) :: cd(*)
      real ( c_double ) :: rv(*)
      real ( c_double ) :: rd(*)
      integer ( c_int ) :: s
    end function Highs_getSolution

    function Highs_getBasis (h, cbs, rbs) result(s) bind (c, name='Highs_getBasis')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ) :: cbs(*)
      integer ( c_int ) :: rbs(*)
      integer ( c_int ) :: s
    end function Highs_getBasis

    function Highs_getObjectiveValue (h) result(ov) bind(c, name='Highs_getObjectiveValue')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      real ( c_double ) :: ov
    end function Highs_getObjectiveValue

    function Highs_getIterationCount (h) result(ic) bind(c, name='Highs_getIterationCount')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ) :: ic
    end function Highs_getIterationCount

    function Highs_addRow (h, lo, up, nz, idx, val) result(s) bind(c, name='Highs_addRow')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      real ( c_double ), VALUE :: lo
      real ( c_double ), VALUE :: up
      integer ( c_int ), VALUE :: nz
      integer ( c_int ) :: idx(*)
      real ( c_double ) :: val(*)
      integer ( c_int ) :: s
    end function Highs_addRow

    function Highs_addRows (h, nnr, lo, up, nnz, st, idx, val) result(s) bind(c, name='Highs_addRows')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: nnr
      real ( c_double ) :: lo(*)
      real ( c_double ) :: up(*)
      integer ( c_int ), VALUE :: nnz
      integer ( c_int ) :: st(*)
      integer ( c_int ) :: idx(*)
      real ( c_double ) :: val(*)
      integer ( c_int ) :: s
    end function Highs_addRows

    function Highs_addCol (h, cc, cl, cu, nnz, idx, val) result(s) bind(c, name='Highs_addCol')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      real ( c_double ) :: cc
      real ( c_double ) :: cl
      real ( c_double ) :: cu
      integer ( c_int ), VALUE :: nnz
      integer ( c_int ) :: idx(*)
      real ( c_double ) :: val(*)
      integer ( c_int ) :: s
    end function Highs_addCol

    function Highs_addCols (h, nnc, cc, cl, cu, nnz, st, idx, val) result(s) bind(c, name='Highs_addCols')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: nnc
      real ( c_double ) :: cc(*)
      real ( c_double ) :: cl(*)
      real ( c_double ) :: cu(*)
      integer ( c_int ), VALUE :: nnz
      integer ( c_int ) :: st(*)
      integer ( c_int ) :: idx(*)
      real ( c_double ) :: val(*)
      integer ( c_int ) :: s
    end function Highs_addCols

    function Highs_changeObjectiveSense (h, sns) result(s) bind(c, name='Highs_changeObjectiveSense')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: sns
      integer ( c_int ) :: s
    end function Highs_changeObjectiveSense
    
    function Highs_changeColCost (h, c, co) result(s) bind(c, name='Highs_changeColCost')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: c
      real ( c_double ), VALUE :: co
      integer ( c_int ) :: s
    end function Highs_changeColCost

    function Highs_changeColsCostBySet (h, nse, set, cost) result(s) bind(c, name='Highs_changeColsCostBySet')
      use iso_c_binding
      type(c_ptr), VALUE :: h
      integer ( c_int ), VALUE :: nse
      integer ( c_int ) :: set(*)
      real ( c_double ) :: cost(*)
      integer ( c_int ) :: s
    end function Highs_changeColsCostBySet


! int Highs_changeColsCostByMask(
!     void *highs,        //!< HiGHS object reference
!     const int *mask,    //!< Full length array with 1 => change; 0 => not
!     const double *cost  //!< Full length array of new costs
! );

! /**
!  * @brief Change the bounds of a column
!  */
! int Highs_changeColBounds(
!     void *highs,         //!< HiGHS object reference
!     const int col,       //!< The index of the column whose bounds are to change
!     const double lower,  //!< The new lower bound
!     const double upper   //!< The new upper bound
! );

! /**
!  * @brief Change the bounds of multiple columns given by an interval
!  */
! int Highs_changeColsBoundsByRange(
!     void *highs,         //!< HiGHS object reference
!     const int from_col,  //!< The index of the first column whose bounds change
!     const int to_col,    //!< One more than the index of the last column whose
!                          //!< bounds change
!     const double
!         *lower,  //!< Array of size to_col-from_col with new lower bounds
!     const double
!         *upper  //!< Array of size to_col-from_col with new upper bounds
! );

! /**
!  * @brief Change the bounds of multiple columns given by a set of indices
!  */
! int Highs_changeColsBoundsBySet(
!     void *highs,                //!< HiGHS object reference
!     const int num_set_entries,  //!< The number of indides in the set
!     const int *set,  //!< Array of size num_set_entries with indices of
!                      //!< columns whose bounds change
!     const double
!         *lower,  //!< Array of size num_set_entries with new lower bounds
!     const double
!         *upper  //!< Array of size num_set_entries with new upper bounds
! );

! /**
!  * @brief Change the cost of multiple columns given by a mask
!  */
! int Highs_changeColsBoundsByMask(
!     void *highs,          //!< HiGHS object reference
!     const int *mask,      //!< Full length array with 1 => change; 0 => not
!     const double *lower,  //!< Full length array of new lower bounds
!     const double *upper   //!< Full length array of new upper bounds
! );

! /**
!  * @brief Change the bounds of a row
!  */
! int Highs_changeRowBounds(
!     void *highs,         //!< HiGHS object reference
!     const int row,       //!< The index of the row whose bounds are to change
!     const double lower,  //!< The new lower bound
!     const double upper   //!< The new upper bound
! );

! /**
!  * @brief Change the bounds of multiple rows given by a set of indices
!  */
! int Highs_changeRowsBoundsBySet(
!     void *highs,                //!< HiGHS object reference
!     const int num_set_entries,  //!< The number of indides in the set
!     const int *set,  //!< Array of size num_set_entries with indices of rows
!                      //!< whose bounds change
!     const double
!         *lower,  //!< Array of size num_set_entries with new lower bounds
!     const double
!         *upper  //!< Array of size num_set_entries with new upper bounds
! );

! /**
!  * @brief Change the cost of multiple rows given by a mask
!  */
! int Highs_changeRowsBoundsByMask(
!     void *highs,          //!< HiGHS object reference
!     const int *mask,      //!< Full length array with 1 => change; 0 => not
!     const double *lower,  //!< Full length array of new lower bounds
!     const double *upper   //!< Full length array of new upper bounds
! );

! /**
!  * @brief Get multiple columns from the model given by an interval
!  */
! int Highs_getColsByRange(
!     void *highs,          //!< HiGHS object reference
!     const int from_col,   //!< The index of the first column to
!                           //!< get from the model
!     const int to_col,     //!< One more than the last column to get
!                           //!< from the model
!     int* num_col,          //!< Number of columns got from the model
!     double *costs,        //!< Array of size num_col with costs
!     double *lower,        //!< Array of size num_col with lower bounds
!     double *upper,        //!< Array of size num_col with upper bounds
!     int* num_nz,           //!< Number of nonzeros got from the model
!     int *matrix_start,    //!< Array of size num_col with start
!                           //!< indices of the columns
!     int *matrix_index,    //!< Array of size num_nz with row
!                           //!< indices for the columns
!     double *matrix_value  //!< Array of size num_nz with row
!                           //!< values for the columns
! );

! /**
!  * @brief Get multiple columns from the model given by a set
!  */
! int Highs_getColsBySet(
!     void *highs,                //!< HiGHS object reference
!     const int num_set_entries,  //!< The number of indides in the set
!     const int *set,             //!< Array of size num_set_entries with indices
!                                 //!< of columns to get
!     int* num_col,                //!< Number of columns got from the model
!     double *costs,              //!< Array of size num_col with costs
!     double *lower,              //!< Array of size num_col with lower bounds
!     double *upper,              //!< Array of size num_col with upper bounds
!     int* num_nz,                 //!< Number of nonzeros got from the model
!     int *matrix_start,          //!< Array of size num_col with start indices
!                                 //!< of the columns
!     int *matrix_index,          //!< Array of size num_nz with row indices
!                                 //!< for the columns
!     double *matrix_value        //!< Array of size num_nz with row values
!                                 //!< for the columns
! );

! /**
!  * @brief Get multiple columns from the model given by a mask
!  */
! int Highs_getColsByMask(
!     void *highs,          //!< HiGHS object reference
!     const int *mask,      //!< Full length array with 1 => get; 0 => not
!     int* num_col,          //!< Number of columns got from the model
!     double *costs,        //!< Array of size num_col with costs
!     double *lower,        //!< Array of size num_col with lower bounds
!     double *upper,        //!< Array of size num_col with upper bounds
!     int* num_nz,           //!< Number of nonzeros got from the model
!     int *matrix_start,    //!<  Array of size num_col with start
!                           //!<  indices of the columns
!     int *matrix_index,    //!<  Array of size num_nz with row indices
!                           //!<  for the columns
!     double *matrix_value  //!<  Array of size num_nz with row values
!                           //!<  for the columns
! );

! /**
!  * @brief Get multiple rows from the model given by an interval
!  */
! int Highs_getRowsByRange(
!     void *highs,          //!< HiGHS object reference
!     const int from_row,   //!< The index of the first row to get from the model
!     const int to_row,     //!< One more than the last row get from the model
!     int* num_row,          //!< Number of rows got from the model
!     double *lower,        //!< Array of size num_row with lower bounds
!     double *upper,        //!< Array of size num_row with upper bounds
!     int* num_nz,           //!< Number of nonzeros got from the model
!     int *matrix_start,    //!< Array of size num_row with start indices of the
!                           //!< rows
!     int *matrix_index,    //!< Array of size num_nz with column indices for the
!                           //!< rows
!     double *matrix_value  //!< Array of size num_nz with column values for the
!                           //!< rows
! );

! /**
!  * @brief Get multiple rows from the model given by a set
!  */
! int Highs_getRowsBySet(
!     void *highs,                //!< HiGHS object reference
!     const int num_set_entries,  //!< The number of indides in the set
!     const int *set,             //!< Array of size num_set_entries with indices
!                                 //!< of rows to get
!     int* num_row,                //!< Number of rows got from the model
!     double *lower,              //!< Array of size num_row with lower bounds
!     double *upper,              //!< Array of size num_row with upper bounds
!     int* num_nz,                 //!< Number of nonzeros got from the model
!     int *matrix_start,          //!< Array of size num_row with start indices
!                                 //!< of the rows
!     int *matrix_index,          //!< Array of size num_nz with column indices
!                                 //!< for the rows
!     double *matrix_value        //!< Array of size num_nz with column
!                                 //!< values for the rows
! );

! /**
!  * @brief Get multiple rows from the model given by a mask
!  */
! int Highs_getRowsByMask(
!     void *highs,          //!< HiGHS object reference
!     const int *mask,      //!< Full length array with 1 => get; 0 => not
!     int* num_row,          //!< Number of rows got from the model
!     double *lower,        //!< Array of size num_row with lower bounds
!     double *upper,        //!< Array of size num_row with upper bounds
!     int* num_nz,           //!< Number of nonzeros got from the model
!     int *matrix_start,    //!< Array of size num_row with start indices
!                           //!< of the rows
!     int *matrix_index,    //!< Array of size num_nz with column indices
!                           //!< for the rows
!     double *matrix_value  //!< Array of size num_nz with column
!                           //!< values for the rows
! );

! /**
!  * @brief Delete multiple columns from the model given by an interval
!  */
! int Highs_deleteColsByRange(
!     void *highs,         //!< HiGHS object reference
!     const int from_col,  //!< The index of the first column
!                          //!< to delete from the model
!     const int to_col     //!< One more than the last column to
!                          //!< delete from the model
! );

! /**
!  * @brief Delete multiple columns from the model given by a set
!  */
! int Highs_deleteColsBySet(
!     void *highs,                //!< HiGHS object reference
!     const int num_set_entries,  //!< The number of indides in the set
!     const int *set  //!< Array of size num_set_entries with indices of columns
!                     //!< to delete
! );

! /**
!  * @brief Delete multiple columns from the model given by a mask
!  */
! int Highs_deleteColsByMask(
!     void *highs,  //!< HiGHS object reference
!     int *mask     //!< Full length array with 1 => delete; 0 => not
! );

! /**
!  * @brief Delete multiple rows from the model given by an interval
!  */
! int Highs_deleteRowsByRange(
!     void *highs,  //!< HiGHS object reference
!     const int
!         from_row,     //!< The index of the first row to delete from the model
!     const int to_row  //!< One more than the last row delete from the model
! );

! /**
!  * @brief Delete multiple rows from the model given by a set
!  */
! int Highs_deleteRowsBySet(
!     void *highs,                //!< HiGHS object reference
!     const int num_set_entries,  //!< The number of indides in the set
!     const int *set  //!< Array of size num_set_entries with indices of columns
!                     //!< to delete
! );

! /**
!  * @brief Delete multiple rows from the model given by a mask
!  */
! int Highs_deleteRowsByMask(
!     void *highs,  //!< HiGHS object reference
!     int *mask     //!< Full length array with 1 => delete; 0 => not
! );

! /**
!  * @brief Returns the number of columns of the current model
!  */
! int Highs_getNumCols(
!     void* highs //!< HiGHS object reference
! );

! /**
!  * @brief Returns the number of rows of the current model
!  */
! int Highs_getNumRows(
!     void* highs //!< HiGHS object reference
! );

! /**
!  * @brief Returns the number of nonzeroes of the current model
!  */
! int Highs_getNumNz(
!     void* highs //!< HiGHS object reference
! );
  
  
  
  
    end interface

end module highs_lp_solver
