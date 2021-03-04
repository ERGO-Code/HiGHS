program fortrantest
  use, intrinsic :: iso_c_binding
  use highs_lp_solver
  implicit none

  ! This illustrates the use of Highs_call, the simple F90 interface to
  ! HiGHS. It's designed to solve the general LP problem
  !
  ! Min c^Tx subject to L <= Ax <= U; l <= x <= u
  !
  ! where A is a matrix with m rows and n columns
  !
  ! The scalar n is numcol
  ! The scalar m is numrow
  !
  ! The vector c is colcost
  ! The vector l is collower
  ! The vector u is colupper
  ! The vector L is rowlower
  ! The vector U is rowupper
  !
  ! The matrix A is represented in packed column-wise form: only its
  ! nonzeros are stored
  !
  ! * The number of nonzeros in A is numnz
  !
  ! * The row indices of the nonnzeros in A are stored column-by-column
  ! in aindex
  !
  ! * The values of the nonnzeros in A are stored column-by-column in
  ! avalue
  !
  ! * The position in aindex/avalue of the index/value of the first
  ! nonzero in each column is stored in astart
  !
  ! Note that astart[0] must be zero
  !
  ! After a successful call to Highs_call, the primal and dual
  ! solution, and the simplex basis are returned as follows
  !
  ! The vector x is colvalue
  ! The vector Ax is rowvalue
  ! The vector of dual values for the variables x is coldual
  ! The vector of dual values for the variables Ax is rowdual
  ! The basic/nonbasic status of the variables x is colbasisstatus
  ! The basic/nonbasic status of the variables Ax is rowbasisstatus
  !
  ! The status of the solution obtained is modelstatus
  !
  ! To solve maximization problems, the values in c must be negated
  !
  ! The use of Highs_call is illustrated for the LP
  !
  ! Min    f  = 2x_0 + 3x_1
  ! s.t.                x_1 <= 6
  !       10 <=  x_0 + 2x_1 <= 14
  !        8 <= 2x_0 +  x_1
  ! 0 <= x_0 <= 3; 1 <= x_1
  
  integer ( c_int ), parameter :: numcol = 2
  integer ( c_int ), parameter :: numrow = 3
  integer ( c_int ), parameter :: numnz = 5
  
  real ( c_double ) colcost(numcol)
  real ( c_double ) collower(numcol)
  real ( c_double ) colupper(numcol)
  real ( c_double ) rowlower(numrow)
  real ( c_double ) rowupper(numrow)
  integer ( c_int ) astart(numcol)
  integer ( c_int ) aindex(numnz)
  real ( c_double ) avalue(numnz)
  
  real ( c_double ) colvalue(numcol)
  real ( c_double ) coldual(numcol)
  real ( c_double ) rowvalue(numrow)
  real ( c_double ) rowdual(numrow)
  integer ( c_int ) colbasisstatus(numcol)
  integer ( c_int ) rowbasisstatus(numrow)
  integer modelstatus
  integer runstatus
  
  ! For the full API test
  type ( c_ptr ) :: highs
  
  ! For the row-wise matrix
  integer ( c_int ) arstart(numrow)
  integer ( c_int ) arindex(numnz)
  real ( c_double ) arvalue(numnz)
  
  real, parameter :: inf = 1e30
  integer, parameter :: runstatus_ok = 0
  integer, parameter :: scaled_model = 0
  integer col, row
  integer simplex_iteration_count, primal_status, dual_status
  double precision objective_function_value
  integer option_type

  integer sense
  integer simplex_scale_strategy
  integer, parameter :: default_simplex_scale_strategy = 2
  integer, parameter :: new_simplex_scale_strategy = 3

  double precision, pointer :: double_null(:)
  integer, pointer :: integer_null(:)
  
  colcost(1) = 2
  colcost(2) = 3
  
  collower(1) = 0
  collower(2) = 1
  
  colupper(1) = 3
  colupper(2) = inf
  
  rowlower(1) = -inf
  rowlower(2) = 10.0
  rowlower(3) = 8.0
  
  rowupper(1) = 6.0
  rowupper(2) = 14.0
  rowupper(3) = inf
  
  astart(1) = 0
  astart(2) = 2
  
  aindex(1) = 1
  aindex(2) = 2
  aindex(3) = 0
  aindex(4) = 1
  aindex(5) = 2
  
  avalue(1) = 1
  avalue(2) = 2
  avalue(3) = 1
  avalue(4) = 2
  avalue(5) = 1

  runstatus = Highs_call( numcol, numrow, numnz,&
       colcost, collower, colupper, rowlower, rowupper,&
       astart, aindex, avalue,&
       colvalue, coldual, rowvalue, rowdual,&
       colbasisstatus, rowbasisstatus, modelstatus)

  if (runstatus .ne. runstatus_ok) then
     write(*, '(a, i1, a, i2)')'Highs_call run status is not ', runstatus, ' but ', runstatus
     stop
  endif
  write(*, '(a, i1, a, i2)')'Run status = ', runstatus, '; Model status = ', modelstatus
      
  if (modelstatus .eq. 9) then
     objective_function_value = 0
     ! Report the column primal and dual values, and basis status
     do col = 1, numcol
        write(*, '(a, i1, a, f10.4, a, f10.4, a, i2)') &
             'Col', col, ' = ', colvalue(col), &
             '; dual = ', coldual(col), &
             '; status = ', colbasisstatus(col)
        objective_function_value = objective_function_value + colvalue(col)*colcost(col)
     enddo
     ! Report the row primal and dual values, and basis status
     do row = 1, numrow
        write(*, '(a, i1, a, f10.4, a, f10.4, a, i2)') &
             'Row', row, ' = ', rowvalue(row), &
             '; dual = ', rowdual(row), &
             '; status = ', rowbasisstatus(row)
     enddo
     write(*, '(a, f10.4)')'Optimal objective value = ', objective_function_value
  endif

  ! Define the constraint matrix row-wise, as it is added to the LP with the rows
  arstart(1) = 0
  arstart(2) = 1
  arstart(3) = 3
  arindex(1) = 1
  arindex(2) = 0
  arindex(3) = 1
  arindex(4) = 0
  arindex(5) = 1
  arvalue(1) = 1
  arvalue(2) = 1
  arvalue(3) = 2
  arvalue(4) = 2
  arvalue(5) = 1

  highs = Highs_create()
  call C_F_POINTER(C_NULL_PTR, double_null, [0])
  call C_F_POINTER(C_NULL_PTR, integer_null, [0])

  ! Add two columns to the empty LP
  runstatus = Highs_addCols(highs, numcol, colcost, collower, colupper, 0, integer_null, integer_null, double_null);
  ! Add three rows to the 2-column LP
  runstatus = Highs_addRows(highs, numrow, rowlower, rowupper, numnz, arstart, arindex, arvalue)

  runstatus = Highs_getObjectiveSense(highs, sense);
  write(*, '(a, i2)')"LP problem has objective sense = ", sense
  call assert(sense .eq. 1, "Objective sense")
  
  sense = -1 * sense
  runstatus = Highs_changeObjectiveSense(highs, sense);
  runstatus = Highs_getObjectiveSense(highs, sense);
  call assert(sense .eq. -1, "Changed Objective sense")
  
  runstatus = Highs_getIntOptionValue(highs, "simplex_scale_strategy", simplex_scale_strategy);
  call assert(simplex_scale_strategy .eq. default_simplex_scale_strategy,&
       "simplex_scale_strategy .eq. default_simplex_scale_strategy")
  runstatus = Highs_setIntOptionValue(highs, "simplex_scale_strategy", new_simplex_scale_strategy)
  runstatus = Highs_getIntOptionValue(highs, "simplex_scale_strategy", simplex_scale_strategy);
  call assert(simplex_scale_strategy .eq. new_simplex_scale_strategy,&
       "simplex_scale_strategy .eq. new_simplex_scale_strategy")

  ! There are some functions to check what type of option value you should provide.
  runstatus = Highs_getOptionType(highs, "simplex_scale_strategy", option_type);
  call assert(runstatus .eq. 0, "getHighsOptionType runstatus")
  call assert(option_type .eq. 1, "getHighsOptionType option_type")
  runstatus = Highs_getOptionType(highs, "bad_option", option_type)
  call assert(runstatus .eq. 2, "getHighsOptionType runstatus")

  if (1 .eq. 1) then
  runstatus = Highs_run(highs);
  ! Get the model status
  modelstatus = Highs_getModelStatus(highs, scaled_model);

  write(*, '(a, i1, a, i2)')'Run status = ', runstatus, '; Model status = ', modelstatus
  write(*, '(a, i1, a, a)')'Run status = ', runstatus, '; Model status = ', Highs_modelStatusToChar(highs, modelstatus)

  runstatus = Highs_getDoubleInfoValue(highs, "objective_function_value"//C_NULL_CHAR, objective_function_value);
  runstatus = Highs_getIntInfoValue(highs, "simplex_iteration_count"//C_NULL_CHAR, simplex_iteration_count);
  runstatus = Highs_getIntInfoValue(highs, "primal_status"//C_NULL_CHAR, primal_status);
  runstatus = Highs_getIntInfoValue(highs, "dual_status"//C_NULL_CHAR, dual_status);

  write(*, '(a, f10.4, a, i6)')"Objective value = ", objective_function_value, "; Iteration count = ", simplex_iteration_count
  if (modelstatus .eq. 9) then
!    print*, "Solution primal status = %s\n", Highs_primalDualStatusToChar(highs, primal_status)
!    printf("Solution dual status = %s\n", Highs_primalDualStatusToChar(highs, dual_status));
    ! Get the primal and dual solution
    runstatus = Highs_getSolution(highs, colvalue, coldual, rowvalue, rowdual);
    ! Get the basis
    runstatus = Highs_getBasis(highs, colbasisstatus, rowbasisstatus);
     ! Report the column primal and dual values, and basis status
     do col = 1, numcol
        write(*, '(a, i1, a, f10.4, a, f10.4, a, i2)') &
             'Col', col, ' = ', colvalue(col), &
             '; dual = ', coldual(col), &
             '; status = ', colbasisstatus(col)
     enddo
     ! Report the row primal and dual values, and basis status
     do row = 1, numrow
        write(*, '(a, i1, a, f10.4, a, f10.4, a, i2)') &
             'Row', row, ' = ', rowvalue(row), &
             '; dual = ', rowdual(row), &
             '; status = ', rowbasisstatus(row)
     enddo
  endif
endif

end program fortrantest

subroutine assert ( logic, message)
  logical logic
  character*(*) message
  if (.not.logic) then
     write(*, '(a, a)')'assert fail for ', message
     stop
  endif  
end subroutine assert

