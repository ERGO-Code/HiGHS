program fortrantest
  use, intrinsic :: iso_c_binding
  use highs_fortran_api
  implicit none

  ! This illustrates the use of Highs_lpCall, the simple F90 interface to
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
  ! After a successful call to Highs_lpCall, the primal and dual
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
  ! The use of Highs_lpCall is illustrated for the LP/MIP example
  !
  ! Min    f  = 2x_0 + 3x_1
  ! s.t.                x_1 <= 6
  !       10 <=  x_0 + 2x_1 <= 14
  !        8 <= 2x_0 +  x_1
  ! 0 <= x_0 <= 3; 1 <= x_1
  integer ( c_int ), parameter :: numcol = 2
  integer ( c_int ), parameter :: numrow = 3
  integer ( c_int ), parameter :: numnz = 5
  integer ( c_int ), parameter :: aformat_colwise = 1
  integer ( c_int ), parameter :: sense = 1
  real ( c_double ), parameter :: offset = 0
  
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
  integer ( c_int ), parameter :: modelstatus_optimal = 7
  integer ( c_int ), parameter :: runstatus_error = -1
  integer ( c_int ), parameter :: runstatus_ok = 0
  integer ( c_int ), parameter :: runstatus_warning = -1
  ! For the full API test
  type ( c_ptr ) :: highs
  
  ! For the row-wise matrix
  integer ( c_int ) arstart(numrow)
  integer ( c_int ) arindex(numnz)
  real ( c_double ) arvalue(numnz)
  
  real, parameter :: inf = 1e30
  integer col, row, el
  integer from_el, to_el
  integer iteration_count, primal_solution_status, dual_solution_status
  double precision objective_function_value
  double precision objective_error
  integer option_type
  integer dummy_info

  integer alt_sense
  integer scale_strategy

  integer, parameter :: default_scale_strategy = 2
  integer, parameter :: new_scale_strategy = 3
  double precision, parameter ::  dual_tolerance = 1d-6
  logical ( c_bool ) write_solution_to_file
  integer ( c_int ) write_solution_style

  double precision, pointer :: double_null(:)
  integer, pointer :: integer_null(:)
  character( c_char ) :: file_name(7)
  
  ! Illustrate the solution of a QP
  !
  ! minimize -x_2 - 3x_3 + (1/2)(2x_1^2 - 2x_1x_3 + 0.2x_2^2 + 2x_3^2)
  !
  ! subject to x_1 + x_3 <= 2; x>=0
  !
  ! Solution x_1 = 0.5; x_2 = 5.0; x_3 = 1.5
  integer ( c_int ), parameter :: qp_numcol = 3
  integer ( c_int ), parameter :: qp_numrow = 1
  integer ( c_int ), parameter :: qp_numnz = 2
  integer ( c_int ), parameter :: qp_hessian_numnz = 4
  integer ( c_int ), parameter :: qformat_triangular = 1

  real ( c_double ) qp_colcost(qp_numcol)
  real ( c_double ) qp_collower(qp_numcol)
  real ( c_double ) qp_colupper(qp_numcol)
  real ( c_double ) qp_rowlower(qp_numrow)
  real ( c_double ) qp_rowupper(qp_numrow)
  integer ( c_int ) qp_astart(qp_numcol)
  integer ( c_int ) qp_aindex(qp_numnz)
  real ( c_double ) qp_avalue(qp_numnz)
  integer ( c_int ) qp_qstart(qp_numcol)
  integer ( c_int ) qp_qindex(qp_hessian_numnz)
  real ( c_double ) qp_qvalue(qp_hessian_numnz)
  
  real ( c_double ) qp_sol(qp_numcol)
  real ( c_double ) qp_colvalue(qp_numcol)
  real ( c_double ) qp_coldual(qp_numcol)
  real ( c_double ) qp_rowvalue(qp_numrow)
  real ( c_double ) qp_rowdual(qp_numrow)
  integer ( c_int ) qp_colbasisstatus(qp_numcol)
  integer ( c_int ) qp_rowbasisstatus(qp_numrow)
  integer qp_modelstatus
  integer qp_runstatus
  logical, parameter :: no_highs_logging = .TRUE.
  
  logical ( c_bool ), parameter :: logical_false = .false.
  logical ( c_bool ), parameter :: logical_true = .true.

  ! Set up the LP/MIP example
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

  qp_sol(1) = 0.5
  qp_sol(2) = 5.0
  qp_sol(3) = 1.5

  !================================================================================
  ! Illustrate use of Highs_lpCall to solve a given LP
  print*, "*********"
  print*, "Section 1"
  print*, "*********"
  runstatus = Highs_lpCall( numcol, numrow, numnz,&
       aformat_colwise, sense, offset, &
       colcost, collower, colupper, rowlower, rowupper,&
       astart, aindex, avalue,&
       colvalue, coldual, rowvalue, rowdual,&
       colbasisstatus, rowbasisstatus, modelstatus)

  if (runstatus .ne. runstatus_ok) then
     write(*, '(a, i1, a, i2)')'Highs_lpCall run status is ', runstatus, ' not ', runstatus_ok
     stop
  endif
  write(*, '(a, i1, a, i2)')'Run status = ', runstatus, '; Model status = ', modelstatus
      
  if (modelstatus .eq. modelstatus_optimal) then
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

  ! Illustrate use of Highs_create() to create a pointer to an
  ! instance of the Highs class, then Highs_passLp to pass LP to
  ! HiGHS, and Highs_run to solve it
  highs = Highs_create()
  if (no_highs_logging) then
     runstatus = Highs_setBoolOptionValue(highs, "output_flag"//C_NULL_CHAR, logical_false)
  endif

  runstatus = Highs_passLp(highs, numcol, numrow, numnz, aformat_colwise, &
  sense, offset, colcost, collower, colupper, rowlower, rowupper, &
  astart, aindex, avalue)
  runstatus = Highs_run(highs)
  modelstatus = Highs_getModelStatus(highs);
  print*, "modelstatus = ", modelstatus
  call assert(runstatus .eq. runstatus_ok, "Highs_run runstatus")
  call assert(modelstatus .eq. modelstatus_optimal, "Highs_run modelstatus optimal")
  write(*, '(a, i1, a, i2)')'Run status = ', runstatus, '; Model status = ', modelstatus
  runstatus = Highs_getDoubleInfoValue(highs, "objective_function_value"//C_NULL_CHAR, objective_function_value);
  runstatus = Highs_getIntInfoValue(highs, "simplex_iteration_count"//C_NULL_CHAR, iteration_count);
  write(*, '(a, f10.4, a, i6)')"Objective value = ", objective_function_value, "; Iteration count = ", iteration_count

  call Highs_destroy(highs)
  !================================================================================
  ! Illustrate use of Highs_addCols and Highs_addRows to build model,
  ! and then Highs_changeObjectiveSense to switch to maximization
  print*, "*********"
  print*, "Section 2"
  print*, "*********"
  highs = Highs_create()
  if (no_highs_logging) then
     runstatus = Highs_setBoolOptionValue(highs, "output_flag"//C_NULL_CHAR, logical_false)
  endif
  ! Create double and integer values equal to NULL pointer
  call C_F_POINTER(C_NULL_PTR, double_null, [0])
  call C_F_POINTER(C_NULL_PTR, integer_null, [0])

  ! Add two columns to the empty LP, but no matrix. After numnz=0, can
  ! just pass arrays rather than NULL
  runstatus = Highs_addCols(highs, numcol, colcost, collower, colupper, 0, integer_null, integer_null, double_null);
  ! Add three rows to the 2-column LP
  runstatus = Highs_addRows(highs, numrow, rowlower, rowupper, numnz, arstart, arindex, arvalue)

  runstatus = Highs_getObjectiveSense(highs, alt_sense);
  write(*, '(a, i2)')"LP problem has objective sense = ", alt_sense
  call assert(alt_sense .eq. sense, "Objective sense")
  
  alt_sense = -1 * alt_sense
  runstatus = Highs_changeObjectiveSense(highs, alt_sense);
  runstatus = Highs_getObjectiveSense(highs, alt_sense);
  call assert(alt_sense .eq. -1, "Changed Objective sense")
  
  ! Get and set option values
  runstatus = Highs_getIntOptionValue(highs, "simplex_scale_strategy"//C_NULL_CHAR, scale_strategy);
  call assert(scale_strategy .eq. default_scale_strategy,&
       "scale_strategy .eq. default_scale_strategy")
  runstatus = Highs_setIntOptionValue(highs, "simplex_scale_strategy"//C_NULL_CHAR, new_scale_strategy)
  runstatus = Highs_getIntOptionValue(highs, "simplex_scale_strategy"//C_NULL_CHAR, scale_strategy);
  call assert(scale_strategy .eq. new_scale_strategy,&
       "scale_strategy .eq. new_scale_strategy")

  runstatus = Highs_setDoubleOptionValue(highs, "primal_feasibility_tolerance"//C_NULL_CHAR, 1d-6);
  call assert(runstatus .eq. runstatus_ok, "setDoubleOptionValue runstatus")
  runstatus = Highs_setDoubleOptionValue(highs, "dual_feasibility_tolerance"//C_NULL_CHAR, dual_tolerance);
  call assert(runstatus .eq. runstatus_ok, "setDoubleOptionValue runstatus")

  ! There are some functions to check what type of option value you should provide.
  runstatus = Highs_getOptionType(highs, "simplex_scale_strategy"//C_NULL_CHAR, option_type);
  call assert(runstatus .eq. runstatus_ok, "getOptionType runstatus = 0")
  call assert(option_type .eq. 1, "getOptionType option_type")
  ! This is what happens if an invalid name is passed
  runstatus = Highs_getOptionType(highs, "bad_option"//C_NULL_CHAR, option_type)
  call assert(runstatus .eq. runstatus_error, "getOptionType runstatus")

  ! Suppress HiGHS output
  print*, "Suppressing all HiGHS output"
  runstatus = Highs_setBoolOptionValue(highs, "output_flag"//C_NULL_CHAR, logical_false)
  ! Solve the LP
  runstatus = Highs_run(highs);
  ! Get the model status
  modelstatus = Highs_getModelStatus(highs);
  write(*, '(a, i1, a, i2)')'Run status = ', runstatus, '; Model status = ', modelstatus

  ! Get solution data
  runstatus = Highs_getDoubleInfoValue(highs, "objective_function_value"//C_NULL_CHAR, objective_function_value);
  runstatus = Highs_getIntInfoValue(highs, "simplex_iteration_count"//C_NULL_CHAR, iteration_count);
  runstatus = Highs_getIntInfoValue(highs, "primal_solution_status"//C_NULL_CHAR, primal_solution_status);
  runstatus = Highs_getIntInfoValue(highs, "dual_solution_status"//C_NULL_CHAR, dual_solution_status);
  ! This is what happens if an invalid name is passed
  runstatus = Highs_getIntInfoValue(highs, "bad_info"//C_NULL_CHAR, dummy_info)
  call assert(runstatus .eq. runstatus_error, "getOptionType runstatus")

  write(*, '(a, f10.4, a, i6)')"Objective value = ", objective_function_value, "; Iteration count = ", iteration_count
  print*, "modelstatus = ", modelstatus
  call assert(modelstatus .eq. modelstatus_optimal, "Optimal => modelstatus = modelstatus_optimal")
  if (modelstatus .eq. modelstatus_optimal) then
     call assert(primal_solution_status .eq. 2, "Optimal => primal_solution_status = 2")
     call assert(dual_solution_status .eq. 2, "Optimal => dual_solution_status = 2")
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

  ! Write out model as MPS for use later
  runstatus = Highs_writeModel(highs, "F90.mps"//C_NULL_CHAR)
  print*, "runstatus = ", runstatus
  call assert(runstatus .ne. runstatus_warning, "Highs_writeModel runstatus")
  
  call Highs_destroy(highs)
  !================================================================================
  ! Illustrate use of Highs_readModel to read model, and
  ! Highs_writeSolution to write the solution
  print*, "*********"
  print*, "Section 3"
  print*, "*********"
  highs = Highs_create()
  if (no_highs_logging) then
     runstatus = Highs_setBoolOptionValue(highs, "output_flag"//C_NULL_CHAR, logical_false)
  endif

  ! Read the LP
  runstatus = Highs_readModel(highs, "F90.mps"//C_NULL_CHAR)
  call assert(runstatus .eq. runstatus_ok, "Highs_readModel runstatus")
  ! Solve the LP
  runstatus = Highs_run(highs);
  ! Get the model status
  modelstatus = Highs_getModelStatus(highs);
  write(*, '(a, i1, a, i2)')'Run status = ', runstatus, '; Model status = ', modelstatus
  runstatus = Highs_getDoubleInfoValue(highs, "objective_function_value"//C_NULL_CHAR, objective_function_value);
  runstatus = Highs_getIntInfoValue(highs, "simplex_iteration_count"//C_NULL_CHAR, iteration_count);
  write(*, '(a, f10.4, a, i6)')"Objective value = ", objective_function_value, "; Iteration count = ", iteration_count

  ! Write solution to the screen
  runstatus = Highs_writeSolutionPretty(highs, ""//C_NULL_CHAR)
  call Highs_destroy(highs)
  
  !================================================================================
  ! Illustrate use of setting bool options and string options
  ! (solution_file) so only run(highs) is required
  print*, "*********"
  print*, "Section 4"
  print*, "*********"
  highs = Highs_create()
  if (no_highs_logging) then
     runstatus = Highs_setBoolOptionValue(highs, "output_flag"//C_NULL_CHAR, logical_false)
  endif

  ! Get and set string options
  runstatus = Highs_getStringOptionValue(highs, "solution_file"//C_NULL_CHAR, file_name)
  print*, "Default solution_file is |", file_name, "|"  
  runstatus = Highs_setStringOptionValue(highs, "solution_file"//C_NULL_CHAR, "F90.sol"//C_NULL_CHAR)
  runstatus = Highs_getStringOptionValue(highs, "solution_file"//C_NULL_CHAR, file_name)
  print*, "New solution_file is |", file_name, "|"  

! Get and set bool options. NB Cannot pass .true. as it's 4-byte
  runstatus = Highs_getBoolOptionValue(highs, "write_solution_to_file"//C_NULL_CHAR, write_solution_to_file)
  print*, "Default write_solution_to_file = ", write_solution_to_file
  write_solution_to_file = .true.
  runstatus = Highs_setBoolOptionValue(highs, "write_solution_to_file"//C_NULL_CHAR, write_solution_to_file)

  runstatus = Highs_getIntOptionValue(highs, "write_solution_style"//C_NULL_CHAR, write_solution_style)
  print*, "Default write_solution_style = ", write_solution_style
  write_solution_style = 1;
  runstatus = Highs_setIntOptionValue(highs, "write_solution_style"//C_NULL_CHAR, write_solution_style)

  ! Report all the deviations from default options
  runstatus = Highs_writeOptionsDeviations(highs, "OptionsDeviations.set"//C_NULL_CHAR)
  
  ! Reset all the options
  runstatus = Highs_resetOptions(highs)
  
  ! Report all the options
  runstatus = Highs_writeOptions(highs, "Options.set"//C_NULL_CHAR)

  if (no_highs_logging) then
     runstatus = Highs_setBoolOptionValue(highs, "output_flag"//C_NULL_CHAR, logical_false)
  endif
  
  ! Read the LP
  runstatus = Highs_readModel(highs, "F90.mps"//C_NULL_CHAR)
  ! Solve the LP
  runstatus = Highs_run(highs);
  ! Get the model status
  modelstatus = Highs_getModelStatus(highs);
  call assert(modelstatus .eq. modelstatus_optimal, "Highs_run modelstatus optimal")
  write(*, '(a, i1, a, i2)')'Run status = ', runstatus, '; Model status = ', modelstatus
  runstatus = Highs_getDoubleInfoValue(highs, "objective_function_value"//C_NULL_CHAR, objective_function_value);
  runstatus = Highs_getIntInfoValue(highs, "simplex_iteration_count"//C_NULL_CHAR, iteration_count);
  write(*, '(a, f10.4, a, i6)')"Objective value = ", objective_function_value, "; Iteration count = ", iteration_count
  call Highs_destroy(highs)
  
  !================================================================================
  ! Illustrate use of Highs_qpCall to solve a given QP

  print*, "**********"
  print*, "QP Example"
  print*, "**********"

  qp_colcost(1) = 0
  qp_colcost(2) = -1
  qp_colcost(3) = -3
  
  qp_collower(1) = 0
  qp_collower(2) = 0
  qp_collower(3) = 0
  
  qp_colupper(1) = inf
  qp_colupper(2) = inf
  qp_colupper(3) = inf
  
  qp_rowlower(1) = -inf
  qp_rowupper(1) = 2
  
  qp_astart(1) = 0
  qp_astart(2) = 1
  qp_astart(3) = 1
  
  qp_aindex(1) = 0
  qp_aindex(2) = 0
  
  qp_avalue(1) = 1
  qp_avalue(2) = 1

  qp_qstart(1) = 0
  qp_qstart(2) = 2
  qp_qstart(3) = 3
  
  qp_qindex(1) = 0
  qp_qindex(2) = 2
  qp_qindex(3) = 1
  qp_qindex(4) = 2
  
  qp_qvalue(1) = 2.0
  qp_qvalue(2) = -1.0
  qp_qvalue(3) = 0.2
  qp_qvalue(4) = 2.0

  runstatus = Highs_qpCall( qp_numcol, qp_numrow, qp_numnz, qp_hessian_numnz,&
       aformat_colwise, qformat_triangular, sense, offset,&
       qp_colcost, qp_collower, qp_colupper, qp_rowlower, qp_rowupper,&
       qp_astart, qp_aindex, qp_avalue,&
       qp_qstart, qp_qindex, qp_qvalue,&
       qp_colvalue, qp_coldual, qp_rowvalue, qp_rowdual,&
       qp_colbasisstatus, qp_rowbasisstatus, modelstatus)

  if (runstatus .ne. runstatus_ok) then
     write(*, '(a, i1, a, i2)')'Highs_lpCall run status is ', runstatus, ' not ', runstatus_ok
     stop
  endif
  write(*, '(a, i1, a, i2)')'Run status = ', runstatus, '; Model status = ', modelstatus
      
  if (modelstatus .eq. modelstatus_optimal) then
     objective_function_value = 0
     ! Report the column primal and dual values, and basis status
     do col = 1, qp_numcol
       write(*, '(a, i1, a, f10.4, a, f10.4, a, i2)') &
             'Col', col, ' = ', qp_colvalue(col), &
             '; dual = ', qp_coldual(col)
     call assert(abs(qp_colvalue(col)-qp_sol(col)) .le. 1e-4, "Solution error")
     enddo
     ! Report the row primal and dual values, and basis status
     do row = 1, qp_numrow
        write(*, '(a, i1, a, f10.4, a, f10.4, a, i2)') &
             'Row', row, ' = ', qp_rowvalue(row), &
             '; dual = ', qp_rowdual(row)
     enddo
     do col = 1, qp_numcol
       objective_function_value = objective_function_value + qp_colvalue(col)*qp_colcost(col)
     enddo
     do col = 1, qp_numcol
        from_el = qp_qstart(col)
        if (col < qp_numcol) then
           to_el = qp_qstart(col+1)-1
        else
           to_el = qp_hessian_numnz
        endif
        objective_function_value = &
             objective_function_value + 0.5 * qp_colvalue(col) * qp_qvalue(from_el+1) * qp_colvalue(col)
        do el = from_el+1, to_el
           row = qp_qindex(el+1)+1
           objective_function_value = &
                objective_function_value + qp_colvalue(col) * qp_qvalue(el+1) * qp_colvalue(row)
        enddo
     enddo
     write(*, '(a, f10.4)')'Optimal objective value = ', objective_function_value
     objective_error = abs(objective_function_value+5.25)
     call assert(objective_error .le. 1e-4, "Objective error")
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

