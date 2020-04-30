program fortrantest
  use, intrinsic :: iso_c_binding
  use highs_lp_solver
  implicit none

  integer ( c_int ), parameter :: n = 2
  integer ( c_int ), parameter :: m = 2
  integer ( c_int ), parameter :: nz = 4

  real ( c_double ) colcost(n)
  real ( c_double ) collower(n)
  real ( c_double ) colupper(n)
  real ( c_double ) rowlower(m)
  real ( c_double ) rowupper(m)
  integer ( c_int ) astart(n)
  integer ( c_int ) aindex(nz)
  real ( c_double ) avalue(nz)

  real ( c_double ) colvalue(n)
  real ( c_double ) coldual(n)
  real ( c_double ) rowvalue(m)
  real ( c_double ) rowdual(m)
  integer ( c_int ) colbasisstatus(n)
  integer ( c_int ) rowbasisstatus(m)
  integer ( c_int ) modelstatus
  integer ( c_int ) runstatus

  colcost(1) = 1
  colcost(2) = -2
  collower(1) = 0
  collower(2) = 0
  colupper(1) = 1000
  colupper(2) = 1000
  rowlower(1) = 0.0
  rowlower(2) = 0.0
  rowupper(1) = 10.0
  rowupper(2) = 10.0
  astart(1) = 0
  astart(2) = 2
  aindex(1) = 0
  aindex(2) = 1
  aindex(3) = 0
  aindex(4) = 1
  avalue(1) = 1
  avalue(2) = -1
  avalue(3) = 3
  avalue(4) = 0.2

  runstatus = Highs_call( n, m, nz,&
  colcost, collower, colupper,&
  rowlower, rowupper,&
  astart, aindex, avalue,&
  colvalue, coldual, rowvalue, rowdual,&
  colbasisstatus, rowbasisstatus, modelstatus)
      

end program fortrantest
