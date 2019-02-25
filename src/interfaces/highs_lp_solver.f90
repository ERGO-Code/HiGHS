module highs_lp_solver
  interface
    subroutine callhighs ( n, m, nz, cc, cl, cu, rl, ru, as, ai, av ) bind ( c )
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
    end subroutine callhighs
  end interface
end module highs_lp_solver
