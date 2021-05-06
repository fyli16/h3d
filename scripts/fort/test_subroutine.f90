module solve_mod
  contains
  subroutine solve(a, b, c, r1, r2)
    implicit none
    real, intent(in) :: a, b, c
    real, intent(inout) :: r1, r2
    real :: term, a2
    term = b**2. - 4.*a*c
    a2 = 2.*a
    if (term>0) then
      term = sqrt(term)
      r1 = (-b+term)/a2
      r2 = (-b-term)/a2
    else
      print *, 'b^2-4ac<0; abort!'
      stop
    endif
  end subroutine solve
end module solve_mod


program test
  use solve_mod
  implicit none
  real :: a, b, c, r1, r2
  a=1; b=5; c=1
  call solve(a, b, c, r1, r2)
  print *, 'roots are: ', r1, ' ', r2
  call solve(a, b, c, r1, r2)
  print *, 'roots are: ', r1, ' ', r2
end program test
    
