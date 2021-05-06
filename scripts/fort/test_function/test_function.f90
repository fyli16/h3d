module ke
  contains
  real function kinetic_energy( v )
    real, dimension(:), intent(in) :: v
    integer i
    real :: ke = 0.0
    !real :: ke
    !ke = 0.0
    do i = 1, size(v)
      ke = ke + v(i)**2
    enddo
    kinetic_energy = .5*ke
  end function kinetic_energy
end module ke

program test
  use ke
  implicit none
  real, dimension(:), allocatable :: v
  real :: ene
  
  v=(/1.0, 2.0/)
  ene = kinetic_energy(v)
  print*, ene
  ene = kinetic_energy(v)
  print*, ene
end program test


