module type_mod
  implicit none
  type meshtype
    integer*8 :: type
    real*8:: number
  end type meshtype
  
  type(meshtype), parameter :: cell=meshtype(0, 0) ! initializer
  type(meshtype), parameter :: node=meshtype(1, 1)
  
  contains
  subroutine print_type()
    print *, cell%type, node%type
    print *, cell%number, node%number
  end subroutine print_type
end module type_mod

program test
  use type_mod
  implicit none
  
  call print_type()
end program test
