program test
  implicit none
  integer :: i
  integer, dimension(3) :: arr1

  arr1 = (/1.0, 2.0, 3.0/)

  do i = 1, 0
    !print*, arr1(i)
    print*, 'hello'
  enddo

  print*, arr1(1), arr1(2), arr1(3)
    
end program test
