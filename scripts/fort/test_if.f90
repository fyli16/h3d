program test
  implicit none
  integer :: i
  
  i=3
  if (i<0) then
    print*, 'i<0'
  else if (i<5) then
    print*, '0<=i<5'
  else if (i<10) then
    print*, '5<=i<10'
  endif

end program test
