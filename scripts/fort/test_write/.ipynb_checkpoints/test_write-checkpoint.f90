program test
  implicit none
  
  real*8 :: term1, term2, term3, term4
  integer :: it, myid, call_pass
  character(len=240) :: filename
  
  myid = 0
  it = 100
  call_pass = 1
  term1 = -3.2223223E-009
  term2 = 6.4399E-010
  term3 = -2.434343E-010
  term4 = 1.9665E-013
  
  write(filename,"(a,i4.4,a)") 'ecal_ez_pass1_', myid, '.dat' 
  print*, filename
  open(unit=101,file=filename,status='unknown')
  !write(int(100+call_pass), '(I6,E10.6,E10.6,E10.6,E10.6)') it, &
  write(int(100+call_pass), '(I6,1x,6E14.6,1x)') it, &
          term1, term2, term3, term4
  close(unit=101)
  
  write(filename,"(a,i4.4,a)") 'ecal_ez_pass2_', myid, '.dat' 
  print*, filename
  
end program test
