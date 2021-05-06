program test
  implicit none
  integer, dimension(3):: a, b
  double precision, dimension(3):: c, d
  integer:: i, j, k
  integer*4:: v1
  integer*8:: v2
  real*8:: r1, twopi
  double precision:: r2
  character(len=5):: cha
  logical :: periods(2)
  integer*8 :: c2(3), c3(4)
  real :: r3(4)
  character :: s*8, timeunit*4  ! string 's' which has a length of 8
  character(len=120) :: data_directory
  character :: data_dir*5
  real, dimension(2,4) :: arr1
  integer:: it
  real:: time, delta_t, tot_t
  real*8:: single_prec, double_prec
  integer*8:: recl_for_double_precision, recl_for_real
  integer, dimension(:), allocatable :: ix
  
  it=1900; delta_t=0.45887; tot_t=190.23897
  


  a = (/1, 2, 32/)
  b = (/(i, i=4,6)/)
  c = (/1, 2, 32/)
  c2 = (/1, 2, 32/)
  !v1 =1_4; v2=2
  cha = "hello"
  periods = .true.;
  twopi = 2.0*acos(-1.0)
  k=32/512
  c3 = 1
  r3 = c3/2.0
  s = 'hello,world'
  i = 3456
  data_directory='data/'
  data_dir='data/'
  double_prec = 0.
  single_prec = 0.
  inquire (IOLENGTH=recl_for_double_precision) double_prec
  inquire (IOLENGTH=recl_for_real) single_prec
  print*, 'single_prec = ', single_prec, ', double_prec = ', double_prec

  print*, 'a =', a
  print*, b
  print*, c
  print*, 'huge(integer*4)=', huge(v1), ', huge(integer*8)=', huge(v2)  ! max range of the data type
  print*, huge(r1), huge(r2) ! this proves 'real*8' and 'double precision' are of the same type
  print*, periods
  print*, c2
  print*, twopi
  write(6,"(A8,I2,A8,I2,A8,I2)") 'c2(1)=', c2(1), 'c2(2)=', c2(2), 'c2(3)=', c2(3)
  print*, k
  print*, size(c3)
  print*, 'mod(1,2) = ', mod(1,2)
  !print*, 'mod(2,0) = ', mod(2,0)
  print*, c3
  print*, s
  write(timeunit, "(i4.4)") 501
  print*, timeunit
  print*, int(i,8), int(i,8)
  print*, trim(data_directory)//'den/den_001.gda'
  print*, trim(adjustl(data_directory))//'den/den_001.gda'
  print*, data_dir//'den/den_001.gda'
  print*, v2
  print*, arr1
  write(6,"(A3,I6,A8,F6.2,A11,F6.2,A9,F6.2)") 'it=', it, 'time=', time, 'delta_t=', delta_t, 'tot_t=', tot_t
  
  v1 = 3
  v2 = 100
  print*, 'mod(v1, v2) = ', mod(v1, v2)
  if (v2>0 .and. mod(200,v2)==0) print*, 'mod(3,v2) = ', mod(3,v2)
  
  allocate(ix(4))
  ix(:)=0
  print*, ix

  !print*, cha


end program test
