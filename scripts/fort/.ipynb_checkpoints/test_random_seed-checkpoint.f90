program test
  implicit none
  integer :: seed_size
  integer, allocatable :: seed(:)
  real :: ranval(4)
  
  call random_seed(size=seed_size)
  print*, seed_size
  
  allocate(seed(seed_size))
  
  call random_seed(get=seed)
  print*, '1 : ', seed

  call random_seed(get=seed)
  print*, '2 : ', seed
  
  seed=31415
  print*, '3 : ', seed
  
  call random_seed(put=seed) 
  print*, '4 : ', seed
  
  call random_seed(get=seed) 
  print*, '5 : ', seed
  
  call random_number(harvest=ranval)
  print*, ranval
  
  
  deallocate(seed)           ! safe

end program test

