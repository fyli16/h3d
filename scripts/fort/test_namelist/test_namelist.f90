program test
    implicit none
    real :: a,b
    integer :: i,j
    
    namelist /input/ a,b,i,j
    open(10,file='inp.f90', status='old')
    read(10,nml=input)
    write(*,nml=input)
    close(10)
    
    print*, a+b
end program test
