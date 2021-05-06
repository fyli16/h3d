program test
    use parameters
    implicit none
    integer :: ny, jb, je, myid
    integer, dimension(2):: dims
    
    ny=4
    nz=224
    dims(1) = 2
    dims(2) = 16
    myid = 0
    
    call decomp1d(nz, dims(2), myid, jb, je)
    print*, jb, je
    !print*, int(dims(2), 8)
    !print*, mod(nz, dims(2))
    
    print*, nz
end program test


