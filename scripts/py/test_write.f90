program test_write
    implicit none
    character (len=160):: filename
    integer :: myid
    myid = 2
    write(filename,"(a,i4.4,a)")'tracking_',myid,'.dat'
    print*, filename
end program test_write