program test_write_files
    implicit none
    real(kind=8), dimension(5):: bx
    character(len=160):: data_directory, filename
    
    bx = [1., 2., 3., 4., 5.]
    data_directory="./data"
    filename = trim(adjustl(data_directory))//'bx_1000.gda'


    
    ! open (1,                                                                         &
    !     file= trim(trim(adjustl(data_directory)))//'/bx_1000.dat', &
    !     form='unformatted',                                                                   &
    !     action='write',access='direct', status='unknown',recl=5)
    
    ! write (1) bx
    print *, bx
    print *, data_directory
    print *, filename

end program test_write_files