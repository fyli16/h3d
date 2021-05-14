program test
    implicit none
    integer :: file_unit(25), file_ref, j

    file_ref = 250
    do j = 1, 25
        file_unit(j) = file_ref + j
    enddo

    print*, file_unit
end program test