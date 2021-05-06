program test_mpi
  use mpi 
  implicit none
  real*8, parameter :: PI25DT = 3.141592653589793238462643
  real*8 :: mypi, pi, h, sum, x, f, a
  integer :: myid, numprocs, i, ierr
  integer*8 :: n

  f(a) = 4.d0 / (1.d0 + a*a)

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  !print *, 'mpi_comm_world=', MPI_COMM_WORLD
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
  print *, 'mpi_comm_world=', MPI_COMM_WORLD
  print *, 'numprocs=', numprocs

  do 
    if (myid .eq. 0) then
        print *, 'Enter the number of intervals: (0 quits)'
        read(*,*) n
    endif 

    call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (n .le. 0) then
      call MPI_FINALIZE(ierr)
      stop
    endif

    h = 1.d0/n 
    sum = 0.d0
    do i = myid+1, n, numprocs
        x = h * (dble(i) - 0.5d0)
        sum = sum + f(x)
    enddo
    mypi = h* sum

    !call MPI_BARRIER(MPI_COMM_WORLD, ierr )
    !if (myid==0) then
    !  print*, 'all process have reached this point'
    !endif 

    call MPI_REDUCE(mypi, pi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (myid .eq. 0) then
        print *, 'pi is ', pi, 'Error is ', abs(pi - PI25DT)
    endif
  enddo 

end program test_mpi