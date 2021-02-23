      program gda_to_pvd
        ! convert binary output file *.gda to *.pvd for ParaView
        use VTR
        use mpi
        type(VTR_file_handle) :: fd

        real(kind=4), dimension(240) :: x,y,z
        real(kind=4), dimension(240,240,240) :: ex,ey,ez,bx,by,bz
        integer :: it, ierr, numprocs, myid
        character(len=128) :: filename
        integer :: fh

        call MPI_INIT(ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
        write(*,*) 'numprocs=',numprocs, 'myid=',myid
        do i=1,240
          x(i)=i*1.0
          y(i)=i*1.0
          z(i)=i*1.0
        enddo
        do it=16000,16000,1000
          write(filename,'(a,i0,a)')'../../data/ex_',it,'.gda'
          call MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
          if (myid .eq. 0) then
            call MPI_File_read(fh, ex, 240*240*240, MPI_REAL, MPI_STATUS_IGNORE, ierr)
            if (ierr .eq. MPI_SUCCESS) then
              write(*,*) ex(1,1,1),ex(240,240,240)
              call VTR_open_file(PREFIX='ex', FD=fd)
              call VTR_write_mesh(FD=fd, X=x, Y=y, Z=z)
              call VTR_write_var(FD=fd, NAME="Ex", FIELD=ex)
              call VTR_close_file(FD=fd)
            else
              write(*,*) 'Read error, code=',ierr
            endif
            call VTR_collect_file(fd)
          endif
          call MPI_File_close(fh,ierr)
        enddo
        call MPI_FINALIZE(ierr)

      end program
