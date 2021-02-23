      program tracking
        !use mpi
        implicit none
        integer, parameter:: tracking_width=14
        integer:: npar, nt, npes
        integer :: ierr, numprocs, myid, io
        double precision, dimension(:,:,:), allocatable :: ptl
        double precision, dimension(:,:), allocatable :: oneline
        integer :: i,j,k, head, len, it, ntot
        character(len=128) :: filename

        !call MPI_INIT(IERR)
        !call MPI_COMM_SIZE(MPI_COMM_WORLD, NUMPROCS, IERR)
        !call MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)

        nt=65000
        npar=1600
        npes=1600
        allocate(ptl(tracking_width,nt,npar))
        allocate(oneline(tracking_width,npar))
        ptl=0
        
        do i=0, npes-1
        !do i=npes-20, npes-1
        !do i=1000, 1030
        !do i=16,16
          write(filename,'(a,I4.4,a)')'tracking_',i,'.dat'
          write(*,*)filename
          open(unit=13,file='../data/'//filename,form='unformatted',status='old',iostat=io)
          if (io/=0) then
            write(*,*)'Error opening file:',filename
            return
          endif
          do
            !inquire(unit=13,RECL=len)
            !write(*,*)'record length:',len
            read(13,iostat=io)it, ntot
            if (io .ne. 0) then
              if (io .eq. -1) then ! end of file
                goto 999
              else
                write(*,*)'Read ntot'
                write(*,*)'it=',it
                write(*,*)'I/O error, code:',io
                goto 999
              endif
            else
              !if (ntot>0) then
              !  write(*,*) 'i=',i,'it=',it, 'ntot=',ntot
              !endif
              if (it>65000) goto 999
              read(13,iostat=io)oneline(:,1:ntot)
              if (io .eq. 0) then
                do j=1, ntot
                  k=int(oneline(8,j))
                  !write(*,*)k
                  ptl(:,it,k)=oneline(:,j)
                enddo
              else
                write(*,*)'Read record: it=',it,'ntot=',ntot
                write(*,*)'I/O error, code:',io
              endif
            endif
          enddo
999       continue
          close(13)
        enddo
        open(unit=14,file='../data/particles.dat',form='unformatted',status='unknown',iostat=io)
        if (io/=0) then
          write(*,*)'Error opening file: particles.dat'
          return
        endif
        write(14)ptl(:,:,1:npar)
        close(14)


        !call MPI_FINALIZE(IERR)
      end program
