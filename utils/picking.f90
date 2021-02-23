      program picking
        use mpi
        implicit none
        integer, parameter:: tracking_width=14, top=200
        integer:: npar, nt, npes
        integer :: ierr, numprocs, myid, io
        double precision, dimension(:,:,:), allocatable :: ptl
        double precision, dimension(:,:), allocatable :: oneline
        integer :: i,j,k, head, len, it, ntot
        character(len=128) :: filename
        double precision, dimension(top) :: energy
        integer, dimension(top) :: tags
        double precision :: e

        nt=65000
        npar=1600
        npes=1600
        allocate(ptl(tracking_width,nt,npar))
        open(unit=14,file='../data/particles.dat',form='unformatted',status='old',iostat=io)
        if (io/=0) then
          write(*,*)'Error opening file: particles.dat'
          return
        endif
        read(14)ptl
        close(14)
        if (.True.) then
          energy=0d0
          tags=0
          do i =1, npar
            ! final velocity
            e=ptl(4,nt,i)**2+ptl(5,nt,i)**2+ptl(6,nt,i)**2
            !write(*,*)i,e
            do j=1, top
              if (e.gt.energy(j))then
                ! insert a particle
                do k=top, j+1,-1
                  energy(k)=energy(k-1)
                  tags(k)=tags(k-1)
                enddo
                energy(j)=e
                tags(j)=i
                exit
              endif
            enddo
          enddo
          do j=1, top
            k=tags(j)
            write(filename,'(a,i4.4,a)')'../data/p_',k,'.dat'
            write(*,*)filename
            write(*,*)'energy=',energy(j)
            open(unit=15,file=filename,form='unformatted',status='unknown',iostat=io)
            write(15)ptl(:,:,k)
            close(15)
          enddo
          ! select particle of highest energy
        else
          ! random selection
          !do i=1,50
          do i=1583,1583
            write(filename,'(a,i4.4,a)')'../data/p_',i,'.dat'
            open(unit=15,file=filename,form='unformatted',status='unknown',iostat=io)
            write(15)ptl(:,:,i)
            close(15)
          enddo
        endif

        !call MPI_FINALIZE(IERR)
      end program
