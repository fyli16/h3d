!VR*****************************************************
!VR: wrapper for user diagnotics
!VR: eventually, it will be transition the whole code to 
!VR: modular layout similar to other modern codes,
!VR: where the "engine" is separate from output, etc
!VR*****************************************************

!---------------------------------------------------------------------
subroutine user_diagnostics
  use parameter_mod, only: time_begin_array, time_end_array, time_elapsed, tracking_mpi
  
  call date_and_time(values=time_begin_array(:,31))
  if (tracking_mpi) then
    call virtual_probes2
  else
    call virtual_probes
  endif
  call date_and_time(values=time_end_array(:,31))
  call accumulate_time_difference(time_begin_array(1,31),time_end_array(1,31),time_elapsed(31))

  call date_and_time(values=time_begin_array(:,32))
  if (tracking_mpi) then
    call track_particles2
  else
    call track_particles
  endif
  call date_and_time(values=time_end_array(:,32))
  call accumulate_time_difference(time_begin_array(1,32),time_end_array(1,32),time_elapsed(32))
end subroutine user_diagnostics


!---------------------------------------------------------------------
subroutine virtual_probes
  use parameter_mod
  implicit none
  integer :: i,j,k,m,bufsize
  double precision, dimension(nprobes) :: factor
  i=mod(it-1, nbufsteps)+1 ! send/receive data every 200 steps
  bufsize=nprobes*nbufsteps
  factor=(/wpiwci**2,wpiwci**2,wpiwci**2,wpiwci,wpiwci,wpiwci/)
  buf(1,i)=ex(2,jb,kb)
  buf(2,i)=ey(2,jb,kb)
  buf(3,i)=ez(2,jb,kb)
  buf(4,i)=bx(2,jb,kb)
  buf(5,i)=by(2,jb,kb)
  buf(6,i)=bz(2,jb,kb)
  buftime(i)=it  ! TBF: in fact, only run 0 needs this
  if (i==nbufsteps .or. it==itfin) then
    if (myid==0) then
      do j=1,nbufsteps
        do k=1,nprobes
          buf2(k,j)=buf(k,j)*factor(k)
        enddo
      enddo
      ! receive data from other processes
      do m=1, npes-1
        call MPI_Recv(buf, bufsize, MPI_DOUBLE, m, 1, MPI_COMM_WORLD, status, ierr)
        do j=1,nbufsteps
          do k=1,nprobes
            buf2(m*nprobes+k,j)=buf(k,j)*factor(k)
          enddo
        enddo
      enddo
      do j=1,nbufsteps
        write(12,'(I6,1x)',advance='no')buftime(j)
        do k=1,nprobes*npes
          write(12,'(E14.6,1x)',advance='no')buf2(k,j)
        enddo
        write(12,*)
      enddo
    else
      ! send data to root 
      call MPI_Send(buf, bufsize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, ierr)
    endif
  endif
end subroutine virtual_probes


!---------------------------------------------------------------------
subroutine virtual_probes2
  use parameter_mod
  implicit none
  integer :: i,j,k,m,bufsize
  double precision, dimension(nprobes) :: factor
  factor=(/wpiwci**2,wpiwci**2,wpiwci**2,wpiwci,wpiwci,wpiwci/)
  write(12,'(I6,1x,6E14.6,1x)')it, ex(2,jb,kb)*factor(1),ey(2,jb,kb)*factor(2),ez(2,jb,kb)*factor(3),&
    bx(2,jb,kb)*factor(4),by(2,jb,kb)*factor(5),bz(2,jb,kb)*factor(6)
end subroutine


!---------------------------------------------------------------------
subroutine track_particles
  use parameter_mod
  implicit none
  double precision, dimension(tracking_width,nspec*maxtags) :: buf_p2
  integer :: iix,iiy,iiz,is, isize, m, n,l,i,k,ii
  if (myid==0) then
    i=mod(it-1, nbufsteps)+1 
    ! copy particles in rank 0
    do n=1,ntot
      k=INT(buf_p1(8,n)) ! this is particle tag
      do l=1,tracking_width
        buf_particle(l,k,i)=buf_p1(l,n)
      enddo
    enddo
    ! receive particles from rank 1 to npes-1
    do m=1, npes-1
      call MPI_Recv(isize, 1, MPI_INTEGER, m, 999, MPI_COMM_WORLD, status, ierr)
      !write(*,*)'receiving',m,isize
      call MPI_Recv(buf_p2, isize*tracking_width, MPI_DOUBLE, m, m, MPI_COMM_WORLD, status, ierr)
      do n=1,isize
        k=INT(buf_p2(8,n))
        do l=1,tracking_width
          buf_particle(l,k,i)=buf_p2(l,n)
        enddo
      enddo
    enddo
    if (i==nbufsteps) then
      if (tracking_binary) then
      !!! binary output
        write(13)buf_particle
      else
      !!! formatted output
      ! TODO: make this following compact
        do ii=1,nbufsteps
          do k=1, nspec*maxtags
            do l=1,7
              write(13,'(E14.6,1x)',advance='no')buf_particle(l,k,ii)
            enddo
            write(13,'(I8)',advance='no')INT(buf_particle(8,k,ii))
            do l=9,14
              write(13,'(E14.6,1x)',advance='no')buf_particle(l,k,ii)
            enddo
            write(13,*)
          enddo
        enddo
      endif
    endif

  else
    !write(*,*)'sending',myid,ntot
    call MPI_Send(ntot, 1, MPI_INTEGER, 0, 999, MPI_COMM_WORLD, ierr)
    call MPI_Send(buf_p1, ntot*tracking_width, MPI_DOUBLE, 0, myid, MPI_COMM_WORLD, ierr)
  endif

end subroutine track_particles


!---------------------------------------------------------------------
subroutine track_particles2
  ! using MPI_IO, one file only
  use parameter_mod
  implicit none
  integer :: n, k, offset
  !offset=(it-1)*tracking_width*maxtags_pe*npes*nspec
  !write(*,*)myid, it, ntot,offset
  !do n=1, ntot
  !  k=INT(buf_p1(8,n)) ! this is particle tag
  !  offset=offset+(k-1)*tracking_width
  !  call MPI_File_write_at(tracking_fh, offset, buf_p1(1,n), tracking_width, MPI_DOUBLE, status, ierr)
  !enddo
  write(13)it, ntot
  write(13)buf_p1(:,1:ntot)

end subroutine track_particles2


!---------------------------------------------------------------------
! wrapper for user disganostic restart framework
! passes the unit to write to
subroutine user_data_write_restart(wunit)
  implicit none
  integer, intent (in) :: wunit
end subroutine user_data_write_restart


!---------------------------------------------------------------------
! wrapper for user disganostic restart framework
! passes the unit to read from
subroutine user_diagnostics_restart(wunit)
  implicit none
  integer, intent (in) :: wunit
end subroutine user_diagnostics_restart
