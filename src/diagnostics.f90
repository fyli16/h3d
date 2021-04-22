! VR: wrapper for user diagnotics
! VR: eventually, it will be transitioned to 
! VR: modular layout similar to other modern codes,
! VR: where the "engine" is separate from output, etc

!---------------------------------------------------------------------
! user diagnostics
!---------------------------------------------------------------------
module m_diagnostics
  use m_parameters
  implicit none 

  contains 
  
  subroutine diagnostics
    integer :: i
    ! convert 'it' to char and broadcast to all ranks,
    ! which will be used in file dumps by rank.
    if (myid==0) then
      my_short_int = it
      call integer_to_character(cycle_ascii, len(cycle_ascii), my_short_int)
      cycle_ascii_new=trim(adjustl(cycle_ascii))
    endif
    call MPI_BCAST(cycle_ascii,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(cycle_ascii_new,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
    
    ! write mesh data
    if ( n_write_mesh>0 .and. mod(it,n_write_mesh)==0 ) then
      ! this block is not executed when MPI_IO_format=.true.
      if (myid==0 .and. .not.MPI_IO_format) then
        call open_files
      endif 
      
      ! calculate par & perp temperatures (needed only for diagnostics)
      call date_and_time(values=time_begin(:,6))
      if (ndim /= 1) then
        call caltemp2_global
      else
        call caltemp2_global_2d
      endif
      call date_and_time(values=time_end(:,6))
      call accumulate_time(time_begin(1,6),time_end(1,6),time_elapsed(6))

      ! write data
      call write_mesh_data
      
      ! this block is not executed when MPI_IO_format=.true.
      if (myid==0 .and. .not.MPI_IO_format) then
        do i = 1, 25
          close(file_unit(i))
        enddo
      endif

    endif 

    ! write history data: energy.dat
    call diag_energy_hist
      
    ! write particles (within a volume)
    if (n_write_particle>0 .and. mod(it,n_write_particle)==0) then
      call write_particle_in_volume
    endif

    ! write probe data
    call date_and_time(values=time_begin(:,31))
    if ( n_write_probes>0 .and. mod(it,n_write_probes)==0 ) then
      if (tracking_mpi) then
        call virtual_probes_mpi
      else
        call virtual_probes
      endif
    endif 
    call date_and_time(values=time_end(:,31))
    call accumulate_time(time_begin(1,31),time_end(1,31),time_elapsed(31))

    ! write particle tracking data
    call date_and_time(values=time_begin(:,32))
    if ( n_write_tracking>0 .and. mod(it,n_write_tracking)>0 ) then
      if (tracking_mpi) then
        call track_particles_mpi
      else
        call track_particles
      endif
    endif 
    call date_and_time(values=time_end(:,32))
    call accumulate_time(time_begin(1,32),time_end(1,32),time_elapsed(32))

    ! write restart files
    if ( n_write_restart>0 .and. it>itstart .and. mod(it,n_write_restart)==0 ) then
      itrestart = it
      do i = 0, nprocs_over_60  
        if ( mod(int(myid,8),nprocs_over_60+1).eq.i) then
          call write_read_restart_files(1.0)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,IERR)
      enddo

      ! write the latest restart dump info into file
      if (myid == 0) then
        open(unit=222, file=trim(restart_directory)//'restart_index.dat', status='unknown')
        write(222,*) restart_index, itrestart
        close(222)
      endif

      ! keep only two sets of restart files (overwrite previous dump by swapping index)
      if (restart_index == 1) then
        restart_index=2
      else
        restart_index=1
      endif

    endif 

  end subroutine diagnostics


  subroutine diag_energy_hist
    if ( myid==0 .and. n_write_energy>0 .and. mod(it,n_write_energy)==0 ) then
      write(11,*) it, efld, bfld, efluid, ethermal, eptcl ! energy.dat
      ! write(14,*) it, time_elapsed(1:40)  ! time.dat 
    endif
  end subroutine diag_energy_hist


  !---------------------------------------------------------------------
  subroutine virtual_probes
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
    if (i==nbufsteps .or. it==itrestart) then
      if (myid==0) then
        do j=1,nbufsteps
          do k=1,nprobes
            buf2(k,j)=buf(k,j)*factor(k)
          enddo
        enddo
        ! receive data from other processes
        do m=1, nprocs-1
          call MPI_Recv(buf, bufsize, MPI_DOUBLE, m, 1, MPI_COMM_WORLD, status, ierr)
          do j=1,nbufsteps
            do k=1,nprobes
              buf2(m*nprobes+k,j)=buf(k,j)*factor(k)
            enddo
          enddo
        enddo
        do j=1,nbufsteps
          write(12,'(I6,1x)',advance='no')buftime(j)
          do k=1,nprobes*nprocs
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
  ! write field probe data by MPI rank
  !---------------------------------------------------------------------
  subroutine virtual_probes_mpi

    integer :: i, j, k, m, bufsize
    double precision, dimension(nprobes) :: factor

    factor = (/wpiwci**2, wpiwci**2, wpiwci**2, wpiwci, wpiwci, wpiwci/)
    write(12,'(I6,1x,6E14.6,1x)') it, &
          ex(2,jb,kb)*factor(1), ey(2,jb,kb)*factor(2), ez(2,jb,kb)*factor(3), &
          bx(2,jb,kb)*factor(4), by(2,jb,kb)*factor(5), bz(2,jb,kb)*factor(6)

  end subroutine virtual_probes_mpi


  !---------------------------------------------------------------------
  subroutine track_particles
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
      ! receive particles from rank 1 to nprocs-1
      do m=1, nprocs-1
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
  ! write tracking particles by MPI rank
  !---------------------------------------------------------------------
  subroutine track_particles_mpi

    write(13) it, ntot
    write(13) buf_p1(:,1:ntot)

  end subroutine track_particles_mpi


  !---------------------------------------------------------------------
  ! wrapper for user disganostic restart framework
  ! passes the unit to write to
  !---------------------------------------------------------------------
  subroutine user_data_write_restart(wunit)
    integer, intent (in) :: wunit
    
  end subroutine user_data_write_restart
    
  !---------------------------------------------------------------------
  ! wrapper for user disganostic restart framework
  ! passes the unit to read from
  !---------------------------------------------------------------------
  subroutine user_diagnostics_restart(wunit)
    integer, intent (in) :: wunit

  end subroutine user_diagnostics_restart

end module m_diagnostics