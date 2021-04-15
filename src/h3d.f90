!***************************************************************************
!                                                                          *
!                                 VERSION 6.0                              *
!                           YURI'S NONUNIFORM MESH                         *
!                           3D IMPLEMENTATION ONLY                         *
!                      UNIFORM LOADING IN PHYSICAL SPACE                   *
!               UNIFORM LOADING IN LOGICAL SPACE NOT YET IMPLEMENTED       *
!                                                                          *
!***************************************************************************

program h3d 
  use parameter_mod
  use functions_mod
  use mesh_mod
  implicit none

  ! Initialize MPI
  call MPI_INIT(IERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)

  ! Read input file
  call read_input
  ! Decompose MPI/simulation domain  
  call domain_decomp
  ! allocate global parameters
  call allocate_arrays  
    
  ! set up mesh 
  call setup_mesh

  ! open history diagnostic files
  call open_hist_files

  ! restart or a fresh start
  call makelist  ! make list of particles (see utils.f90)
  if (restart) then ! restart from previous run, read wave-particle parameters
    call init_restart 
  else ! fresh start, initialize wave-particle parameters
    call init_wavepart 
  endif 

  ! main simulation loops
  call sim_loops

  ! shutdown the program
  call shutdown

end program h3d


!---------------------------------------------------------------------
! initialize data from a previous run
!---------------------------------------------------------------------
subroutine init_restart
  use parameter_mod
  use mesh_mod
  implicit none 

  integer*8 :: ixe, iye, ize, i, j, k

  if (myid == 0) then
    open(unit=222, file=trim(restart_directory)//'restart_index.dat', status='old')
    read(222,*) restart_index, itfin ! itfin is the final restart dump of previous run
    print*, " "
    print*, " Restart from set # ", restart_index
    print*, " Restart from iteration # ", itfin
    print*, " "
    close(222)
  endif

  call MPI_BCAST(restart_index,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(itfin        ,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)

  do i = 0, nprocs_over_60 
    if (mod(int(myid,8), nprocs_over_60+1) .eq. i) then 
        call restart_read_write(-1.0)  ! read restart data
        call MPI_BCAST(itfin,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
    endif
  enddo

  ! determine start, finish number of steps
  itstart = itfin
  it = itstart
  itfinish = (tmax-t_stopped)/dtwci + itstart
    
  if (restart_index == 1) then
    restart_index=2
  else
    restart_index=1
  endif

  ! Uniform mesh - Same as is in version 5.0
  yb = (jb-1)*hy  ! where is hy, hz defined before this?
  ye = je    *hy
  zb = (kb-1)*hz
  ze = ke    *hz
    
  ! Nonuniform mesh
  zb = meshZ%xn(kb+1)
  ze = meshZ%xn(ke+2)
  do ipe = 0,nprocs-1
    zbglobal(ipe) = meshZ%xn(kbglobal(ipe)+1)
    zeglobal(ipe) = meshZ%xn(keglobal(ipe)+2)
  enddo

  yb = meshY%xn(jb+1)
  ye = meshY%xn(je+2)
  do ipe = 0,nprocs-1
    ybglobal(ipe) = meshY%xn(jbglobal(ipe)+1)
    yeglobal(ipe) = meshY%xn(jeglobal(ipe)+2)
  enddo
            
  volume_fraction = (ye-yb)*(ze-zb)/(ymax*zmax)
  
  xb = 0.
  xe = xmax
  xb_logical = MESH_UNMAP(meshX,xb)
  xe_logical = MESH_UNMAP(meshX,xe)
  yb_logical = MESH_UNMAP(meshY,yb)
  ye_logical = MESH_UNMAP(meshY,ye)
  zb_logical = MESH_UNMAP(meshZ,zb)
  ze_logical = MESH_UNMAP(meshZ,ze)
            
  do i = 1, nspec
    npm = npx(i)*npy(i)*npz(i)*nprocs
    dfac(i) = real(ny*nz*nx)/real(npm)
    do ixe=1,nx2
        do iye=jb-1,je+1
          do ize=kb-1,ke+1
              qp_cell(ixe,iye,ize,i) = meshX%dxc(ixe)*meshY%dxc(iye+1)*meshZ%dxc(ize+1)*dfac(i)*frac(i)
          enddo
        enddo
    enddo
  enddo
          
  do i = 1, nxmax
    xc_uniform(i) = hx*(i-1.5)
    xv_uniform(i) = hx*(i-2.0)
  enddo

  do j = 1, nymax
    yc_uniform(j) = hy*(j-0.5)
    yv_uniform(j) = hy*(j-1.0)
  enddo

  do k = 1, nzmax
    zc_uniform(k) = hz*(k-0.5)
    zv_uniform(k) = hz*(k-1.0)
  enddo
    
  if (myid ==0) then
    my_short_int = it
    call integer_to_character(cycle_ascii,len(cycle_ascii),my_short_int)
    if (cycle_ascii=='') cycle_ascii='0'
    cycle_ascii_new=trim(adjustl(cycle_ascii))
    print*, " cycle = ", cycle_ascii_new
  endif
  
  call MPI_BCAST(cycle_ascii    ,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(cycle_ascii_new,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

  return 

end subroutine init_restart


!---------------------------------------------------------------------
! diagnostic data output
!---------------------------------------------------------------------
subroutine data_output
  use parameter_mod
  implicit none 

  integer :: i

  ! convert 'it' to char and broadcast to all ranks
  if (myid==0) then
    my_short_int = it
    call integer_to_character(cycle_ascii, len(cycle_ascii), my_short_int)
    cycle_ascii_new=trim(adjustl(cycle_ascii))
  endif
  call MPI_BCAST(cycle_ascii,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(cycle_ascii_new,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  
  ! this block is not executed when MPI_IO_format=.true.
  if (myid==0 .and. .not.MPI_IO_format) then
    call open_files
  endif 
  
  ! calculate par & perp temperatures (needed only for diagnostics)
  call date_and_time(values=time_begin_array(:,6))
  if (ndim /= 1) then
    call caltemp2_global
  else
    call caltemp2_global_2d
  endif
  call date_and_time(values=time_end_array(:,6))
  call accumulate_time(time_begin_array(1,6),time_end_array(1,6),time_elapsed(6))

  ! write data
  call dataout
  
  ! this block is not executed when MPI_IO_format=.true.
  if (myid==0 .and. .not.MPI_IO_format) then
    do i = 1, 25
      close(file_unit(i))
    enddo
  endif

  ! write particles (within a given volume)
  if (mod(it,n_write_particle)==0) then
    if (myid == 0) then
      my_short_int = it
      call integer_to_character(cycle_ascii, len(cycle_ascii), my_short_int)
    endif
    call MPI_BCAST(cycle_ascii,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
    call particle_in_volume_write
  endif

  ! write restart files
  if (mod(int(it,8),n_write_restart)==0 .and. (it>itstart)) then
    if (myid == 0) print*, 'Writing restart files ...'
    itfin = it
    do i = 0, nprocs_over_60  
      if ( mod(int(myid,8),nprocs_over_60+1).eq.i) then
        call restart_read_write(1.0)
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    enddo

    ! write the latest restart dump info into file
    if (myid == 0) then
      open(unit=222, file=trim(restart_directory)//'restart_index.dat', status='unknown')
      write(222,*) restart_index, itfin
      close(222)
    endif

    ! keep only two sets of restart files (overwrite previous dump by swapping index)
    if (restart_index == 1) then
      restart_index=2
    else
      restart_index=1
    endif

  endif 

end subroutine data_output


!---------------------------------------------------------------------
! simulation loop
!---------------------------------------------------------------------
subroutine sim_loops
  use parameter_mod
  implicit none 

  integer :: i

  if (myid==0) print*, "Executing main simulation loops:"

  ! time stamp just before entering the simulation loop
  call date_and_time(values=now)
  clock_time_init = now(5)*3600.+now(6)*60.+now(7)+now(8)*0.001
  clock_time_old = clock_time_init
  
  ! main simulation loop
  time_elapsed=0.; time_begin_array=0; time_end_array=0

  do while(it <= itfinish)
  
    call date_and_time(values=time_begin_array(:,1)) ! time one-whole-loop

    ! print time-step info
    call date_and_time(values=now)
    clock_time = now(5)*3600.+now(6)*60.+now(7)+now(8)*0.001
    if (myid == 0.and.mod(it,n_print) == 0) then
      write(6,"(A5,I7,A1,I7,A11,F8.3,A14,F8.3,A12,F8.3)") 'it = ', it, '/', itfinish, ',   time = ', time, &
                    ',   delta_t = ', real(clock_time-clock_time_old), &
                    ',   tot_t = ', real(clock_time-clock_time_init)
      clock_time_old = clock_time
    endif

    ! compute resistivity (which could depend on local parameters such as current)
    call date_and_time(values=time_begin_array(:,2))
    if (ndim /= 1) then
      call eta_calc       ! Dietmar's resistivity
    else
      call eta_calc_2d    ! Dietmar's resistivity
    endif
    ntot = 0 ! for particle tracking
    call trans ! trans computes density and v's; it also calls parmov (particle push)
    call date_and_time(values=time_end_array(:,2))
    call accumulate_time(time_begin_array(1,2),time_end_array(1,2),time_elapsed(2))

    ! write history data: energy.dat, time.dat
    if (myid==0) then
      write(11,*) it, efld, bfld, efluidt, ethermt, eptclt ! energy.dat
      write(14,*) it, time_elapsed(1:40)  ! time.dat
    endif

    ! sort particles
    call date_and_time(values=time_begin_array(:,4))
    if (mod(it,10_8) == 0) call sortit    !  sort the particles
    call date_and_time(values=time_end_array(:,4))
    call accumulate_time(time_begin_array(1,4),time_end_array(1,4),time_elapsed(4))

    ! call field solver
    call date_and_time(values=time_begin_array(:,5))
    if (.not.testorbt) then
      if (ndim /=1) then 
        call field
      else
        call field_2d
      endif
    endif        
    call date_and_time(values=time_end_array(:,5))
    call accumulate_time(time_begin_array(1,5),time_end_array(1,5),time_elapsed(5))

    ! inject_wave ?
    ! if (it == 21000) call inject_wave
    ! if (mod(it,100) == 0) call kick

    ! VR: call user diagnostics
    call date_and_time(values=time_begin_array(:,30))
    call user_diagnostics
    call date_and_time(values=time_end_array(:,30))
    call accumulate_time(time_begin_array(1,30),time_end_array(1,30),time_elapsed(30))
      
    ! data output (mesh data, particles, restart files, etc.)
    if (.not.testorbt.and.mod(it,n_write_data)==0) then        
      call data_output
    endif

    time = time + dtwci
    it = it + 1

    call date_and_time(values=time_end_array(:,1))
    call accumulate_time(time_begin_array(1,1),time_end_array(1,1),time_elapsed(1))

  enddo 

end subroutine sim_loops


!---------------------------------------------------------------------
! shutdown the simulation and exit
!---------------------------------------------------------------------
subroutine shutdown
  use parameter_mod
  implicit none 

  if (myid == 0) then
    close(unit=11) ! energy.dat
    close(unit=12) ! probe.dat
    close(unit=13) ! tracking.dat
    close(unit=14) ! time.dat
  endif

  if (tracking_mpi) close(unit=13)
      
  if (myid==0) then
    print*, " "
    print*, " "
    print*, " *** Run completed *** "
    print*, " "
    print*, " "
    print*, " subroutine trans             (s)          =",time_elapsed(2)
    print*, " subroutine sort              (s)          =",time_elapsed(4)
    print*, " subroutine field             (s)          =",time_elapsed(5)
    print*, " subroutine caltemp2          (s)          =",time_elapsed(6)
    print*, " subroutine user_diagnostics  (s)          =",time_elapsed(30)
    print*, " total time                   (s)          =",time_elapsed(1)
    print*, " "
    print*, " "
    print*, " In subroutine caltemp2,"
    print*, "   subroutine xreal           (s)          =",time_elapsed(24)
    print*, "   subroutine xrealbcc        (s)          =",time_elapsed(25)
    print*, "   interpolation              (s)          =",time_elapsed(26)
    print*, "   total caltemp2             (s)          =",time_elapsed(23)
    print*, " "
    print*, " "
    print*, " In subroutine trans," 
    print*, "   subroutine parmov          (s)          =",time_elapsed(7)
    print*, "   subroutine energy          (s)          =",time_elapsed(8)
    print*, "   total trans                (s)          =",time_elapsed(20)
    print*, " "
    print*, " "
    print*, " In subroutine field,"
    print*, "   subroutine pressgrad       (s)          =",time_elapsed(9)
    print*, "   subroutine bcalc           (s)          =",time_elapsed(10)
    print*, "   subroutine ecalc           (s)          =",time_elapsed(11)
    print*, "   subroutine focalc          (s)          =",time_elapsed(12)
    print*, "   total field                (s)          =",time_elapsed(21)
    print*, " "
    print*, " "
    print*, " In subroutine parmov,"
    print*, "   particle pushing           (s)          =",time_elapsed(13)
    print*, "   particle exchange          (s)          =",time_elapsed(14)
    print*, "   moment calculation         (s)          =",time_elapsed(15)
    print*, "   total parmov               (s)          =",time_elapsed(19)
    print*, " "
    print*, " "
    print*, " In subroutine bcalc,"
    print*, "   subroutine ecalc           (s)          =",time_elapsed(16)
    print*, "   subroutine xrealbcc        (s)          =",time_elapsed(17)
    print*, "   total bcalc                (s)          =",time_elapsed(22)
    print*, " "
    print*, " "
    print*, " In subroutine user_diagnostics,"
    print*, "   sub virtual_probes         (s)          =",time_elapsed(31)
    print*, "   sub track_particle         (s)          =",time_elapsed(32)
    print*, "   total user_diagnose        (s)          =",time_elapsed(30)
    print*, " "
    print*, " "
    print*, " In subroutine ecalc (called from subroutine bcalc and field),"
    print*, "   subroutine xrealbcc        (s)          =",time_elapsed(18)
    print*, " "
    print*, " "
    print*, "communication time/total time (%)          =" &
                ,100.*(time_elapsed(24)+time_elapsed(25)+time_elapsed(14) &
                +time_elapsed(17)+time_elapsed(18))/time_elapsed(1)
    print*, " "
    print*, " "
    print*, "Further breakdown of communication time "
    print*, "  particle exchage in subroutine parmov (%) =" &
                ,100.*time_elapsed(14)/time_elapsed(1)
    print*, "  subroutine xrealbcc                   (%) =" &
                ,100.*(time_elapsed(25)+time_elapsed(17)+time_elapsed(18))/time_elapsed(1)
    print*, "  subroutine xreal                      (%) =" &
                ,100.*time_elapsed(24)/time_elapsed(1)
  endif

  call MPI_FINALIZE(IERR)
  stop

end subroutine shutdown