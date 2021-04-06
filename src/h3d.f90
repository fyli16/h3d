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
    use functions_f90
    use mesh2d
    implicit none

    integer, dimension(8) :: curr_time
    integer*8 :: itstart, itfinish
    real*8, dimension(:,:,:), allocatable :: uniform_mesh      
    ! VR : allocating a global mesh can not work on large runs with small amount of memory per rank
    ! double precision, dimension(:,:,:), allocatable:: nonuniform_mesh_global

    ! Initialize MPI
    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
    if (myid==0) then
      write(6,*) 
      write(6,*) "Init mpi ..."
      write(6,*) "  Total number of processors = ", numprocs
    endif

    ! Read input file
    call read_input()
    
    ! Decompose MPI following input  
    call mpi_decomposition()
    
    ! Set global parameters
    call set_parameters()  
      
    ! set up mesh 
    call setup_mesh()
    allocate ( uniform_mesh(nxmax, jb-1:je+1, kb-1:ke+1) )
    ! VR: allocate (nonuniform_mesh_global(nxmax,0:ny+1,0:nz+1))
 
    ! open history diagnostic files
    call open_hist_diag_files()

    ! notime=0 will output detailed timing
    if (notime==0) call open_timing_diag_files

    ! make list of particles (see utils.f90)
    call makelist

    ! restart or a fresh start
    if (restart) then
      call init_restart()
    else  ! fresh start
      time = 0.0
      restart_index = 1
      call init_wave
    endif

    ! time stamp just before entering the simulation loop
    call date_and_time(values=curr_time)
    clock_time_init=( curr_time(5)*3600.+curr_time(6)*60.+curr_time(7)+curr_time(8)*0.001)
    clock_time_old = clock_time_init

    ! march forward in time
    itstart = itfin + 1
    it = itstart

    ! change how itfinish is computed
    itfinish = (tmax-t_stopped)/dtwci + itstart - 1
    if (myid == 0) write(6,*) 't_stopped = ', t_stopped
    if (myid == 0) write(6,*) 'itstart, itfinish = ', itstart, ' ', itfinish
    
    time_elapsed=0.; time_begin_array=0; time_end_array=0

    ! main simulation loop
    do while(it <= itfinish)
      call one_simulation_loop(itstart, itfinish, uniform_mesh)
    enddo  

    ! close files that are open
    call close_files()

    ! print misc timing info, finalize MPI, and shutdown the program
    call shutdown()

  end program h3d


!---------------------------------------------------------------------
subroutine setup_mesh()
  use parameter_mod
  use mesh2d
  implicit none
  
  integer :: i

  if (myid==0) then
    write(6,*)
    write(6,*) "Setting up mesh ..."
  endif

  ! Initialize nonuniform mesh
  call MESH_INIT(meshX,xaa,xbb,xmax,nax,nbx,nx) ! initialize x-mesh
  call MESH_INIT(meshY,yaa,ybb,ymax,nay,nby,ny) ! initialize y-mesh
  call MESH_INIT(meshZ,zaa,zbb,zmax,naz,nbz,nz) ! initialize z-mesh

  call MESH_INDEX(meshX,CELL,ixv_2_c_map)
  call MESH_INDEX(meshY,CELL,iyv_2_c_map)
  call MESH_INDEX(meshZ,CELL,izv_2_c_map)
  call MESH_INDEX(meshX,NODE,ixv_2_v_map)
  call MESH_INDEX(meshY,NODE,iyv_2_v_map)
  call MESH_INDEX(meshZ,NODE,izv_2_v_map)
  call MESH_INDEX(meshX,CELL,ixc_2_c_map,CELL)
  call MESH_INDEX(meshY,CELL,iyc_2_c_map,CELL)
  call MESH_INDEX(meshZ,CELL,izc_2_c_map,CELL)
  call MESH_INDEX(meshX,NODE,ixc_2_v_map,CELL)
  call MESH_INDEX(meshY,NODE,iyc_2_v_map,CELL)
  call MESH_INDEX(meshZ,NODE,izc_2_v_map,CELL)

  ! write mesh properties into a file
  if (myid == 0) then
    open(unit=10,file='mesh_vertices.dat',status='unknown',form='formatted')
    write(10,*) meshX%nl+1, meshY%nl+1, meshZ%nl+1

    do i=2,meshX%nl+2 
      write(10,*) meshX%xn(i)
    enddo
    do i=2,meshX%nl+2 
      write(10,*) meshX%dxc(i)
    enddo

    do i=2,meshY%nl+2 
      write(10,*) meshY%xn(i)
    enddo
    do i=2,meshY%nl+2 
      write(10,*) meshY%dxc(i)
    enddo

    do i=2,meshZ%nl+2 
      write(10,*) meshZ%xn(i)
    enddo
    do i=2,meshZ%nl+2 
      write(10,*) meshZ%dxc(i)
    enddo
    
    close(unit=10)
  endif

  return
end subroutine setup_mesh


!---------------------------------------------------------------------
! initialize data for a restart run
subroutine init_restart()
  use parameter_mod
  use mesh2d
  implicit none 

  integer*8 :: ixe, iye, ize, i, j, k, itstart
  integer :: is, iwrite

  if (myid == 0) then
    open(unit=222,file=trim(adjustl(restart_directory))//'restart_index.dat' ,status='old')
    read(222,*) restart_index, itfin
    close(222)
  endif
  call MPI_BCAST(restart_index,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(itfin        ,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)

  if (myid == 0) then
    write(6,*) " "
    write(6,*) " Restart from set # ", restart_index
    write(6,*) " Restart at iteration # ", itfin
    write(6,*) " "
  endif

  do iwrite = 0, npes_over_60 
    if (mod( int(myid,8), npes_over_60 + 1) .eq. iwrite) then
        call restart_read_write(-1.0, itstart)  ! read restart data
        call MPI_BCAST(itfin,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
    endif
  enddo
    
  if (restart_index == 1) then
    restart_index=2
  else
    restart_index=1
  endif

  ! Uniform mesh - Same as is in version 5.0
  if (myid==0) then
    write(6,*) " "
    write(6,*) "hx, hy, hz =", hx, hy, hz
    write(6,*) " "
  endif
  yb = (jb-1)*hy  ! where is hy, hz defined before this?
  ye = je    *hy
  zb = (kb-1)*hz
  ze = ke    *hz
    
  ! Nonuniform mesh
  zb = meshZ%xn(kb+1)
  ze = meshZ%xn(ke+2)
  do ipe = 0,npes-1
    zbglobal(ipe) = meshZ%xn(kbglobal(ipe)+1)
    zeglobal(ipe) = meshZ%xn(keglobal(ipe)+2)
  enddo

  yb = meshY%xn(jb+1)
  ye = meshY%xn(je+2)
  do ipe = 0,npes-1
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
            
  do is = 1, nspec
    npm = npx(is)*npy(is)*npz(is)*npes
    dfac(is) = real(ny*nz*nx)/real(npm)
    do ixe=1,nx2
        do iye=jb-1,je+1
          do ize=kb-1,ke+1
              qp_cell(ixe,iye,ize,is) = meshX%dxc(ixe)*meshY%dxc(iye+1)*meshZ%dxc(ize+1)*dfac(is)*frac(is)
          enddo
        enddo
    enddo
  enddo
          
  do i=1,nxmax
    xc_uniform(i) = hx*(i-1.5)
    xv_uniform(i) = hx*(i-2.0)
  enddo
  do j = 1, nymax
    yc_uniform(j) = hy*(j-0.5)
    yv_uniform(j) = hy*(j-1.0)
  enddo
  do k=1,nzmax
    zc_uniform(k) = hz*(k-0.5)
    zv_uniform(k) = hz*(k-1.0)
  enddo
    
  if (myid ==0) then
      my_short_int = it
      call integer_to_character(cycle_ascii,len(cycle_ascii),my_short_int)
      cycle_ascii_new=trim(adjustl(cycle_ascii))
      write(6,*) " cycle = ", cycle_ascii_new
  endif
  
  call MPI_BCAST(cycle_ascii    ,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(cycle_ascii_new,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

  return 

end subroutine init_restart


!---------------------------------------------------------------------
subroutine data_output(uniform_mesh)
  use parameter_mod
  implicit none 

  integer :: j
  integer*8 :: numvars, irecnum
  real*8 :: uniform_mesh(nxmax,jb-1:je+1,kb-1:ke+1)

  if (myid ==0) then
    my_short_int=it
    call integer_to_character(cycle_ascii, len(cycle_ascii), my_short_int)
    cycle_ascii_new=trim(adjustl(cycle_ascii))
    write(6,*) " cycle = ",cycle_ascii_new
  endif
  call MPI_BCAST(cycle_ascii,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(cycle_ascii_new,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  
  if (myid == 0 .and. .not. MPI_IO_format) call open_files
  
  call date_and_time(values=time_begin_array(:,6))
  if (ndim /= 1) then
    call caltemp2_global
  else
    call caltemp2_global_2d
  endif
  call date_and_time(values=time_end_array(:,6))
  call accumulate_time_difference(time_begin_array(1,6),time_end_array(1,6),time_elapsed(6))

  numvars = 18
  irecnum = 1
  ! VR: in the old version, "nonuniform_mesh_global" was passed to dataout
  call dataout(bx, by, bz, den, ex, ey, ez, vix, viy, viz, tpar, tperp,                   &
      p_xx,p_xy,p_xz,p_yy,p_yz,p_zz,fox,foy,foz, vxs,vys,vzs,                 &
      nxmax,nymax,nzmax,file_unit,myid,                                       &
      numvars,irecnum,kb,ke,numprocs,wpiwci,jb,je,ny,nz,nylmax,nzlmax,nspecm, &
      eta, eta_times_b_dot_j, eta_par,            &
      uniform_mesh, trim(data_directory), trim(adjustl(cycle_ascii)), MPI_IO_format)
  if (myid == 0 .and. .not. MPI_IO_format) then
    do j = 1, 25
      close(file_unit(j))
    enddo
  endif

  if (t_begin*wpiwci <= time .and. time <= t_end*wpiwci .and. mod(int(it,8),n_write_particle)==0 .or. it==1) then
    if (myid == 0) then
      my_short_int = it
      call integer_to_character(cycle_ascii, len(cycle_ascii), my_short_int)
      write(6,*) " calling particle_in_volume_write with cycle_ascii = ", cycle_ascii
    endif
    call MPI_BCAST(cycle_ascii,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
    call particle_in_volume_write
  endif

  return 

end subroutine data_output


!---------------------------------------------------------------------
subroutine one_simulation_loop(itstart, itfinish, uniform_mesh)
  use parameter_mod
  implicit none 

  integer, dimension(8) :: curr_time
  integer :: iwrite, i
  integer*8 :: itstart, itfinish
  real*8 :: uniform_mesh(nxmax,jb-1:je+1,kb-1:ke+1)

  call get_cleanup_status(len(cleanup_status))

  if (cleanup_status=='CLEANUP_STATUS=EXIT') then
    call shutdown()
  endif

  call date_and_time(values=time_begin_array(:,1))
  call date_and_time(values=curr_time)
  clock_time=( curr_time(5)*3600.+curr_time(6)*60.+curr_time(7)+curr_time(8)*0.001)
  if (notime == 0) then
    write(file_unit_time,"(i4,' begin    ',f15.3)") it, real(clock_time-clock_time_init)
  endif
  if (myid == 0.and.mod(it,10_8) == 0) then
    write(6,*) 'it=', it, ', delta_time=', real(clock_time - clock_time_old), ', tot_time=', real(clock_time-clock_time_init)
    clock_time_old = clock_time
  endif

  ! determine whether or not to print diagnostic information
  if (mod(it,n_print) == 0) then
    print_info=.true.
  else
    print_info=.false.
  endif

  ! determine whether or not to write data into files
  if (mod(it,n_write_data) == 0 ) then
    write_data=.true.
  else
    write_data=.false.
  endif

  call date_and_time(values=time_begin_array(:,3))

  ! VR: we do not need injection from the boundary 
  ! do is=1,nspec
  !   call particle_newinject_linked_list(is)
  ! enddo

  call date_and_time(values=time_end_array(:,3))
  call date_and_time(values=time_begin_array(:,2))

  !VR: compute values of resistivity (that could depend on local parameters, such as
  !VR: current
  if (ndim /= 1) then
    call etacalc       ! Dietmar's resistivity
  else
    call etacalc_2d    ! Dietmar's resistivity
  endif
     
  if (notime == 0) then
    call date_and_time(values=curr_time)
    clock_time=( curr_time(5)*3600.+curr_time(6)*60.+curr_time(7)+curr_time(8)*0.001)
    write(file_unit_time,"(i4,' trans    ',f15.3)") it, real(clock_time-clock_time_init)
  endif

  ! VR: trans computes density and v's 
  ! VR: note that trans calls parmov, i.e. it also does a particle push
  ntot=0 ! for particle tracking
  call trans
  call date_and_time(values=time_end_array(:,2))

  ! write output to energy.dat
  if (myid==0) then
    write(11,*) it, efld, bfld, efluidt, ethermt, eptclt ! energy.dat
    write(14,*) it, time_elapsed(1:40)  ! time.dat
  endif

  if (notime == 0) then
    call date_and_time(values=curr_time)
    clock_time=( curr_time(5)*3600.+curr_time(6)*60.+curr_time(7)+curr_time(8)*0.001)
    write(file_unit_time,"(i4,' sortit   ',f15.3)") it,real(clock_time-clock_time_init)
  endif

  ! VR: sort particles
  call date_and_time(values=time_begin_array(:,4))
  if (mod(it,10_8) == 0) call sortit    !  sort the particles
  call date_and_time(values=time_end_array(:,4))

  if (notime == 0) then
    call date_and_time(values=curr_time)
    clock_time=( curr_time(5)*3600.+curr_time(6)*60.+curr_time(7)+curr_time(8)*0.001)
    write(file_unit_time,"(i4,' field    ',f15.3)") it,real(clock_time-clock_time_init)
  endif

  ! VR: call field solver
  call date_and_time(values=time_begin_array(:,5))
  if (.not.testorbt) then
    if (ndim /=1) then 
      call field
    else
      call field_2d
    endif
  endif        
  call date_and_time(values=time_end_array(:,5))

  ! ??
  ! if (it == 21000) call inject_wave
  ! if (mod(it,100) == 0) call kick
  if (notime == 0) then
    call date_and_time(values=curr_time)
    clock_time=( curr_time(5)*3600.+curr_time(6)*60.+curr_time(7)+curr_time(8)*0.001)
    write(file_unit_time,"(i4,' diagnose ',f15.3)") it, real(clock_time-clock_time_init)
  endif

  ! VR: call user diagnostics
  call date_and_time(values=time_begin_array(:,30))
  call user_diagnostics
  call date_and_time(values=time_end_array(:,30))
  call accumulate_time_difference(time_begin_array(1,30),time_end_array(1,30),time_elapsed(30))
    
  ! VR: this seems to be data output region
  if (.not.testorbt.and.(write_data.or.it==itfinish.or.it==1)) then        
    call data_output(uniform_mesh)
  endif

  ! write restart files at specified steps
  if (write_restart .and. mod(int(it,8), n_write_restart)==0) then
    if (myid == 0) then
      WRITE(6,*)
      WRITE(6,*) 'Writing restart files ...'
      WRITE(6,*)
    endif

    itfin = it

    do iwrite = 0, npes_over_60  
      if (mod( int(myid,8) ,npes_over_60 + 1).eq.iwrite) then
        call restart_read_write(1.0, itstart)
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    enddo

    if (myid == 0) then
    open(unit=222,file=trim(restart_directory)//'restart_index.dat' ,status='unknown')
      write(222,*) restart_index, itfin
      close(222)
    endif

    if (restart_index == 1) then
      restart_index=2
    else
      restart_index=1
    endif

  endif

  if (notime == 0) then
    call date_and_time(values=curr_time)
    clock_time=( curr_time(5)*3600.+curr_time(6)*60.+curr_time(7)+curr_time(8)*0.001)
    write(file_unit_time,"(i4,' end      ',f15.3)") it,real(clock_time-clock_time_init)
  endif

  time = time + dt
  it = it + 1

  call date_and_time(values=time_end_array(:,1))

  do i = 1, 5
    call accumulate_time_difference(time_begin_array(1,i),time_end_array(1,i),time_elapsed(i))
  enddo

  return

end subroutine one_simulation_loop


!---------------------------------------------------------------------
subroutine close_files()
  use parameter_mod, only : myid, tracking_mpi
  implicit none 

  if (myid == 0) then
    close(unit=11) ! energy.dat
    close(unit=12) ! probe.dat
    close(unit=13) ! tracking.dat
    close(unit=14) ! time.dat
  endif

  if (tracking_mpi) then
    close(unit=13)
  endif

  return 
end subroutine close_files


!---------------------------------------------------------------------
! shutdown the simulation and exit
subroutine shutdown()
  use parameter_mod
  implicit none 

  if (notime == 0) close(file_unit_time)
      
  if (myid==0) then
    write(6,*) " "
    write(6,*) " "
    write(6,*) " *** RUN COMPLETED *** RUN COMPLETED *** RUN COMPLETED "
    write(6,*) " "
    write(6,*) " "
    write(6,*) " subroutine trans             (s)          =",time_elapsed(2)
    write(6,*) " subroutine injectpar         (s)          =",time_elapsed(3)
    write(6,*) " subroutine sort              (s)          =",time_elapsed(4)
    write(6,*) " subroutine field             (s)          =",time_elapsed(5)
    write(6,*) " subroutine caltemp2          (s)          =",time_elapsed(6)
    write(6,*) " subroutine user_diagnostics  (s)          =",time_elapsed(30)
    write(6,*) " total time                   (s)          =",time_elapsed(1)
    write(6,*) " "
    write(6,*) " "
    write(6,*) " In subroutine caltemp2,"
    write(6,*) "   subroutine xreal           (s)          =",time_elapsed(24)
    write(6,*) "   subroutine xrealbcc        (s)          =",time_elapsed(25)
    write(6,*) "   interpolation              (s)          =",time_elapsed(26)
    write(6,*) "   total caltemp2             (s)          =",time_elapsed(23)
    write(6,*) " "
    write(6,*) " "
    write(6,*) " In subroutine trans," 
    write(6,*) "   subroutine parmov          (s)          =",time_elapsed(7)
    write(6,*) "   subroutine energy          (s)          =",time_elapsed(8)
    write(6,*) "   total trans                (s)          =",time_elapsed(20)
    write(6,*) " "
    write(6,*) " "
    write(6,*) " In subroutine field,"
    write(6,*) "   subroutine pressgrad       (s)          =",time_elapsed(9)
    write(6,*) "   subroutine bcalc           (s)          =",time_elapsed(10)
    write(6,*) "   subroutine ecalc           (s)          =",time_elapsed(11)
    write(6,*) "   subroutine focalc          (s)          =",time_elapsed(12)
    write(6,*) "   total field                (s)          =",time_elapsed(21)
    write(6,*) " "
    write(6,*) " "
    write(6,*) "   In subroutine parmov,"
    write(6,*) "     particle pushing         (s)          =",time_elapsed(13)
    write(6,*) "     particle exchange        (s)          =",time_elapsed(14)
    write(6,*) "     moment calculation       (s)          =",time_elapsed(15)
    write(6,*) "     total parmov             (s)          =",time_elapsed(19)
    write(6,*) " "
    write(6,*) " "
    write(6,*) "   In subroutine bcalc,"
    write(6,*) "     subroutine ecalc         (s)          =",time_elapsed(16)
    write(6,*) "     subroutine xrealbcc      (s)          =",time_elapsed(17)
    write(6,*) "     total bcalc              (s)          =",time_elapsed(22)
    write(6,*) " "
    write(6,*) " "
    write(6,*) "   In subroutine user_diagnostics,"
    write(6,*) "     sub virtual_probes       (s)          =",time_elapsed(31)
    write(6,*) "     sub track_particle       (s)          =",time_elapsed(32)
    write(6,*) "     total user_diagnose      (s)          =",time_elapsed(30)
    write(6,*) " "
    write(6,*) " "
    write(6,*) "     In subroutine ecalc (called from subroutine bcalc and field),"
    write(6,*) "       subroutine xrealbcc    (s)          =",time_elapsed(18)
    write(6,*) " "
    write(6,*) " "
    write(6,*) "communication time/total time (%)          =" &
    ,100.*(time_elapsed(24)+time_elapsed(25)+time_elapsed(14)+time_elapsed(17) &
    +time_elapsed(18))/time_elapsed(1)
    write(6,*) " "
    write(6,*) " "
    write(6,*) "Further breakdown of communication time "
    write(6,*) "  particle exchage in subroutine parmov (%) =" &
    ,100.*time_elapsed(14)/time_elapsed(1)
    write(6,*) "  subroutine xrealbcc                   (%) =" &
    ,100.*(time_elapsed(25)+time_elapsed(17)+time_elapsed(18))/time_elapsed(1)
    write(6,*) "  subroutine xreal                      (%) =" &
    ,100.*time_elapsed(24)/time_elapsed(1)
  endif

  call MPI_FINALIZE(IERR)

  stop

end subroutine shutdown