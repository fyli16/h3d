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

  real*8, dimension(:,:,:), allocatable :: uniform_mesh     
  ! VR : allocating a global mesh can not work on large runs with small amount of memory per rank 
  ! real*8, dimension(:,:,:), allocatable:: nonuniform_mesh_global

  ! Initialize MPI
  call MPI_INIT(IERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)

  ! Read input file
  call read_input
  ! Decompose MPI/simulation domain  
  call domain_decomposition
  ! allocate global parameters
  call allocate_global_arrays  
    
  ! set up mesh 
  call setup_mesh()
  allocate ( uniform_mesh(nxmax, jb-1:je+1, kb-1:ke+1) )
  ! VR: allocate (nonuniform_mesh_global(nxmax,0:ny+1,0:nz+1))

  ! open history diagnostic files
  call open_hist_diag_files()

  ! restart or a fresh start
  call makelist  ! make list of particles (see utils.f90)
  if (restart) then
    call init_restart()
  else  ! fresh start
    call init_wave
  endif

  ! time stamp just before entering the simulation loop
  call date_and_time(values=curr_time)
  clock_time_init = curr_time(5)*3600.+curr_time(6)*60.+curr_time(7)+curr_time(8)*0.001
  clock_time_old = clock_time_init

  ! determine start, finish number of steps
  itstart = itfin
  it = itstart
  itfinish = (tmax-t_stopped)/dtwci + itstart
  if (myid==0) then 
    print*, " "
    print*, 't_stopped = ', t_stopped
    print*, 'itstart = ', itstart
    print*, 'itfinish = ', itfinish
    print*, " "
  endif 
  
  ! main simulation loop
  time_elapsed=0.; time_begin_array=0; time_end_array=0
  do while(it <= itfinish)
    call one_simulation_loop(uniform_mesh)
  enddo  

  ! shutdown the program
  call shutdown()

end program h3d


!---------------------------------------------------------------------
! set up uniform mesh
!---------------------------------------------------------------------
subroutine setup_mesh()
  use parameter_mod
  use mesh2d
  implicit none
  
  integer :: i

  if (myid==0) then
    write(6,*) " "
    write(6,*) "Setting up mesh ..."
  endif

  ! Initialize nonuniform mesh
  ! where meshX, meshY, meshZ are declared in 'mesh2d'
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

    do i = 2, meshX%nl+2 
      write(10,*) meshX%xn(i)
    enddo
    do i = 2, meshX%nl+2 
      write(10,*) meshX%dxc(i)
    enddo

    do i = 2, meshY%nl+2 
      write(10,*) meshY%xn(i)
    enddo
    do i = 2, meshY%nl+2 
      write(10,*) meshY%dxc(i)
    enddo

    do i = 2, meshZ%nl+2 
      write(10,*) meshZ%xn(i)
    enddo
    do i = 2, meshZ%nl+2 
      write(10,*) meshZ%dxc(i)
    enddo
    
    close(unit=10)
  endif

  return
end subroutine setup_mesh


!---------------------------------------------------------------------
! initialize data for a restart run
!---------------------------------------------------------------------
subroutine init_restart()
  use parameter_mod
  use mesh2d
  implicit none 

  integer*8 :: ixe, iye, ize, i, j, k
  integer :: is, iwrite

  if (myid == 0) then
    open(unit=222, file=trim(restart_directory)//'restart_index.dat', status='old')
    read(222,*) restart_index, itfin ! itfin is the final restart dump of previous run
    write(6,*) " "
    write(6,*) " Restart from set # ", restart_index
    write(6,*) " Restart from iteration # ", itfin
    write(6,*) " "
    close(222)
  endif

  call MPI_BCAST(restart_index,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(itfin        ,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)

  do iwrite = 0, npes_over_60 
    if (mod(int(myid,8), npes_over_60+1) .eq. iwrite) then 
        call restart_read_write(-1.0)  ! read restart data
        call MPI_BCAST(itfin,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
    endif
  enddo
    
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
    write(6,*) " cycle = ", cycle_ascii_new
  endif
  
  call MPI_BCAST(cycle_ascii    ,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(cycle_ascii_new,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

  return 

end subroutine init_restart


!---------------------------------------------------------------------
! diagnostic data output
!---------------------------------------------------------------------
subroutine data_output(uniform_mesh)
  use parameter_mod
  implicit none 

  integer :: j
  integer*8 :: numvars, irecnum
  real*8 :: uniform_mesh(nxmax, jb-1:je+1, kb-1:ke+1)

  ! ??
  if (myid==0) then
    my_short_int = it
    call integer_to_character(cycle_ascii, len(cycle_ascii), my_short_int)
    if (cycle_ascii=='') cycle_ascii='0'
    cycle_ascii_new=trim(adjustl(cycle_ascii))
    write(6,*) " cycle = ", cycle_ascii_new
  endif
  call MPI_BCAST(cycle_ascii,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(cycle_ascii_new,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
  
  if (myid==0 .and. .not.MPI_IO_format) call open_files
  
  call date_and_time(values=time_begin_array(:,6))
  if (ndim /= 1) then
    call caltemp2_global
  else
    call caltemp2_global_2d
  endif
  call date_and_time(values=time_end_array(:,6))
  call accumulate_time(time_begin_array(1,6),time_end_array(1,6),time_elapsed(6))

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

  if (mod(it,n_write_particle)==0) then
    if (myid == 0) then
      my_short_int = it
      call integer_to_character(cycle_ascii, len(cycle_ascii), my_short_int)
      if (cycle_ascii=='') cycle_ascii='0'
      write(6,*) " calling particle_in_volume_write at cycle = ", cycle_ascii
    endif
    call MPI_BCAST(cycle_ascii,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
    call particle_in_volume_write
  endif

  return 

end subroutine data_output


!---------------------------------------------------------------------
! simulation loop
!---------------------------------------------------------------------
subroutine one_simulation_loop(uniform_mesh)
  use parameter_mod
  implicit none 

  integer :: iwrite
  real*8 :: uniform_mesh(nxmax,jb-1:je+1,kb-1:ke+1)

  call date_and_time(values=time_begin_array(:,1)) ! time one-whole-loop

  ! print time-step info
  call date_and_time(values=curr_time)
  clock_time = curr_time(5)*3600.+curr_time(6)*60.+curr_time(7)+curr_time(8)*0.001
  if (myid == 0.and.mod(it,10_8) == 0) then
    write(6,*) 'it=', it, ', time=', time, ', delta_t=', real(clock_time-clock_time_old), &
                ', tot_t=', real(clock_time-clock_time_init)
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
  call trans ! trans computes density and v's; it also calls parmov, i.e., it does particle push
  call date_and_time(values=time_end_array(:,2))
  call accumulate_time(time_begin_array(1,2),time_end_array(1,2),time_elapsed(2))

  ! write output to energy.dat
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
    
  ! write data
  if (.not.testorbt.and.mod(it,n_write_data)==0) then        
    call data_output(uniform_mesh)
  endif

  ! write restart files
  if (mod(int(it,8),n_write_restart)==0 .and. (it>itstart)) then
    if (myid == 0) write(6,*) 'Writing restart files ...'
    itfin = it
    do iwrite = 0, npes_over_60  
      if (mod( int(myid,8) ,npes_over_60 + 1).eq.iwrite) then
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

    ! keep only two sets of restart files, 
    ! so overwrite the previous dump by swapping file index 
    if (restart_index == 1) then
      restart_index=2
    else
      restart_index=1
    endif

  endif

  time = time + dtwci
  it = it + 1

  call date_and_time(values=time_end_array(:,1))
  call accumulate_time(time_begin_array(1,1),time_end_array(1,1),time_elapsed(1))

  return

end subroutine one_simulation_loop


!---------------------------------------------------------------------
! shutdown the simulation and exit
!---------------------------------------------------------------------
subroutine shutdown()
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
    write(6,*) " "
    write(6,*) " "
    write(6,*) " *** RUN COMPLETED *** "
    write(6,*) " "
    write(6,*) " "
    write(6,*) " subroutine trans             (s)          =",time_elapsed(2)
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
    write(6,*) " In subroutine parmov,"
    write(6,*) "   particle pushing           (s)          =",time_elapsed(13)
    write(6,*) "   particle exchange          (s)          =",time_elapsed(14)
    write(6,*) "   moment calculation         (s)          =",time_elapsed(15)
    write(6,*) "   total parmov               (s)          =",time_elapsed(19)
    write(6,*) " "
    write(6,*) " "
    write(6,*) " In subroutine bcalc,"
    write(6,*) "   subroutine ecalc           (s)          =",time_elapsed(16)
    write(6,*) "   subroutine xrealbcc        (s)          =",time_elapsed(17)
    write(6,*) "   total bcalc                (s)          =",time_elapsed(22)
    write(6,*) " "
    write(6,*) " "
    write(6,*) " In subroutine user_diagnostics,"
    write(6,*) "   sub virtual_probes         (s)          =",time_elapsed(31)
    write(6,*) "   sub track_particle         (s)          =",time_elapsed(32)
    write(6,*) "   total user_diagnose        (s)          =",time_elapsed(30)
    write(6,*) " "
    write(6,*) " "
    write(6,*) " In subroutine ecalc (called from subroutine bcalc and field),"
    write(6,*) "   subroutine xrealbcc        (s)          =",time_elapsed(18)
    write(6,*) " "
    write(6,*) " "
    write(6,*) "communication time/total time (%)          =" &
                ,100.*(time_elapsed(24)+time_elapsed(25)+time_elapsed(14) &
                +time_elapsed(17)+time_elapsed(18))/time_elapsed(1)
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