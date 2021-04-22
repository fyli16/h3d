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
  use m_parameters
  use functions_mod
  use m_mesh
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
  use m_parameters
  use m_mesh
  implicit none 

  integer*8 :: ixe, iye, ize, i, j, k

  if (myid == 0) then
    open(unit=222, file=trim(restart_directory)//'restart_index.dat', status='old')
    read(222,*) restart_index, itrestart 
    print*, " "
    print*, " Restart from set # ", restart_index
    print*, " Restart from iteration # ", itrestart
    print*, " "
    close(222)
  endif

  call MPI_BCAST(restart_index,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(itrestart    ,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)

  do i = 0, nprocs_over_60 
    if (mod(int(myid,8), nprocs_over_60+1) .eq. i) then 
        call write_read_restart_files(-1.0)  ! read restart data
        call MPI_BCAST(itrestart,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
    endif
  enddo

  ! determine start, finish number of steps
  itstart = itrestart
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
! main simulation loops
!---------------------------------------------------------------------
subroutine sim_loops
  use m_parameters
  use m_diagnostics
  implicit none 

  integer :: i

  ! initialize time arrays, and get a time stamp
  ! just before entering the simulation loop
  time_elapsed=0.; time_begin=0; time_end=0
  call get_time(clock_init)
  clock_old = clock_init
  
  ! main simulation loop
  if (myid==0) print*, "Executing main simulation loops:"

  do while(it <= itfinish)
    call date_and_time(values=time_begin(:,1)) ! time one-whole-loop

    ! print time & step info
    if (myid==0 .and. mod(it,n_print)==0) then
      call get_time(clock_now)
      write(6,"(A5,I7,A2,I7,A11,F8.3,A14,F8.3,A12,F8.3)") 'it = ', it, '/', itfinish, &
                    ',   time = ', time, &
                    ',   delta_t = ', real(clock_now - clock_old), &
                    ',   tot_t = ', real(clock_now - clock_init)
      clock_old = clock_now
    endif

    ! compute resistivity 
    ! (which could depend on local parameters such as current)
    call date_and_time(values=time_begin(:,2))
    if (ndim /= 1) then
      call eta_calc       
    else
      call eta_calc_2d    
    endif
    call date_and_time(values=time_end(:,2))
    call accumulate_time(time_begin(1,2),time_end(1,2),time_elapsed(2))

    ! compute density and v's, and push particles
    call date_and_time(values=time_begin(:,3))
    ntot = 0 ! for particle tracking
    call trans
    call date_and_time(values=time_end(:,3))
    call accumulate_time(time_begin(1,3),time_end(1,3),time_elapsed(3))

    ! sort particles
    call date_and_time(values=time_begin(:,4))
    if (mod(it,n_sort) == 0) call sortit  
    call date_and_time(values=time_end(:,4))
    call accumulate_time(time_begin(1,4),time_end(1,4),time_elapsed(4))

    ! call field solver
    call date_and_time(values=time_begin(:,5))
    if (ndim /=1) then 
      call field
    else
      call field_2d
    endif     
    call date_and_time(values=time_end(:,5))
    call accumulate_time(time_begin(1,5),time_end(1,5),time_elapsed(5))

    ! inject_wave ?
    ! if (it == 21000) call inject_wave
    ! if (mod(it,100) == 0) call kick

    ! call user diagnostics
    call date_and_time(values=time_begin(:,30))
    call diagnostics
    call date_and_time(values=time_end(:,30))
    call accumulate_time(time_begin(1,30),time_end(1,30),time_elapsed(30))

    ! time/step increment
    time = time + dtwci
    it = it + 1

    call date_and_time(values=time_end(:,1))
    call accumulate_time(time_begin(1,1),time_end(1,1),time_elapsed(1))

  enddo 

end subroutine sim_loops


!---------------------------------------------------------------------
! shutdown simulation and exit
!---------------------------------------------------------------------
subroutine shutdown
  use m_parameters
  implicit none 

  if (myid == 0) then
    close(unit=11) ! energy.dat
    close(unit=12) ! probe.dat
    close(unit=13) ! tracking.dat
    ! close(unit=14) ! time.dat
  endif

  if (tracking_mpi) close(unit=13)
      
  if (myid==0) then
    print*, " "
    print*, " "
    print*, " *** Run completed *** "
    print*, " "
    print*, " "
    print*, " total time                   (s)          =",time_elapsed(1)
    print*, "   subroutine eta_cal         (s)          =",time_elapsed(2)
    print*, "   subroutine trans           (s)          =",time_elapsed(3)
    print*, "   subroutine sort            (s)          =",time_elapsed(4)
    print*, "   subroutine field           (s)          =",time_elapsed(5)
    print*, "   subroutine caltemp2        (s)          =",time_elapsed(6)
    print*, "   subroutine diagnostics     (s)          =",time_elapsed(30)
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
    print*, "   total trans                (s)          =",time_elapsed(3)
    print*, " "
    print*, " "
    print*, " In subroutine field,"
    print*, "   subroutine pressgrad       (s)          =",time_elapsed(9)
    print*, "   subroutine bcalc           (s)          =",time_elapsed(10)
    print*, "   subroutine ecalc           (s)          =",time_elapsed(11)
    print*, "   subroutine focalc          (s)          =",time_elapsed(12)
    print*, "   total field                (s)          =",time_elapsed(5)
    print*, " "
    print*, " "
    print*, " In subroutine parmov,"
    print*, "   particle pushing           (s)          =",time_elapsed(13)
    print*, "   particle exchange          (s)          =",time_elapsed(14)
    print*, "   moment calculation         (s)          =",time_elapsed(15)
    print*, "   total parmov               (s)          =",time_elapsed(7)
    print*, " "
    print*, " "
    print*, " In subroutine bcalc,"
    print*, "   subroutine ecalc           (s)          =",time_elapsed(16)
    print*, "   subroutine xrealbcc        (s)          =",time_elapsed(17)
    print*, "   total bcalc                (s)          =",time_elapsed(22)
    print*, " "
    print*, " "
    print*, " In subroutine diagnostics,"
    print*, "   sub virtual_probes         (s)          =",time_elapsed(31)
    print*, "   sub track_particle         (s)          =",time_elapsed(32)
    print*, "   total user_diagnose        (s)          =",time_elapsed(30)
    print*, " "
    print*, " "
    print*, " In subroutine ecalc (called from subroutine bcalc and field),"
    print*, "   subroutine xrealbcc        (s)          =",time_elapsed(18)
    print*, " "
    print*, " "
    print*, "communication time / total time (%)          =", &
                100.*(time_elapsed(24)+time_elapsed(25)+time_elapsed(14) &
                +time_elapsed(17)+time_elapsed(18))/time_elapsed(1)
    print*, " "
    print*, " "
    print*, "Further breakdown of communication time "
    print*, "  particle exchage in subroutine parmov (%) =", &
                100.*time_elapsed(14)/time_elapsed(1)
    print*, "  subroutine xrealbcc                   (%) =", &
                100.*(time_elapsed(25)+time_elapsed(17)+time_elapsed(18))/time_elapsed(1)
    print*, "  subroutine xreal                      (%) =", &
                100.*time_elapsed(24)/time_elapsed(1)
  endif

  call MPI_FINALIZE(IERR)
  stop

end subroutine shutdown