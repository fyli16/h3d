!---------------------------------------------------------------------
!                                 H3D (V6.0)                               
!                           YURI'S NONUNIFORM MESH                         
!                           3D IMPLEMENTATION ONLY                         
!                      UNIFORM LOADING IN PHYSICAL SPACE                   
!               UNIFORM LOADING IN LOGICAL SPACE NOT YET IMPLEMENTED       
!---------------------------------------------------------------------

program h3d 

  use m_parameter
  use m_init

  use m_eta
  use m_particle
  use m_field
  use m_diag

  implicit none

  ! read input deck
  call init_input

  ! MPI domain decomposition 
  call init_decomp

  ! allocate global arrays
  call init_arrays

  ! initialize mesh 
  call init_mesh

  ! restart or a fresh start
  if (restart) then 
    call init_restart 
  else
    call init_wavepart 
  endif 

  ! main simulation loops
  call run_sim

  ! close program and exit
  call close_sim


  contains 

  !---------------------------------------------------------------------
  ! read in input file/parameters
  !---------------------------------------------------------------------
  subroutine init_input
    namelist /datum/ &
    tmax, dtwci, restart, &   ! global info
    MPI_IO_format, &
    nx, ny, nz, xmax, ymax, zmax, npx, npy, npz, node_conf, periods, &  ! simulation domain
    xaa, xbb, nax, nbx, yaa, ybb, nay, nby, zaa, zbb, naz, nbz, &
    uniform_load_logical, &
    n_subcycles, nskipx, nskipy, nskipz, iterb, &  ! field solver
    nspec, n_sort, qspec, wspec, frac, denmin, wpiwci, beta_spec, beta_e, &  ! plasma setup
    ieta, resis, netax, netay, etamin, etamax, eta_par, eta_zs, &
    anisot, gamma, ave1, ave2, phib, smoothing, smooth_pass, &
    dB_B0, num_wave_cycles, &  ! init waves
    n_print, n_diag_mesh, n_diag_energy, n_diag_probe, & ! diagnostics
    n_diag_tracking, n_write_restart, n_diag_particle, &  
    tracking_binary, tracking_mpi, xbox_l, xbox_r, ybox_l, ybox_r, zbox_l, zbox_r, &
    fxsho, nxcel, rcorr, ishape, teti  ! others

    ! Initialize MPI
    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)

    ! convert number 'myid' to character
    ! these characters will be used in dumping restart files by rank
    my_short_int = myid
    call integer_to_character(myid_char, len(myid_char), my_short_int)

    ! print head info
    if (myid==0) then
      print*, " "
      print*, "***************************************************************************"
      print*, "*           |||     |||    ======\\     |||=====\\                        *"
      print*, "*           |||     |||           ||    |||      \\                       *"
      print*, "*           |||=====|||     ======||    |||       ||                      *"
      print*, "*           |||     |||           ||    |||      //                       *"
      print*, "*           |||     |||    ======//     |||=====//                        *"
      print*, "***************************************************************************"
      print*, " "
    endif 

    ! read in input deck
    if (myid == 0) then
      print*, " "
      print*, "Reading input file"
      print*, "-------------------------------------------------"
      open(5, file='input.f90', form='formatted', status='old')
      read(5, nml=datum, iostat=ierr)
      ! write(6, datum)
    endif

    ! Broadcast input parameters (read in at rank 0) to all other ranks
    ! global sim. info
    call MPI_BCAST(tmax                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(dtwci                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(restart                ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(MPI_IO_format          ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    ! sim. domain
    call MPI_BCAST(nx                     ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ny                     ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nz                     ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xmax                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ymax                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zmax                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(npx                    ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(npy                    ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(npz                    ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(node_conf              ,2     ,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(periods                ,2     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xaa                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xbb                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nax                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nbx                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(yaa                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ybb                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nay                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nby                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zaa                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zbb                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(naz                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nbz                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(uniform_load_logical,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
    ! field solver
    call MPI_BCAST(n_subcycles            ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nskipx                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nskipy                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nskipz                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(iterb                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    ! plasma setup
    call MPI_BCAST(nspec                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_sort                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(qspec                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(wspec                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(frac                   ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(denmin                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(wpiwci                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(beta_spec              ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(beta_e                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ieta                   ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(resis                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(netax                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(netay                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(etamin                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(etamax                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(eta_par                ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(eta_zs                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(anisot                 ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(gamma                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ave1                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ave2                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(phib                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(smoothing              ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(smooth_pass            ,1     ,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
    ! init waves
    call MPI_BCAST(dB_B0                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(num_wave_cycles        ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    ! diagnostic control
    call MPI_BCAST(n_print                ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_diag_mesh            ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_diag_energy          ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_diag_probe           ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_diag_tracking        ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_diag_particle        ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_write_restart        ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(tracking_binary        ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(tracking_mpi           ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xbox_l                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xbox_r                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ybox_l                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ybox_r                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zbox_l                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zbox_r                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    ! others 
    call MPI_BCAST(fxsho                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nxcel                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(rcorr                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ishape                 ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(teti                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

    ! The unit of dt is 1/wci in input file, 
    ! but now converted to 1/wpi inside the code
    dt = dtwci * wpiwci

    ! steps of iteration 
    ! (could be modified later in 'init_restart' if restart=.true.)
    it = 0; itrestart = 0; 
    itstart = it; itfinish = tmax/dtwci

    ! field subcycling
    ! n_subcycles = max(n_subcycles, 1_8)

    ! set output directories
    data_directory = 'data/'
    restart_directory = 'restart/'
    restart_index_suffix(1) = '.1'
    restart_index_suffix(2) = '.2'

    ! print some simulation info
    if (myid==0) then
      print*, " "
      if (restart) then 
        print*, "*** Restart run ***"
      else 
        print*, "*** New run ***"
      endif 

      print*, " "
      print*, "--- global simulation info ---"
      print*, "tmax, dtwci = ", tmax, dtwci
      print*, "node_conf   = ", node_conf
      print*
      print*, "--- mesh info ---"
      print*, "xmax, ymax, zmax = ", xmax, ymax, zmax
      print*, "nx,   ny,   nz   = ", nx, ny, nz
      print*, "xaa, xbb = ", xaa, xbb
      print*, "nax, nbx = ", nax, nbx 
      print*, "yaa, ybb = ", yaa, ybb
      print*, "nay, nby = ", nay, nby 
      print*, "zaa, zbb = ", zaa, zbb
      print*, "naz, nbz = ", naz, nbz 
      print*
      print*, "--- plasma info ---"
      print*, "nspec = ", nspec
      print*, "npx   = ", npx
      print*, "npy   = ", npy
      print*, "npz   = ", npz
      
      

    endif  

  end subroutine init_input


  !---------------------------------------------------------------------
  ! initialize MPI domain decomposition
  !---------------------------------------------------------------------
  subroutine init_decomp
    integer :: i

    if (myid==0) then
      print*, " "
      print*, "MPI decompsition"
      print*, "-------------------------------------------------"
      print*, "Total number of processors = ", nprocs
    endif 

    ! specify MPI decomposition (along y, z only; no decomposition along x) 
    if (nz==1 .and. ny==1) then ! only nx>1 and 1 rank will be used  
      ndim=0; dims(1)=1; dims(2)=1
    else if (nz == 1) then ! ny>1 and decomposition only occurs in y
      ndim=1; dims(1)=node_conf(1); dims(2)=1
    else ! ny>1, nz>1, and decomposition in both y and z
      ndim=2; dims(1)=node_conf(1); dims(2)=node_conf(2)
    endif

    ! npy(z) now means number of particles in each rank along y(z)
    npy=npy/dims(1); npz=npz/dims(2)

    ! print MPI decomposition information
    if (myid == 0) then
      do i = 1, ndim
        print*, "Dimension = ", i, " Dims = ", dims(i)
      enddo
    endif

    REORDER = .TRUE. ! reorder ranking or not
    ! Makes a new communicator with topology information attached
    call MPI_CART_CREATE(MPI_COMM_WORLD, NDIM, DIMS, PERIODS, REORDER, COMM2D, IERR) 
    ! Determines calling rank in the new communicator
    call MPI_COMM_RANK(COMM2D, MYID, IERR)
    ! Retrieves Cartesian topology info of the new communicator
    call MPI_CART_GET(COMM2D, NDIM, DIMS, PERIODS, COORDS, IERR)
    ! splits mesh cell elements between nprocs processors
    call MPE_DECOMP1D(NY, DIMS(1), COORDS(1), JB, JE)
    call MPE_DECOMP1D(NZ, DIMS(2), COORDS(2), KB, KE)
    ! print decomposition info (debug purpose)
    ! print*, 'myid=', myid, 'jb, je =', jb, je, 'kb, ke = ',kb, ke, 'coords =', coords

    ! max number of cells. why adding 2?
    nxmax = nx + 2; nymax = ny + 2; nzmax = nz + 2
    ! local max number of cells
    nylmax = je - jb + 1 ; nzlmax = ke - kb + 1  
    if (myid == 0) then
      print*, "Local array size in x-direction = ", nx
      print*, "Local array size in y-direction = ", nylmax
      print*, "Local array size in z-direction = ", nzlmax
    endif

  end subroutine init_decomp

  !---------------------------------------------------------------------
  ! main simulation loops
  !---------------------------------------------------------------------
  subroutine run_sim
    integer :: i

    ! initialize time arrays, and get a time stamp
    ! just before entering the simulation loop
    time_elapsed=0.; time_begin=0; time_end=0
    call get_time(clock_init)
    clock_old = clock_init
    
    ! main simulation loop
    if (myid==0) then 
      print*, " " 
      print*, "Executing main simulation loops:"
      print*, "-------------------------------------------------"
    endif 

    do while(it <= itfinish)

      call date_and_time(values=time_begin(:,1))

      ! print time & step info
      if (myid==0 .and. mod(it,n_print)==0) then
        call get_time(clock_now)
        write(6,"(A5,I7,A2,I7,A11,F8.3,A14,F8.3,A12,F8.3)") 'it = ', it, '/', itfinish, &
                      ',   time = ', time, &
                      ',   delta_t = ', real(clock_now - clock_old), &
                      ',   tot_t = ', real(clock_now - clock_init)
        clock_old = clock_now
      endif

      ! calculate resistivity (eta)
      ! (which could depend on local parameters such as current)
      call date_and_time(values=time_begin(:,2))
      if (ndim /= 1) then
        call cal_eta       
      else
        call cal_eta_2d    
      endif
      call date_and_time(values=time_end(:,2))
      call add_time(time_begin(1,2),time_end(1,2),time_elapsed(2))

      ! calculate density and v's, and push particles
      call date_and_time(values=time_begin(:,3))
      ntot = 0 ! for particle tracking
      call trans
      call date_and_time(values=time_end(:,3))
      call add_time(time_begin(1,3),time_end(1,3),time_elapsed(3))

      ! sort particles
      call date_and_time(values=time_begin(:,4))
      if (mod(it,n_sort) == 0) call sortit  
      call date_and_time(values=time_end(:,4))
      call add_time(time_begin(1,4),time_end(1,4),time_elapsed(4))

      ! call field solver
      call date_and_time(values=time_begin(:,5))
      if (ndim /=1) then 
        call field
      else
        call field_2d
      endif     
      call date_and_time(values=time_end(:,5))
      call add_time(time_begin(1,5),time_end(1,5),time_elapsed(5))

      ! call user diagnostics
      call date_and_time(values=time_begin(:,6))
      call diag
      call date_and_time(values=time_end(:,6))
      call add_time(time_begin(1,6),time_end(1,6),time_elapsed(6))

      ! time/step increment
      time = time + dtwci
      it = it + 1

      call date_and_time(values=time_end(:,1))
      call add_time(time_begin(1,1),time_end(1,1),time_elapsed(1))

    enddo 

  end subroutine run_sim


  !---------------------------------------------------------------------
  ! close simulation and exit
  !---------------------------------------------------------------------
  subroutine close_sim

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
      print*, "*** Run completed *** "
      print*, " "
      print*, " "
      print*, "total time                        (s)          =",time_elapsed(1)
      print*, "   sub cal_eta                    (s)          =",time_elapsed(2)
      print*, "   sub trans                      (s)          =",time_elapsed(3)
      print*, "   sub sort                       (s)          =",time_elapsed(4)
      print*, "   sub field                      (s)          =",time_elapsed(5)
      print*, "   sub diag                       (s)          =",time_elapsed(6)
      print*, " "
      print*, " "
      print*, "In trans.parmov,"
      print*, "   sub push                       (s)          =",time_elapsed(33)
      print*, "   sub particle_boundary          (s)          =",time_elapsed(34)
      print*, "   sub moment_calculation         (s)          =",time_elapsed(35)
      print*, "   total parmov                   (s)          =",time_elapsed(31)
      print*, " "
      print*, " "
      print*, "In field,"
      print*, "   sub pressgrad                  (s)          =",time_elapsed(51)
      print*, "   sub bcalc                      (s)          =",time_elapsed(52)
      print*, "   sub ecalc                      (s)          =",time_elapsed(53)
      print*, "   sub focalc                     (s)          =",time_elapsed(54)
      print*, "   total field                    (s)          =",time_elapsed(5)
      print*, " "
      print*, " "
      print*, "In diag,"
      print*, "   sub diag_mesh                  (s)          =",time_elapsed(61)
      print*, "   sub diag_energy                (s)          =",time_elapsed(62)
      print*, "   sub diag_particle              (s)          =",time_elapsed(63)
      print*, "   sub diag_probe                 (s)          =",time_elapsed(64)
      print*, "   sub diag_tracking              (s)          =",time_elapsed(65)
      print*, "   sub write_restart              (s)          =",time_elapsed(66)
      print*, "   total diag                     (s)          =",time_elapsed(6)
    endif

    call MPI_FINALIZE(IERR)
    stop

  end subroutine close_sim

end program h3d