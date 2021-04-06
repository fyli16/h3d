module initialize
  use parameter_mod
  implicit none

  integer*4:: input_error

  contains
  !---------------------------------------------------------------------
  ! MPI initialization
  subroutine init_mpi()  
    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
    if (myid==0) then
      write(6,*) 
      write(6,*) "Init mpi ..."
      write(6,*) "  Total number of processors = ", numprocs
    endif

    ! get_simulation_id
    ! if (myid==0) call get_sim_id(sim_id)

    ! time stamp
    !  call date_and_time(values=wall_clock_begin)
    initial_time = MPI_Wtime()
    return
  end subroutine init_mpi


  !---------------------------------------------------------------------
  subroutine read_input()
    external get_environment_variable1, get_environment_variable2
    ! external get_environment_variable

    namelist /datum/ tmax, t_begin, t_end, dtwci, dt, restart, &   ! global info
    restrt_write, quota, MPI_IO_format, &
    nx, ny, nz, xmax, ymax, zmax, npx, npy, npz, &  ! simulation domain
    nodey, nodez, &
    xaa, xbb, nax, nbx, yaa, ybb, nay, nby, zaa, zbb, naz, nbz, &
    uniform_loading_in_logical_grid, &
    buffer_zone, moat_zone, profile_power, &
    n_subcycles, nskipx, nskipy, nskipz, iterb, testorbt, norbskip, &  ! field solver
    nspec, qspec, wspec, frac, denmin, wpiwci, btspec, bete, &  ! plasma setup
    ieta, resis, netax, netay, netaz, etamin, etamax, eta_par, &
    anisot, gama, ave1, ave2, phib, smoothing, smooth_coef, &
    dB_B0, num_cycles, &  ! init waves
    nprint, nwrtdata, nwrtrestart, nwrtparticle, &  ! diagnostics
    xbox_l, xbox_r, ybox_l, ybox_r, zbox_l, zbox_r, &
    Yee, global, harris, fxsho, nxcel, &  ! others
    rcorr, ishape, teti, post_process
      
    time_elapsed=0.; time_begin_array=0; time_end_array=0
    buffer_zone=0.  ! set to 0 anyway despite contained in input
    notime=1 ! notime=0 will output detailed timing
    !tracking_binary=.false.

    ! get the i/o directory names from the environment variable
    call get_environment_variable1(data_directory, len(data_directory))
    ! call get_environment_variable(data_directory, len(data_directory), 'DATA_DIRECTORY')
    data_directory = trim(adjustl(data_directory))//'/'
    call get_environment_variable2(restart_directory, len(restart_directory))
    ! call get_environment_variable(restart_directory, len(restart_directory), 'RESTART_DIRECTORY')
    restart_directory = trim(adjustl(restart_directory))//'/'
    restart_index_suffix(1) = '.1'
    restart_index_suffix(2) = '.2'

    my_short_int = myid
    call integer_to_character(myid_char,len(myid_char), my_short_int)
    if (myid_char == '') myid_char='0'

    ! read input deck
    if (myid == 0) then
      write(6,*) " "
      write(6,*) "Reading input file ..."
      open(5, file='input.f90', form='formatted', status='old')
      read(5, nml=datum, iostat=input_error)
      write(6, datum)
    endif
    iwt=0; t_stopped=0. ! set default values

    ! Broadcast info to all ranks
    ! global sim. info
    call MPI_BCAST(tmax                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(t_begin                ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(t_end                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(dtwci                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(dt                     ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(restart                ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(restrt_write           ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(quota                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
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
    call MPI_BCAST(nodey                  ,1     ,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nodez                  ,1     ,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
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
    call MPI_BCAST(uniform_loading_in_logical_grid,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(buffer_zone            ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(moat_zone              ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(profile_power          ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    ! field solver
    call MPI_BCAST(n_subcycles            ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nskipx                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nskipy                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nskipz                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(iterb                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(testorbt               ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(norbskip               ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    ! plasma setup
    call MPI_BCAST(nspec                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(qspec                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(wspec                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(frac                   ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(denmin                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(wpiwci                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(btspec                 ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(bete                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ieta                   ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(resis                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(netax                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(netay                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(netaz                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(etamin                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(etamax                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(eta_par                ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(anisot                 ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(gama                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ave1                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ave2                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(phib                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(smoothing              ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(smooth_coef            ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    ! init waves
    call MPI_BCAST(dB_B0                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(num_cycles             ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    ! diagnostic control
    call MPI_BCAST(nprint                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nwrtdata               ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nwrtrestart            ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nwrtparticle           ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xbox_l                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xbox_r                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ybox_l                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ybox_r                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zbox_l                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zbox_r                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    ! others 
    call MPI_BCAST(Yee                    ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(global                 ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(harris                 ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(fxsho                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nxcel                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(rcorr                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ishape                 ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(teti                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(post_process           ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

    ! The unit of dt is 1/wci in input file, but converted to 1/wpi here
    dt = dtwci * wpiwci
    ! field subcycling
    ! n_subcycles=max(n_subcycles,1_8)

    if (myid==0) then
      write(6,*)
      if (restart) then
        write(6,*) "*** Restarted from "//trim(adjustl(restart_directory))
      else
        write(6,*) "*** New run *** "
      endif
    endif 
    return
  end subroutine read_input


  !---------------------------------------------------------------------
  subroutine mpi_decomposition()
    implicit none
    integer :: i

    ! set MPI Cartesian geometry, define stride vector types, obtain new
    ! ID for the processors, perform 2D decomposition of the
    ! computational mesh, and find nearest neighbors (in y and z directions)
    ! specify decomposition (along y, z only; no decomposition along x) 
    ndim = 2
    dims(1) = nodey
    dims(2) = nodez

    ! create a division of processors in a cartesian grid
    ! where DIMS is an input/output parameter and an integer array of size ndims specifying 
    ! the number of nodes in each dimension, and a value of 0 indicates that MPI_Dims_create 
    ! should fill in a suitable value.
    call MPI_DIMS_CREATE(NUMPROCS, NDIM, DIMS, IERR)

    ! now npy means number of particles in each core along y
    npy = npy/dims(1)
    npz = npz/dims(2)
    if (myid == 0) then
      write(6,*)
      write(6,*) "MPI decompsition ..."
      do i = 1, ndim
        write(6,*) "  Dimension = ", i, " Dims = ", dims(i)
      enddo
    endif

    PERIODS = .TRUE. ! logical array of size ndims specifying whether the grid is periodic (true) or not (false) in each dimension
    REORDER = .TRUE. ! ranking may be reordered (true) or not (false) (logical)
    ! Makes a new communicator to which topology information has been attached
    call MPI_CART_CREATE(MPI_COMM_WORLD, NDIM, DIMS, PERIODS, REORDER, COMM2D, IERR) 
    ! Determines the rank of the calling process in the new communicator
    call MPI_COMM_RANK(COMM2D, MYID, IERR)
    ! Retrieves Cartesian topology information associated with the new communicator
    call MPI_CART_GET(COMM2D, NDIM, DIMS, PERIODS, COORDS, IERR)
    ! splits N elements between numprocs processors
    call MPE_DECOMP1D(NY, DIMS(1), COORDS(1), JB, JE)
    call MPE_DECOMP1D(NZ, DIMS(2), COORDS(2), KB, KE)
    ! print domain decomposition info
    ! write(6,*) 'myid=', myid, 'jb, je =', jb, je, 'kb, ke = ',kb, ke, 'coords =', coords

    nxmax  = nx + 2  
    nymax  = ny + 2
    nzmax  = nz + 2
    nylmax = je - jb + 1 
    nzlmax = ke - kb + 1  
    if (myid == 0) then
      write(6,*) "  Local array size in y-direction = ", nylmax
      write(6,*) "  Local array size in z-direction = ", nzlmax
    endif
    return
  end subroutine mpi_decomposition


  !---------------------------------------------------------------------
  ! Set global parameters
  subroutine set_parameters()
    implicit none
    integer :: i, j, k
    double_prec = 0.
    single_prec = 0.
    inquire (IOLENGTH=recl_for_double_precision) double_prec
    inquire (IOLENGTH=recl_for_real) single_prec

    if (myid==0) then
      write(6,*)
      write(6,*) "Setting up global parameters ..."
    endif

    ! nparbuf = nxmax*(nylmax+2)*(nzlmax+2)
    nspecm = nspec  ! nspecm is just a mirror of nspec
    npes = numprocs  ! npes is a copy of numprocs
    npes_over_60 = npes / 512  ! if numprocs > 512

    ! estimate on particle storage requirement
    nplmax = 0  ! max number of local particles in each process
    do i = 1, nspec
        nplmax = nplmax + npx(i)*npy(i)*npz(i)  ! notice npy, npz are already divided by nodey, nodez respectively
    enddo
    nplmax = 5*nplmax  ! pad storage requirement by a factor of 5

    ! number of tags used to track particles
    ! maxtags was initialized as 100
    maxtags_pe = maxtags/npes/nspec
    if (maxtags_pe==0) then
        maxtags_pe = 1  ! at least tracking one particle per species per pe
        maxtags = maxtags_pe*npes
    endif

    ! gathered enough info, now allocate arrays
    allocate (zbglobal(0:npes-1), zeglobal(0:npes-1), ybglobal(0:npes-1), yeglobal(0:npes-1)                     &
              ,kbglobal(0:npes-1), keglobal(0:npes-1), jbglobal(0:npes-1), jeglobal(0:npes-1)                     &
              ,nsendp(0:npes-1), nrecvp(0:npes-1), myid_stop(0:npes-1))
    allocate ( x(nplmax),y(nplmax),z(nplmax),vx(nplmax),vy(nplmax),vz(nplmax),link(nplmax),porder(nplmax)     &
              ,qp(nplmax))
    allocate (ptag(nplmax))
    allocate ( ninj(nspecm), ninj_global(nspecm),nescape(nspecm),nescape_global(nspecm),npart(nspecm)         &
              ,npart_global(nspecm),qleft(nspecm),qrite(nspecm))
    allocate ( nescape_xy(nspecm),nescape_yx(nspecm),nescape_xz(nspecm),nescape_zx(nspecm)                    &
              ,nescape_yz(nspecm),nescape_zy(nspecm)                                                          &
              ,nescape_xy_global(nspecm),nescape_yx_global(nspecm),nescape_xz_global(nspecm)                  &
              ,nescape_zx_global(nspecm),nescape_yz_global(nspecm),nescape_zy_global(nspecm))
    allocate ( x0(nspecm),x1(nspecm),tx0(nspecm),vpar(nspecm),vper(nspecm),vbal(nxmax,nspecm),bbal(nxmax))
    allocate ( dfac(nspecm),nskip(nspecm),ipleft(nspecm),iprite(nspecm),ipsendleft(nspecm),ipsendrite(nspecm) &
              ,iprecv(nspecm),ipsendtop(nspecm),ipsendbot(nspecm),ipsendlefttop(nspecm),ipsendleftbot(nspecm) &
              ,ipsendritetop(nspecm),ipsendritebot(nspecm),ipsend(nspecm))        
    allocate ( idmap_yz(0:ny+1,0:nz+1), idmap(0:nzmax), idfft(nzmax), kvec(nzlmax), jvec(nylmax) )

    do i = 1, nspecm
      qleft(i) = 0
      qrite(i) = 0
    enddo

    myid_stop(myid) = 0  

    !  Use CART_SHIFT to determine processor to immediate left (NBRLEFT) and right (NBRRITE) of processor MYID
    !  Since code is aperiodic in z, need to manually set the left boundary for processor 0 and right boundary for npes-1
    if (ndim == 2) then
        call MPI_CART_SHIFT(COMM2D,0,1,NBRLEFT,NBRRITE,IERR)
        call MPI_CART_SHIFT(COMM2D,1,1,NBRBOT ,NBRTOP ,IERR)
    else if (ndim == 1) then
        call MPI_CART_SHIFT(COMM2D,0,1,NBRLEFT,NBRRITE,IERR)
        NBRTOP = MYID
        NBRBOT = MYID
    else if (ndim == 0) then
        NBRLEFT = MYID
        NBRRITE = MYID
        NBRTOP = MYID
        NBRBOT = MYID
    endif

    call MPI_SENDRECV(NBRTOP    ,1,MPI_INTEGER ,NBRRITE,0,&
                      NBRLEFTTOP,1,MPI_INTEGER ,NBRLEFT,0,&
                      mpi_comm_world,status,ierr)
    call MPI_SENDRECV(NBRTOP    ,1,MPI_INTEGER ,NBRLEFT,0,&
                  NBRRITETOP,1,MPI_INTEGER ,NBRRITE,0,&
                  mpi_comm_world,status,ierr)
    call MPI_SENDRECV(NBRBOT    ,1,MPI_INTEGER ,NBRRITE,0,&
                  NBRLEFTBOT,1,MPI_INTEGER ,NBRLEFT,0,&
                  mpi_comm_world,status,ierr)
    call MPI_SENDRECV(NBRBOT    ,1,MPI_INTEGER ,NBRLEFT,0,&
                  NBRRITEBOT,1,MPI_INTEGER ,NBRRITE,0,&
                  mpi_comm_world,status,ierr) 

    ! recv, send id
    if (mod(coords(1),2) == 0.and.mod(coords(2),2) == 0) then
        isendid(1)=1
    else
        isendid(1)=0
    endif

    if (mod(coords(1)+1,2) == 0.and.mod(coords(2),2) == 0) then
        irecvid(1,1)=nbrrite
        irecvid(2,1)=-1
        irecvid(3,1)=nbrleft
        irecvid(4,1)=-1
    else if (mod(coords(1),2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,1)=-1
        irecvid(2,1)=nbrtop
        irecvid(3,1)=-1
        irecvid(4,1)=nbrbot
    else if (mod(coords(1)+1,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,1)=nbrritetop
        irecvid(2,1)=nbrlefttop
        irecvid(3,1)=nbrleftbot
        irecvid(4,1)=nbrritebot
    endif
    
    if (mod(coords(1)+1,2) == 0.and.mod(coords(2)  ,2) == 0) then
        isendid(2)=1
    else
        isendid(2)=0
    endif
    if (mod(coords(1)  ,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,2)=nbrrite
        irecvid(2,2)=-1
        irecvid(3,2)=nbrleft
        irecvid(4,2)=-1
    else if (mod(coords(1)+1,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,2)=-1
        irecvid(2,2)=nbrtop
        irecvid(3,2)=-1
        irecvid(4,2)=nbrbot
    else if (mod(coords(1)  ,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,2)=nbrritetop
        irecvid(2,2)=nbrlefttop
        irecvid(3,2)=nbrleftbot
        irecvid(4,2)=nbrritebot
    endif
    
    if (mod(coords(1)  ,2) == 0.and.mod(coords(2)+1,2) == 0) then
        isendid(3)=1
    else
        isendid(3)=0
    endif
    if (mod(coords(1)+1,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,3)=nbrrite
        irecvid(2,3)=-1
        irecvid(3,3)=nbrleft
        irecvid(4,3)=-1
    else if (mod(coords(1)  ,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,3)=-1
        irecvid(2,3)=nbrtop
        irecvid(3,3)=-1
        irecvid(4,3)=nbrbot
    else if (mod(coords(1)+1,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,3)=nbrritetop
        irecvid(2,3)=nbrlefttop
        irecvid(3,3)=nbrleftbot
        irecvid(4,3)=nbrritebot
    endif
    
    if (mod(coords(1)+1,2) == 0.and.mod(coords(2)+1,2) == 0) then
        isendid(4)=1
    else
        isendid(4)=0
    endif
    if (mod(coords(1)  ,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,4)=nbrrite
        irecvid(2,4)=-1
        irecvid(3,4)=nbrleft
        irecvid(4,4)=-1
    else if (mod(coords(1)+1,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,4)=-1
        irecvid(2,4)=nbrtop
        irecvid(3,4)=-1
        irecvid(4,4)=nbrbot
    else if (mod(coords(1)  ,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,4)=nbrritetop
        irecvid(2,4)=nbrlefttop
        irecvid(3,4)=nbrleftbot
        irecvid(4,4)=nbrritebot
    endif

    ! gather jb,je,kb,ke of each rank into *global (where *=jb,je,kb,ke)
    call MPI_ALLGATHER(jb,1,MPI_INTEGER8,jbglobal,1,MPI_INTEGER8,MPI_COMM_WORLD,IERR)
    call MPI_ALLGATHER(je,1,MPI_INTEGER8,jeglobal,1,MPI_INTEGER8,MPI_COMM_WORLD,IERR)
    call MPI_ALLGATHER(kb,1,MPI_INTEGER8,kbglobal,1,MPI_INTEGER8,MPI_COMM_WORLD,IERR)
    call MPI_ALLGATHER(ke,1,MPI_INTEGER8,keglobal,1,MPI_INTEGER8,MPI_COMM_WORLD,IERR)

    !VR: again, this is much simpler than the commented block
    do i = 0, numprocs-1
        do k = kbglobal(i), keglobal(i)
            do j = jbglobal(i), jeglobal(i)
              idmap_yz(j, k) = i
            enddo
        enddo
    enddo

    !VR: fill in ghost cells in idmap     
    idmap_yz(1:ny,0)    = idmap_yz(1:ny,nz)
    idmap_yz(1:ny,nz+1) = idmap_yz(1:ny,1)

    idmap_yz(0,1:nz)    = idmap_yz(ny,1:nz)
    idmap_yz(ny+1,1:nz) = idmap_yz(1,1:nz)

    idmap_yz(0,0)       = idmap_yz(ny,nz)
    idmap_yz(0,nz+1)    = idmap_yz(ny,1)
    idmap_yz(ny+1,0)    = idmap_yz(1,nz)
    idmap_yz(ny+1,nz+1) = idmap_yz(1,1)

    call MPI_TYPE_VECTOR(int(nzlmax+2,4), int(nx+2,4), int((nx+2)*(nylmax+2),4), MPI_DOUBLE_PRECISION, stridery, IERR)
    call MPI_TYPE_COMMIT(stridery, IERR)
    call MPI_TYPE_VECTOR(int(nylmax+2,4), int(nx+2,4), int(nx+2,4), MPI_DOUBLE_PRECISION, STRIDERZ, IERR)
    call MPI_TYPE_COMMIT(STRIDERZ, IERR)

    nptotp = 0  ! total number of particles per processor
    do i = 1, nspec
        nptotp = nptotp + npx(i)*npy(i)*npz(i)
    enddo
        
    if (.not.testorbt) norbskip=1

    call allocate_global_arrays
    return
  end subroutine set_parameters 
  
end module initialize