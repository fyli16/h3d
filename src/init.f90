module m_init
  use m_parameter
  use m_utils
  use m_io
  use m_mesh
  implicit none

  contains 

  !---------------------------------------------------------------------
  ! initialize simulation
  !---------------------------------------------------------------------
  subroutine init_sim
    ! read input deck
    call init_input
    
    ! MPI domain decomposition 
    call init_mpi_decomp
    
    ! allocate global arrays
    call allocate_arrays
    
    ! set up mesh 
    call setup_mesh

    ! open history diagnostic files
    call open_hist_files

    ! restart or a fresh start
    call makelist  ! see utils.f90
    if (restart) then 
      call init_restart 
    else ! fresh start
      call init_wavepart 
    endif 

    return 
  end subroutine init_sim


  !---------------------------------------------------------------------
  ! read in input file/parameters
  !---------------------------------------------------------------------
  subroutine init_input
    integer :: input_error

    namelist /datum/ &
    tmax, dtwci, restart, &   ! global info
    MPI_IO_format, &
    nx, ny, nz, xmax, ymax, zmax, npx, npy, npz, node_conf, periods, &  ! simulation domain
    xaa, xbb, nax, nbx, yaa, ybb, nay, nby, zaa, zbb, naz, nbz, &
    uniform_load_logical, &
    n_subcycles, nskipx, nskipy, nskipz, iterb, &  ! field solver
    nspec, n_sort, qspec, wspec, frac, denmin, wpiwci, btspec, bete, &  ! plasma setup
    ieta, resis, netax, netay, etamin, etamax, eta_par, eta_zs, &
    anisot, gamma, ave1, ave2, phib, smoothing, &
    dB_B0, num_wave_cycles, &  ! init waves
    n_print, n_write_mesh, n_write_energy, n_write_probes, & ! diagnostics
    n_write_tracking, n_write_restart, n_write_particle, &  
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
      print*, "*                                 H3D (V6.0)                               *"
      print*, "*                           YURI'S NONUNIFORM MESH                         *"
      print*, "*                           3D IMPLEMENTATION ONLY                         *"
      print*, "*                      UNIFORM LOADING IN PHYSICAL SPACE                   *"
      print*, "*               UNIFORM LOADING IN LOGICAL SPACE NOT YET IMPLEMENTED       *"
      print*, "***************************************************************************"
      print*, " "
      print*, " "
    endif 

    ! read in input deck
    if (myid == 0) then
      write(6,*) " "
      write(6,*) "Reading input file ..."
      open(5, file='input.f90', form='formatted', status='old')
      read(5, nml=datum, iostat=input_error)
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
    call MPI_BCAST(btspec                 ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(bete                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
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
    ! init waves
    call MPI_BCAST(dB_B0                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(num_wave_cycles        ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    ! diagnostic control
    call MPI_BCAST(n_print                ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_write_mesh           ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_write_energy         ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_write_probes         ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_write_tracking       ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_write_restart        ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_write_particle       ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
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

    ! field subcycling
    ! n_subcycles = max(n_subcycles, 1_8)

    ! set output directories
    data_directory = 'data/'
    restart_directory = 'restart/'
    restart_index_suffix(1) = '.1'
    restart_index_suffix(2) = '.2'

    ! print some simulation info
    if (myid==0) then
      write(6,*)
      if (.not. restart) write(6,*) "  *** New run ***"
      if (restart) write(6,*) "  *** Restart run ***"
      write(6,*)
      write(6,*) "  xmax, ymax, zmax = ", xmax, ymax, zmax
      write(6,*) "  nx, ny, nz = ", nx, ny, nz
      write(6,*) "  tmax, dtwci = ", tmax, dtwci
      write(6,*) "  npx = ", npx
      write(6,*) "  npy = ", npy
      write(6,*) "  npz = ", npz
      write(6,*) "  node_conf = ", node_conf
      write(6,*) "  xaa, xbb = ", xaa, xbb
      write(6,*) "  nax, nbx = ", nax, nbx 
      write(6,*) "  yaa, ybb = ", yaa, ybb
      write(6,*) "  nay, nby = ", nay, nby 
      write(6,*) "  zaa, zbb = ", zaa, zbb
      write(6,*) "  naz, nbz = ", naz, nbz 
    endif  

  end subroutine init_input


  !---------------------------------------------------------------------
  ! initialize MPI domain decomposition
  !---------------------------------------------------------------------
  subroutine init_mpi_decomp
    integer :: i

    ! set MPI Cartesian geometry, define stride vector types, 
    ! obtain new ID for the processors, perform 2D decomposition of the
    ! computational mesh, and find nearest neighbors (in y and z directions).
    ! specify decomposition (along y, z only; no decomposition along x) 
    if (nz==1 .and. ny==1) then ! only nx>1 and 1 rank will be used  
      ndim=1; dims(1)=1; dims(2)=1
    else if (nz == 1) then ! ny>1 and decomposition only occurs in y
      ndim=1; dims(1)=node_conf(1); dims(2)=1
    else ! ny>1, nz>1, and decomposition in both y and z
      ndim=2; dims=node_conf
    endif
    ! this is not needed because both components of dims have been explicitly specified
    ! call MPI_DIMS_CREATE(nprocs,NDIM,DIMS,IERR)

    ! npy now means number of particles in each rank along y
    ! npz now means number of particles in each rank along z
    npy=npy/dims(1); npz=npz/dims(2)

    ! print MPI decomposition information
    if (myid == 0) then
      write(6,*)
      write(6,*) "MPI decompsition ..."
      write(6,*) "  Total number of processors = ", nprocs
      do i = 1, ndim
        write(6,*) "  Dimension = ", i, " Dims = ", dims(i)
      enddo
    endif

    ! PERIODS = .TRUE. ! logical array of size ndims specifying whether the grid is periodic or not
    REORDER = .TRUE. ! ranking may be reordered (true) or not (false) (logical)
    ! Makes a new communicator to which topology information has been attached
    call MPI_CART_CREATE(MPI_COMM_WORLD, NDIM, DIMS, PERIODS, REORDER, COMM2D, IERR) 
    ! Determines the rank of the calling process in the new communicator
    call MPI_COMM_RANK(COMM2D, MYID, IERR)
    ! Retrieves Cartesian topology information (especially 'COORDS') associated with the new communicator
    call MPI_CART_GET(COMM2D, NDIM, DIMS, PERIODS, COORDS, IERR)
    ! splits N elements between nprocs processors
    call MPE_DECOMP1D(NY, DIMS(1), COORDS(1), JB, JE)
    call MPE_DECOMP1D(NZ, DIMS(2), COORDS(2), KB, KE)
    ! print domain decomposition info
    ! write(6,*) 'myid=', myid, 'jb, je =', jb, je, 'kb, ke = ',kb, ke, 'coords =', coords

    ! max number of cells. why adding 2?
    nxmax = nx + 2; nymax = ny + 2; nzmax = nz + 2

    ! local max number of cells
    nylmax = je - jb + 1 ; nzlmax = ke - kb + 1  
    if (myid == 0) then
      write(6,*) "  Local array size in x-direction = ", nx
      write(6,*) "  Local array size in y-direction = ", nylmax
      write(6,*) "  Local array size in z-direction = ", nzlmax
    endif

  end subroutine init_mpi_decomp


  !---------------------------------------------------------------------
  ! init waves and particles
  !---------------------------------------------------------------------
  subroutine init_wavepart
    use m_cal_eta
    use m_field
    use m_particle

    integer*8 :: ibp1, ibp2, i, remake, field_subcycle
    real*8 :: rxe, rye, rze, fxe, fye, fze, dtxi, dtyi, dtzi, &
              x_p, y_p, z_p, x_p_logical, y_p_logical, z_p_logical, &
              r_c, q_p, dtsav
    integer*8 :: ip, ipb1, ipb2, is, ixe, iye, ize, j, k
    real*8 :: vxa, vya, vza, vmag, th, ranval(4)
    real*8 :: x_pos, y_pos, z_pos, B0, VA, mi
    real*8 :: bx_, by_, bz_, ex_, ey_, ez_
    real*8 :: kx, ky, kz, kxmin, kymin, kzmin, dvx_, dvy_, dvz_, sin_factor
    real*8 :: w1e, w2e, w3e, w4e, w5e, w6e, w7e, w8e
    real*8 :: vix1, viy1, viz1, vix2, viy2, viz2, vix3, viy3, viz3, vix4, viy4, viz4, &
              vix5, viy5, viz5, vix6, viy6, viz6, vix7, viy7, viz7, vix8, viy8, viz8
    integer :: ixep1,iyep1,izep1
    integer :: tag, tag0

    real*8 :: loaded_percentage, print_percentage
    integer :: seed_size
    integer, allocatable :: seed(:)

    dtxi = one/meshX%dt ! dtxi=nx
    dtyi = one/meshY%dt ! dtyi=ny
    dtzi = one/meshZ%dt ! dtzi=nz

    ! initialize the random number generator with a different seed for different processors
    call random_seed(size=seed_size) ! obtain seed size
    allocate(seed(seed_size))
    ! VR: arbirary, but reproducible numbers for phases
    if (myid==0) then
      ! method 1: with a given seed
      seed = 31415  
      ! method 2: with seed generated by system
      ! call random_seed() 
      ! call random_seed(get=seed)
    endif
    call MPI_BCAST(seed, seed_size, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)  
    call random_seed(put=myid*seed) ! set current seed

    it = 0; itrestart = 0; 
    itstart = it; itfinish = tmax/dtwci

    ! some constants of mesh
    nx1 = nx+1; nx2 = nx+2
    ny1 = ny+1; ny2 = ny+2
    nz1 = nz+1; nz2 = nz+2
    hx = xmax/nx; hy = ymax/ny; hz = zmax/nz
    hxi = one/hx; hyi = one/hy; hzi = one/hz

    xb = zero; xe = xmax

    yb = meshY%xn(jb+1); ye = meshY%xn(je+2)
    do ipe = 0, nprocs-1
      ybglobal(ipe)=meshY%xn(jbglobal(ipe)+1)
      yeglobal(ipe)=meshY%xn(jeglobal(ipe)+2)
    enddo

    zb = meshZ%xn(kb+1); ze = meshZ%xn(ke+2)
    do ipe = 0, nprocs-1
      zbglobal(ipe)=meshZ%xn(kbglobal(ipe)+1)
      zeglobal(ipe)=meshZ%xn(keglobal(ipe)+2)
    enddo

    volume_fraction = (ye-yb)*(ze-zb)/(ymax*zmax)

    xb_logical = mesh_unmap(meshX,xb)
    xe_logical = mesh_unmap(meshX,xe)
    yb_logical = mesh_unmap(meshY,yb)
    ye_logical = mesh_unmap(meshY,ye)
    zb_logical = mesh_unmap(meshZ,zb)
    ze_logical = mesh_unmap(meshZ,ze)

    do is = 1, nspec
      npm = npx(is)*npy(is)*npz(is)*nprocs
      dfac(is)=real(ny*nz*nx)/real(npm)
      do ixe = 1, nx2 
        do iye = jb-1, je+1
          do ize = kb-1, ke+1
            qp_cell(ixe,iye,ize,is) = meshX%dxc(ixe) * meshY%dxc(iye+1) * meshZ%dxc(ize+1) * dfac(is)*frac(is)
          enddo
        enddo
      enddo
    enddo

    !---------------------------------------------------------------------
    ! initialie wave perturbation on the mesh 
    !---------------------------------------------------------------------
    if (myid==0) then
      print*, " "
      print*, "Initializing wave on the mesh ..."
    endif 

    B0 = one/wpiwci  ! RMS amplitude of the wave: B0=RMS(B)  
    mi=0.
    do j =1, nspec
      mi = mi + frac(j)*wspec(j)
    enddo
    VA = one/wpiwci/sqrt(mi) ! Alfven speed      

    kxmin = two*pi/xmax
    kymin = two*pi/ymax
    kzmin = two*pi/zmax
    kx = zero
    ky = zero
    kz = num_wave_cycles*kzmin

    bx = zero; by = zero; bz = zero
    ex = zero; ey = zero; ez = zero

    do k = kb-1, ke+1
      z_pos = meshZ%xc(k+1)
      do j = jb-1, je+1  
        y_pos = meshY%xc(j+1)
        do i = 1, nx2
          x_pos = meshX%xc(i) ! x has different indexing than y/z               

          ! single Alfven wave
          bx_ =  dB_B0*B0*sin(kz*z_pos)
          by_ = -dB_B0*B0*cos(kz*z_pos)
          bz_ = B0
          ex_ = zero  ! why e component is zero?
          ey_ = zero
          ez_ = zero
          dvx_ = -VA*bx_/B0 
          dvy_ = -VA*by_/B0 
          dvz_ = zero

          bx(i,j,k) = bx_
          by(i,j,k) = by_
          bz(i,j,k) = bz_
          ex(i,j,k) = ex_
          ey(i,j,k) = ey_
          ez(i,j,k) = ez_

          ! use vix to temporary store values of V on the grid
          vix(i,j,k) = dvx_
          viy(i,j,k) = dvy_
          viz(i,j,k) = dvz_
        enddo
      enddo
    enddo
    
    !--------------------------------------------------------------------- 
    ! load particles
    !---------------------------------------------------------------------
    if (myid==0) then
      print*, " "
      print*, "Initializing particles ..." 
    endif

    do is = 1, nspec
      ninj(is) = 0
      ninj_global(is) = 0
      npart(is) = 0
      tx0(is) = btspec(is)/(two*wpiwci**2)/wspec(is)
      x0(is) = zero
      x1(is) = xmax
      call MPI_ALLREDUCE(npart(is),npart_global(is),1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
    enddo
    te0 = bete/(two*wpiwci**2)
    vbal = one

    do is = 1, nspec
      tag0 = maxtags_pe*nspec*myid + (is-1)*maxtags_pe
      tag = 1 
      ipb1 = 1

      ! Nonuniform mesh
      if (uniform_load_logical) then
        ipb2 = npx(is)*npy(is)*npz(is)
      else
        ipb2 = npx(is)*npy(is)*npz(is)*nprocs*volume_fraction
      endif

      npm = npx(is)*npy(is)*npz(is)*nprocs
      dfac(is) = real(ny*nz)*(x1(is)-x0(is))/(hx*real(npm))
      vpar(is) = sqrt(btspec(is)/(wspec(is)*wpiwci*wpiwci))
      vper(is) = vpar(is)*sqrt(anisot(is))

      if (myid==0) then
        write(6,*) "  species #", is
        write(6,*) "  frac = ", frac(is)
        write(6,*) "  npx = ", npx(is)
        write(6,*) "  npy = ", npy(is)
        write(6,*) "  npz = ", npz(is)
        write(6,*) "  dfrac = ", dfac(is)
        write(6,*) " "
      endif

      ! actual particle loading
      print_percentage = zero
      do ip = ipb1, ipb2
        call random_number(harvest=ranval)
        if (uniform_load_logical) then
          x_p_logical  = xb_logical+(XE_LOGICAL-xb_logical)*ranval(1)
          y_p_logical  = yb_logical+(YE_LOGICAL-yb_logical)*ranval(2)
          z_p_logical  = zb_logical+(ZE_LOGICAL-zb_logical)*ranval(3)
          x_p          = mesh_map(meshX,x_p_logical)
          y_p          = mesh_map(meshY,y_p_logical)
          z_p          = mesh_map(meshZ,z_p_logical)
          ixe          = dtxi*x_p_logical+1.50000000000d+00
          iye          = dtyi*y_p_logical+1.50000000000d+00
          ize          = dtzi*z_p_logical+1.50000000000d+00
          q_p          = meshX%dxc(ixe) * meshY%dxc(iye) * meshZ%dxc(ize) * dfac(is)*frac(is)
        else
          x_p  = X0(IS)+(X1(IS)-X0(IS))*ranval(1)
          y_p  = YB+(YE-YB)*ranval(2)
          z_p  = ZB+(ZE-ZB)*ranval(3)
          q_p  = hx*hy*hz*dfac(is)*frac(is)
        endif

        np     = ipstore
        x(np)  = x_p
        y(np)  = y_p
        z(np)  = z_p
        qp(np) = q_p
        if (tag<=maxtags_pe .and. tag0+tag<=maxtags) then ! only track first maxtags particles
          ptag(np)=tag0+tag
        else
          ptag(np)=0 ! do not track
        endif
        tag = tag+1

        ! Nonuniform mesh - using mesh_unmap
        rxe=dtxi*mesh_unmap(meshX,x(np))+1.50000000000d+00
        rye=dtyi*mesh_unmap(meshY,y(np))+1.50000000000d+00
        rze=dtzi*mesh_unmap(meshZ,z(np))+1.50000000000d+00
        ixe=rxe
        iye=rye
        ize=rze
        iye=iye-1  ! integer index in y direction starts at 0
        ize=ize-1  ! integer index in z direction starts at 0

        fxe=rxe-ixe
        fye=rye-iye
        fze=rze-ize
        ixep1 = ixe+1
        iyep1 = iye+1
        izep1 = ize+1

        w1e=(1.-fxe)*(1.-fye)*(1.-fze)
        w2e=fxe*(1.-fye)*(1.-fze)
        w3e=(1.-fxe)*fye*(1.-fze)
        w4e=fxe*fye*(1.-fze)
        w5e=(1.-fxe)*(1.-fye)*fze
        w6e=fxe*(1.-fye)*fze
        w7e=(1.-fxe)*fye*fze
        w8e=fxe*fye*fze
        
        vix1=vix(ixe  ,iye  ,ize  )
        vix2=vix(ixep1,iye  ,ize  )
        vix3=vix(ixe  ,iyep1,ize  )
        vix4=vix(ixep1,iyep1,ize  )
        vix5=vix(ixe  ,iye  ,izep1)
        vix6=vix(ixep1,iye  ,izep1)
        vix7=vix(ixe  ,iyep1,izep1)
        vix8=vix(ixep1,iyep1,izep1)
        viy1=viy(ixe  ,iye  ,ize  )
        viy2=viy(ixep1,iye  ,ize  )
        viy3=viy(ixe  ,iyep1,ize  )
        viy4=viy(ixep1,iyep1,ize  )
        viy5=viy(ixe  ,iye  ,izep1)
        viy6=viy(ixep1,iye  ,izep1)
        viy7=viy(ixe  ,iyep1,izep1)
        viy8=viy(ixep1,iyep1,izep1)
        viz1=viz(ixe  ,iye  ,ize  )
        viz2=viz(ixep1,iye  ,ize  )
        viz3=viz(ixe  ,iyep1,ize  )
        viz4=viz(ixep1,iyep1,ize  )
        viz5=viz(ixe  ,iye  ,izep1)
        viz6=viz(ixep1,iye  ,izep1)
        viz7=viz(ixe  ,iyep1,izep1)
        viz8=viz(ixep1,iyep1,izep1)
        
        dvx_ = w1e*vix1+w2e*vix2+w3e*vix3+w4e*vix4 + w5e*vix5+w6e*vix6+w7e*vix7+w8e*vix8  
        dvy_ = w1e*viy1+w2e*viy2+w3e*viy3+w4e*viy4 + w5e*viy5+w6e*viy6+w7e*viy7+w8e*viy8  
        dvz_ = w1e*viz1+w2e*viz2+w3e*viz3+w4e*viz4 + w5e*viz5+w6e*viz6+w7e*viz7+w8e*viz8  
                  
        ! interpolate V at the particle position from pre-computed values at the grid
        ! dvx_ = -dB_B0*VA*sin(kz*z_p)
          
        call random_number(harvest=ranval)
        vmag = sqrt(-log(one-ranval(1)))
        th = two*pi*ranval(2)
        vza = vpar(is)*vmag*cos(th) + dvz_
        
        vmag = sqrt(-log(one-ranval(3)))
        th = two*pi*ranval(4)
        vxa = vper(is)*vmag*cos(th) + dvx_
        vya = vper(is)*vmag*sin(th) + dvy_

        vx(np) = vxa
        vy(np) = vya
        vz(np) = vza

        ipstore = link(np)
        link(np) = iphead(ixe,iye,ize,is)
        iphead(ixe,iye,ize,is) = np

        loaded_percentage = 100.0*real(ip-ipb1)/(ipb2-ipb1)
        
        if (myid==0 .and. loaded_percentage>=print_percentage) then
            write(6,"(A,F5.1,A)") "   loaded ", loaded_percentage," % of particles"
            print_percentage = print_percentage + 20.0d0
        endif

      enddo
    enddo

    ! what's doing here?
    if (ndim /= 1) then
      call xrealbcc(ex,1_8,nx,ny,nz)
      call xrealbcc(ey,1_8,nx,ny,nz)
      call xrealbcc(ez,1_8,nx,ny,nz)
    else
      call xrealbcc_pack_e_2d(ex,ey,ez,1_8,nx,ny,nz)
    endif
    
    ! set the friction force and resitivity (the former can be zero at t=0)
    if(myid==0) then
      print*, " "
      print*, '  setting friction force and resistivity'
    endif

    do k = kb-1,ke+1
      do j = jb-1,je+1
        do i = 1,nx2
          fox(i,j,k) = zero
          foy(i,j,k) = zero
          foz(i,j,k) = zero
          eta(i,j,k) = resis
        enddo
      enddo
    enddo

    ! what's done here
    if(myid==0) then
      print*, " "
      print*, "  setting dt=0 temporarily and call 'trans'"
    endif
    dtsav = dt
    dt    = zero
    call trans  ! because dt=0, no actual push is done (see 'parmov')
    dt    = dtsav

    ! advance field if n_subcyles>=1
    ! since currently n_subcycles==0, this block is skipped
    do field_subcycle = 1, n_subcycles 
      if (ndim /= 1) then
        call field
      else
        call field_2d
      endif
    enddo

    ! calculate resistivity (Dietmar's resistivity)
    if (ndim /= 1) then
        call cal_eta  
    else
        call cal_eta_2d
    endif
    
    deallocate(seed) 

  end subroutine init_wavepart


  !---------------------------------------------------------------------
  ! initialize restart from a previous run
  !---------------------------------------------------------------------
  subroutine init_restart

    integer*8 :: ixe, iye, ize, i, j, k

    ! read in restart data set and corresponding step of iteration
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

    ! what's the purpose of nprocs_over_60
    do i = 0, nprocs_over_60 
      if (mod(myid,nprocs_over_60+1) .eq. i) then 
          call write_read_restart_files(-1.0)  ! read restart data
          call MPI_BCAST(itrestart,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
      endif
    enddo

    ! determine start, finish number of steps
    itstart = itrestart; it = itstart
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
    xb_logical = mesh_unmap(meshX,xb)
    xe_logical = mesh_unmap(meshX,xe)
    yb_logical = mesh_unmap(meshY,yb)
    ye_logical = mesh_unmap(meshY,ye)
    zb_logical = mesh_unmap(meshZ,zb)
    ze_logical = mesh_unmap(meshZ,ze)
              
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
      
    if (myid==0) then
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

end module m_init 