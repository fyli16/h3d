! Declare global parameters and global arrays in this module
module parameter_mod  
  use mpi
  implicit none
  save

  integer :: my_short_int, i_source, i_destination, i_tag, i_length, i_i
  integer, dimension(8,128) :: time_begin_array, time_end_array
  real*8, dimension(128):: time_elapsed
  integer*8 :: nxmax, nymax, nzmax, nspecm, npes, nvar, nylmax, nzlmax, npm, npes_over_60
  integer :: numprocs, ndim, dims(2), nodey, nodez, ierr, comm2d, myid, req(8)     &
            ,nbrtop, nbrbot, nbrritetop, nbrlefttop, nbrritebot, nbrleftbot        &
            ,nbrleft, nbrrite, ipe, stridery, striderz, iseed(1), coords(2)
  integer :: status(mpi_status_size),status1(mpi_status_size),status2(mpi_status_size),status_array(mpi_status_size,8)
  real*8 :: zb,ze,yb,ye,teti,volume_fraction,cell_volume_ratio,zb_logical,ze_logical,yb_logical,ye_logical    &
            ,xb_logical,xe_logical,xb,xe,smooth_coef
  real*8, dimension(:), allocatable :: zbglobal,zeglobal,ybglobal,yeglobal,xc_uniform,yc_uniform,zc_uniform   &
                                    ,xv_uniform,yv_uniform,zv_uniform
  integer*8, dimension(:), allocatable :: kbglobal, keglobal, jbglobal, jeglobal, nsendp,nrecvp,ixc_2_c_map,iyc_2_c_map,izc_2_c_map   &
                                          ,ixc_2_v_map,iyc_2_v_map,izc_2_v_map     &
                                          ,ixv_2_c_map,iyv_2_c_map,izv_2_c_map     &
                                          ,ixv_2_v_map,iyv_2_v_map,izv_2_v_map     
  real*8, dimension(:,:,:), allocatable :: ex,ey,ez,bx,by,bz,fox,foy,foz,eta,curlex,curley,curlez,&
                                          bx_av, by_av, bz_av, &
                                          bxs,bys,bzs,den,deno,denh,dpedx,dpedy,dpedz,vix,viy,viz,vixo,viyo,   &
                                          vizo,pe,curlbx,curlby,curlbz,                               &
                                          eta_times_b_dot_j
  real*8, dimension(:,:,:,:), allocatable :: dns, dnsh, vxs, vys, vzs, tpar, tperp, qp_cell
  real*8, dimension(:,:,:,:), allocatable :: p_xx,p_xy,p_xz,p_yy,p_yz,p_zz
  real*8, dimension(:,:), allocatable :: ainjxz,ainjzx,deavxz,deavzx,vxavxz,vyavxz,vzavxz,vxavzx,      &
                                        vyavzx,vzavzx,vxcaxz,vycaxz,vzcaxz,vxcazx,vycazx,vzcazx,      &
                                        ainjyz,ainjzy,deavyz,deavzy,vxavyz,vyavyz,vzavyz,vxavzy,      &
                                        vyavzy,vzavzy,vxcayz,vycayz,vzcayz,vxcazy,vycazy,vzcazy,      &
                                        ainjxy,ainjyx,deavxy,deavyx,vxavxy,vyavxy,vzavxy,vxavyx,      &
                                        vyavyx,vzavyx,vxcaxy,vycaxy,vzcaxy,vxcayx,vycayx,vzcayx
  real*8, dimension(:), allocatable :: x, y, z, vx, vy, vz, qp
  integer, dimension(:), allocatable :: ptag ! tag used to trace particles
  integer*8, dimension(:), allocatable :: link,porder
  integer*8 :: nplmax, ipstore, np, n_subcycles
  integer*8, dimension(:,:,:,:), allocatable:: iphead, iptemp
  integer*8, dimension(:), allocatable ::  ninj,ninj_global,nescape,nescape_global,npart,npart_global
  integer*8, dimension(:),allocatable :: nescape_yz,nescape_zy,nescape_xy                 &
                                    ,nescape_yx,nescape_xz,nescape_zx                 &
                                    ,nescape_yz_global,nescape_zy_global              &
                                    ,nescape_xy_global,nescape_yx_global              &
                                    ,nescape_xz_global,nescape_zx_global
  real*8, dimension(:), allocatable :: qleft,qrite
  real*8, dimension(:), allocatable:: x0,x1,tx0,vpar,vper,bbal
  real*8, dimension(:,:), allocatable :: vbal
  real*8, dimension(5) :: rcorr
  integer*8, dimension(5) :: ishape
  real*8, dimension(5) :: btspec, qspec, wspec, frac, anisot
  real*8 :: denmin, resis, wpiwci, bete, fxsho,ave1,ave2,phib,demin2, &
            xmax, ymax, zmax, dt, gama, dtwci, wall_clock_elapsed,tmax,buffer_zone,  &
            xaa, xbb, yaa, ybb, zaa, zbb, t_stopped=0.
  integer*8 :: nax,nbx,nay,nby,naz,nbz
  integer*8, dimension(8) :: wall_clock_begin,wall_clock_end
  integer*8 :: eta_par, nparbuf
  integer*8, dimension(5) :: npx, npy, npz
  integer*8 :: iterb, norbskip, nxcel, netax, netay, netaz, nspec, nx, ny, nz, n_print, &
            n_write_data, n_write_particle, n_write_restart, nskipx,nskipy,nskipz
  real*8 :: etamin,etamax,moat_zone
  integer*8 :: ieta, profile_power
  logical :: testorbt, restart, uniform_loading_in_logical_grid, MPI_IO_format, smoothing, &
              global, harris, Yee, print_info, write_data, write_restart
  real*8 ::  hx, hy, hz, hxi, hyi, hzi, efld, bfld, efluidt, ethermt, eptclt, time, te0
  integer :: it, notime
  integer*8 :: nsteps0, itfin, iwt=0, nx1, nx2, ny1, ny2, nz1, nz2, iopen, file_unit(25), file_unit_time,            &
                file_unit_tmp,file_unit_read(20),nptot,npleaving,npentering,iclock_speed, nptotp
  real*8 :: clock_time_init, clock_time_old, clock_time, clock_time1
  real*8, dimension(:), allocatable :: dfac
  integer*8, dimension(:), allocatable :: nskip,ipleft,iprite,ipsendleft,ipsendrite,iprecv,ipsendtop,ipsendbot     &
                                      ,ipsendlefttop,ipsendleftbot,ipsendritetop,ipsendritebot,ipsend
  integer*8:: idum
  integer*8, dimension(:), allocatable:: idmap
  integer*8, dimension(:,:), allocatable:: idmap_yz
  integer*8 :: kb, ke, jb, je, nsendtotp, nrecvtotp, nsendtot, nrecvtot
  integer*8, dimension(:), allocatable :: idfft, kvec, jvec, myid_stop
  integer*8 :: ihstb, ihste, isendid(4), irecvid(4,4)
  real*4 :: single_prec
  real*8 :: double_prec, xtmp1m, xtmp2m, xbox_l, xbox_r, ybox_l, ybox_r, zbox_l, zbox_r, t_begin, t_end
  real*8, dimension(:,:), allocatable :: buf, buf2, buf_p1
  real*8, dimension(:,:,:), allocatable :: buf_particle
  integer, dimension(:), allocatable :: buftime
  integer, parameter :: nprobes=6, nbufsteps=100, tracking_width=14
  integer :: maxtags=100, maxtags_pe, ntot 
  logical :: tracking_binary, tracking_mpi
  integer :: tracking_fh

  integer*8 :: recl_for_real, recl_for_double_precision
  logical :: periods(2), reorder
  integer*8 :: restart_index
  character(len=2) :: restart_index_suffix(2)
  character(len=160) :: data_directory, restart_directory, cycle_ascii, cycle_ascii_new, &
                      myid_char, cleanup_status
  real*8 :: dB_B0, num_cycles ! for init_wave
  real*8, parameter :: zero=0.0d0, one=1.0d0, two=2.0d0, one_half=0.5d0, pi=acos(-1.)

  contains
  !---------------------------------------------------------------------
  subroutine read_input()
    implicit none

    integer :: input_error

    namelist /datum/ &
    tmax, t_begin, t_end, dtwci, restart, &   ! global info
    MPI_IO_format, &
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
    n_print, n_write_data, write_restart, n_write_restart, n_write_particle, &  ! diagnostics
    tracking_binary, tracking_mpi, xbox_l, xbox_r, ybox_l, ybox_r, zbox_l, zbox_r, notime, &
    Yee, global, harris, fxsho, nxcel, &  ! others
    rcorr, ishape, teti

    ! ??
    my_short_int = myid
    call integer_to_character(myid_char, len(myid_char), my_short_int)
    if (myid_char=='') myid_char='0'

    ! read input deck
    if (myid == 0) then
      write(6,*) " "
      write(6,*) "Reading input file ..."
      open(5, file='input.f90', form='formatted', status='old')
      read(5, nml=datum, iostat=input_error)
      write(6, datum)
    endif

    ! Broadcast input parameters (read in at rank 0) to all other ranks
    ! global sim. info
    call MPI_BCAST(tmax                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(t_begin                ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(t_end                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(dtwci                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(restart                ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(write_restart          ,1     ,MPI_LOGICAL        ,0,MPI_COMM_WORLD,IERR)
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
    call MPI_BCAST(n_print                ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_write_data           ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_write_restart        ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_write_particle       ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(tracking_binary        ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(tracking_mpi           ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(Yee                    ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xbox_l                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xbox_r                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ybox_l                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ybox_r                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zbox_l                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zbox_r                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(notime                 ,1     ,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
    ! others 
    call MPI_BCAST(Yee                    ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(global                 ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(harris                 ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(fxsho                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nxcel                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(rcorr                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ishape                 ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(teti                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

    ! The unit of dt is 1/wci in input file, 
    ! but now converted to 1/wpi inside the code
    dt = dtwci * wpiwci
    ! field subcycling
    ! n_subcycles=max(n_subcycles,1_8)
    data_directory = 'data/'
    restart_directory = 'restart/'
    restart_index_suffix(1) = '.1'
    restart_index_suffix(2) = '.2'

    if (myid==0) then
      write(6,*)
      if (restart) then
        write(6,*) "*** Restarted run ***"
      else
        write(6,*) "*** New run ***"
      endif
    endif 

    return
  end subroutine read_input


  !---------------------------------------------------------------------
  subroutine domain_decomposition()
    implicit none

    integer :: i

    ! set MPI Cartesian geometry, define stride vector types, obtain new
    ! ID for the processors, perform 2D decomposition of the
    ! computational mesh, and find nearest neighbors (in y and z directions).
    ! specify decomposition (along y, z only; no decomposition along x) 
    ndim = 2; dims(1) = nodey; dims(2) = nodez
    ! create division of processors in a cartesian grid
    call MPI_DIMS_CREATE(NUMPROCS, NDIM, DIMS, IERR)

    ! now npy means number of particles in each core along y
    npy = npy/dims(1); npz = npz/dims(2)

    if (myid == 0) then
      write(6,*)
      write(6,*) "MPI decompsition ..."
      do i = 1, ndim
        write(6,*) "  Total number of processors = ", numprocs
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

    nxmax  = nx + 2  ! why adding 2?
    nymax  = ny + 2
    nzmax  = nz + 2
    nylmax = je - jb + 1 
    nzlmax = ke - kb + 1  
    if (myid == 0) then
      write(6,*) "  Local array size in x-direction = ", nx
      write(6,*) "  Local array size in y-direction = ", nylmax
      write(6,*) "  Local array size in z-direction = ", nzlmax
    endif

    return
  end subroutine domain_decomposition


  !---------------------------------------------------------------------
  ! Set global parameters
  subroutine allocate_global_arrays()
    implicit none

    integer :: i, j, k

    double_prec = 0.
    single_prec = 0.
    inquire (IOLENGTH=recl_for_double_precision) double_prec
    inquire (IOLENGTH=recl_for_real) single_prec

    if (myid==0) then
      write(6,*) " "
      write(6,*) "Setting up global arrays ..."
    endif

    ! nparbuf = nxmax*(nylmax+2)*(nzlmax+2)
    nspecm = nspec  ! nspecm is just a mirror copy of nspec
    npes = numprocs  ! npes is just a mirror copy of numprocs
    npes_over_60 = npes/512  ! if numprocs > 512?

    ! estimate on particle storage requirement
    nptotp = 0  ! total number of particles per rank
    do i = 1, nspec
      nptotp = nptotp + npx(i)*npy(i)*npz(i)
    enddo
    nplmax = 2* nptotp  ! pad storage requirement by a factor of 2

    ! number of tags used to track particles per species per rank
    ! maxtags was initialized as 100
    maxtags_pe = maxtags/npes/nspec
    if (maxtags_pe==0) then
        maxtags_pe = 1  ! at least tracking one particle per species per rank
        maxtags = maxtags_pe * npes
    endif

    ! gathered enough info, now allocate arrays
    allocate (zbglobal(0:npes-1), zeglobal(0:npes-1), ybglobal(0:npes-1), yeglobal(0:npes-1)   &
              ,kbglobal(0:npes-1), keglobal(0:npes-1), jbglobal(0:npes-1), jeglobal(0:npes-1)  &
              ,nsendp(0:npes-1), nrecvp(0:npes-1), myid_stop(0:npes-1))
    allocate (x(nplmax), y(nplmax), z(nplmax), vx(nplmax), vy(nplmax), vz(nplmax),   &
            link(nplmax), porder(nplmax), qp(nplmax), ptag(nplmax))
    allocate ( ninj(nspecm), ninj_global(nspecm),nescape(nspecm),nescape_global(nspecm),npart(nspecm)  &
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
        
    if (.not.testorbt) norbskip=1

    allocate ( ex       (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),ey       (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,ez       (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),bx       (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,by       (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),bz       (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,bx_av    (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),by_av    (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,bz_av    (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)  &
              ,fox      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),foy      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,foz      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),eta      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,curlex   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),curley   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,curlez   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),bxs      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,bys      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),bzs      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,den      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),deno     (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,denh     (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,dpedx    (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),dpedy    (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,dpedz    (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),vix      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,viy      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),viz      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,vixo     (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),viyo     (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,vizo     (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),pe       (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,curlbx   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),curlby   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
              ,curlbz   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)                                                  &
              ,eta_times_b_dot_j(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax))
    allocate ( dns(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm),vxs   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm) &
              ,dnsh(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm) &
              ,vys(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm),vzs   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm) &
              ,tpar(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm),tperp(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm) &
              ,qp_cell(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm)) 
    allocate ( p_xx(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm),p_xy(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm)&
              ,p_xz(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm),p_yy(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm)&
              ,p_yz(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm),p_zz(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm) )
    allocate ( ainjxz(nxmax,kb-1:kb+nzlmax),ainjzx(nxmax,kb-1:kb+nzlmax),deavxz(nxmax,kb-1:kb+nzlmax)          &
              ,deavzx(nxmax,kb-1:kb+nzlmax),vxavxz(nxmax,kb-1:kb+nzlmax),vyavxz(nxmax,kb-1:kb+nzlmax)          &
              ,vzavxz(nxmax,kb-1:kb+nzlmax),vxavzx(nxmax,kb-1:kb+nzlmax),vyavzx(nxmax,kb-1:kb+nzlmax)          &
              ,vzavzx(nxmax,kb-1:kb+nzlmax),vxcaxz(nxmax,kb-1:kb+nzlmax),vycaxz(nxmax,kb-1:kb+nzlmax)          &
              ,vzcaxz(nxmax,kb-1:kb+nzlmax),vxcazx(nxmax,kb-1:kb+nzlmax),vycazx(nxmax,kb-1:kb+nzlmax)          &
              ,vzcazx(nxmax,kb-1:kb+nzlmax))
    allocate ( ainjyz(jb-1:jb+nylmax,kb-1:kb+nzlmax),ainjzy(jb-1:jb+nylmax,kb-1:kb+nzlmax)                     &
              ,deavyz(jb-1:jb+nylmax,kb-1:kb+nzlmax),deavzy(jb-1:jb+nylmax,kb-1:kb+nzlmax)                     &
              ,vxavyz(jb-1:jb+nylmax,kb-1:kb+nzlmax),vyavyz(jb-1:jb+nylmax,kb-1:kb+nzlmax)                     &
              ,vzavyz(jb-1:jb+nylmax,kb-1:kb+nzlmax),vxavzy(jb-1:jb+nylmax,kb-1:kb+nzlmax)                     &
              ,vyavzy(jb-1:jb+nylmax,kb-1:kb+nzlmax),vzavzy(jb-1:jb+nylmax,kb-1:kb+nzlmax)                     &
              ,vxcayz(jb-1:jb+nylmax,kb-1:kb+nzlmax),vycayz(jb-1:jb+nylmax,kb-1:kb+nzlmax)                     &
              ,vzcayz(jb-1:jb+nylmax,kb-1:kb+nzlmax),vxcazy(jb-1:jb+nylmax,kb-1:kb+nzlmax)                     &
              ,vycazy(jb-1:jb+nylmax,kb-1:kb+nzlmax),vzcazy(jb-1:jb+nylmax,kb-1:kb+nzlmax))
    allocate ( ainjxy(nxmax,jb-1:jb+nylmax),ainjyx(nxmax,jb-1:jb+nylmax),deavxy(nxmax,jb-1:jb+nylmax)          &
              ,deavyx(nxmax,jb-1:jb+nylmax),vxavxy(nxmax,jb-1:jb+nylmax),vyavxy(nxmax,jb-1:jb+nylmax)          &
              ,vzavxy(nxmax,jb-1:jb+nylmax),vxavyx(nxmax,jb-1:jb+nylmax),vyavyx(nxmax,jb-1:jb+nylmax)          &
              ,vzavyx(nxmax,jb-1:jb+nylmax),vxcaxy(nxmax,jb-1:jb+nylmax),vycaxy(nxmax,jb-1:jb+nylmax)          &
              ,vzcaxy(nxmax,jb-1:jb+nylmax),vxcayx(nxmax,jb-1:jb+nylmax),vycayx(nxmax,jb-1:jb+nylmax)          &
              ,vzcayx(nxmax,jb-1:jb+nylmax))
    allocate (iphead(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm),iptemp(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm))
    allocate (xc_uniform(nxmax),yc_uniform(nymax),zc_uniform(nzmax),xv_uniform(nxmax),yv_uniform(nymax),zv_uniform(nzmax))
    allocate (ixc_2_c_map(nx+1),iyc_2_c_map(ny+1),izc_2_c_map(nz+1))
    allocate (ixc_2_v_map(nx+1),iyc_2_v_map(ny+1),izc_2_v_map(nz+1))
    allocate (ixv_2_c_map(nx+1),iyv_2_c_map(ny+1),izv_2_c_map(nz+1))
    allocate (ixv_2_v_map(nx+1),iyv_2_v_map(ny+1),izv_2_v_map(nz+1))
    allocate (buf(nprobes,nbufsteps),buf2(nprobes*npes,nbufsteps),buftime(nbufsteps))
    allocate (buf_particle(tracking_width,nspec*maxtags,nbufsteps))
    allocate (buf_p1(tracking_width,nspec*maxtags))

    return
  end subroutine allocate_global_arrays 

end module parameter_mod
