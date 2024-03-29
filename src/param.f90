! Declare global parameters and global arrays in this module
module m_parameter  
  use mpi
  implicit none

  save

  ! global 
  logical :: restart, uniform_load_logical, MPI_IO_format 
  logical :: periods(2), reorder=.true.

  integer :: it, itstart, itfinish, itrestart, my_short_int 
  integer :: i_source, i_tag, i_length, i_i
  integer :: time_begin(8,128), time_end(8,128), ierr
  real*8 :: time_elapsed(128)
  real*8 :: clock_init, clock_old, clock_now, clock_time1

  integer :: nprocs, ndim, dims(2), node_conf(2), comm2d, myid, req(8), & 
            nbrtop, nbrbot, nbrritetop, nbrlefttop, nbrritebot, nbrleftbot, &      
            nbrleft, nbrrite, ipe, stridery, striderz, iseed(1), coords(2)
  character(len=160) :: myid_char

  integer :: status(mpi_status_size), status1(mpi_status_size), &
            status2(mpi_status_size), status_array(mpi_status_size,8)

  integer*8:: recl_for_single, recl_for_double
  real :: single_prec
  real*8 :: double_prec

  integer*8 :: restart_index=1
  character(len=2) :: restart_index_suffix(2)
  character(len=160) :: data_directory, restart_directory, cycle_ascii, cycle_ascii_new

  ! mesh 
  integer*8 :: nxmax, nymax, nzmax, nvar, nylmax, nzlmax
  
  real*8 :: xb, xe, yb, ye, zb, ze, volume_fraction, cell_volume_ratio, &
      zb_logical, ze_logical, yb_logical, ye_logical, xb_logical, xe_logical

  real*8, dimension(:), allocatable :: zbglobal, zeglobal, ybglobal, yeglobal, &
                                      xc_uniform, yc_uniform, zc_uniform, &
                                      xv_uniform, yv_uniform, zv_uniform

  integer*8, dimension(:), allocatable :: kbglobal, keglobal, jbglobal, jeglobal, nsendp, nrecvp, &
                          ixc_2_c_map, iyc_2_c_map, izc_2_c_map, & 
                          ixc_2_v_map, iyc_2_v_map, izc_2_v_map, &
                          ixv_2_c_map, iyv_2_c_map, izv_2_c_map, &
                          ixv_2_v_map, iyv_2_v_map, izv_2_v_map  

  real*8, dimension(:,:,:), allocatable :: uniform_mesh, ex, ey, ez, bx, by, bz, fox, foy, foz, &
      eta, curlex, curley, curlez, bx_av, by_av, bz_av, bxs, bys, bzs, den, deno, denh, &
      dpedx, dpedy, dpedz, vix, viy, viz, vixo, viyo, vizo, pe, curlbx, curlby, curlbz, eta_times_b_dot_j

  real*8, dimension(:,:,:,:), allocatable :: dns, dnsh, vxs, vys, vzs, tpar, tperp, qp_cell

  real*8, dimension(:,:,:,:), allocatable :: p_xx, p_xy, p_xz, p_yy, p_yz, p_zz

  real*8, dimension(:,:), allocatable :: ainjxz,ainjzx,deavxz,deavzx,vxavxz,vyavxz,vzavxz,vxavzx,      &
                                         vyavzx,vzavzx,vxcaxz,vycaxz,vzcaxz,vxcazx,vycazx,vzcazx,      &
                                         ainjyz,ainjzy,deavyz,deavzy,vxavyz,vyavyz,vzavyz,vxavzy,      &
                                         vyavzy,vzavzy,vxcayz,vycayz,vzcayz,vxcazy,vycazy,vzcazy,      &
                                         ainjxy,ainjyx,deavxy,deavyx,vxavxy,vyavxy,vzavxy,vxavyx,      &
                                         vyavyx,vzavyx,vxcaxy,vycaxy,vzcaxy,vxcayx,vycayx,vzcayx


  ! particles
  real*8, dimension(:), allocatable :: x, y, z, vx, vy, vz, qp

  integer*8, dimension(:), allocatable :: ptag, link, porder

  integer*8 :: nplmax, ipstore, np, npm 

  integer*8, dimension(:,:,:,:), allocatable:: iphead, iptemp
  integer*8, dimension(:), allocatable ::  ninj, ninj_global, nescape,nescape_global, npart, npart_global
  integer*8, dimension(:),allocatable :: nescape_yz,nescape_zy,nescape_xy                 &
                                        ,nescape_yx,nescape_xz,nescape_zx                 &
                                        ,nescape_yz_global,nescape_zy_global              &
                                        ,nescape_xy_global,nescape_yx_global              &
                                        ,nescape_xz_global,nescape_zx_global

  real*8, dimension(:), allocatable:: x0,x1,tx0,vpar,vper

  real*8, dimension(5) :: beta_spec, qspec, wspec, frac, anisot

  real*8 :: denmin, resis, wpiwci, beta_elec, &
            xmax, ymax, zmax, dt, gamma, dtwci, wall_clock_elapsed, tmax,  &
            xaa, xbb, yaa, ybb, zaa, zbb, t_stopped=0.

  integer*8 :: nax, nbx, nay, nby, naz, nbz
  integer*8, dimension(5) :: ppcx, ppcy, ppcz, nplx, nply, nplz
  integer :: n_sub_b, nspec, n_sort
  integer*8 :: nx, ny, nz 

  logical :: smoothing 
  integer*8 :: smooth_pass 
  
  ! field solver
  integer*8 :: ieta, netax, netay, eta_par, eta_zs, mask_zs 
  real*8 :: mask_r, mask_B0_fac
  logical :: mask
  real*8 :: etamin, etamax

  real*8 ::  hx, hy, hz, hxi, hyi, hzi, dtxi, dtyi, dtzi, &
            efld, bfld, efluid, ethermal, eptcl, time, te0

  integer*8 :: nsteps0, iwt=0, nx1, nx2, ny1, ny2, nz1, nz2, iopen, file_unit(25), &
              file_unit_read(20), nptot, npleaving, npentering, iclock_speed, nptotp

  real*8, dimension(:), allocatable :: dfac
  integer*8, dimension(:), allocatable :: nskip, ipleft, iprite, ipsendleft, ipsendrite, iprecv, &
            ipsendtop, ipsendbot, ipsendlefttop,ipsendleftbot,ipsendritetop,ipsendritebot,ipsend
  integer*8:: idum
  integer*8, dimension(:), allocatable:: idmap
  integer*8, dimension(:,:), allocatable:: idmap_yz
  integer*8 :: kb, ke, jb, je, nsendtotp, nrecvtotp, nsendtot, nrecvtot
  integer*8, dimension(:), allocatable :: idfft, kvec, jvec
  integer*8 :: ihstb, ihste, isendid(4), irecvid(4,4)
  real*8 :: xtmp1m, xtmp2m, xbox_l, xbox_r, ybox_l, ybox_r, zbox_l, zbox_r
  real*8, dimension(:,:), allocatable :: buf, buf2, buf_p1
  real*8, dimension(:,:,:), allocatable :: buf_particle
  integer, dimension(:), allocatable :: buftime
  integer, parameter :: nprobes=6, nbufsteps=100, tracking_width=14
  integer :: maxtags=100, maxtags_pe, ntot 
  logical :: tracking_binary, tracking_mpi
  
  ! wave loading (inside box)
  real*8 :: dB_B0, wave_cycles, sign_cos
  integer :: wave_upramp, wave_flat, wave_downramp

  ! wave injection (max. 4 waves)
  logical :: inj_waves_b, inj_waves_bv, inj_waves_e
  integer, dimension(4) :: inj_z_pos, inj_wave_pol, inj_wave_radius
  real*8, dimension(4) :: inj_dB_B0, inj_wave_cycles, inj_sign_cos, &
        inj_t_upramp, inj_t_flat, inj_t_downramp 

  integer :: seed_size
  integer, allocatable :: seed(:)

  ! diagnostics
  integer :: n_print, n_diag_mesh, n_diag_energy, n_diag_probe, n_diag_tracking, &
            n_diag_particle, n_write_restart
  integer :: probe_x
  ! integer :: n_debug_ez=100

  ! some constant parameters
  real*8, parameter :: zero=0.0d0, one=1.0d0, two=2.0d0, one_half=0.5d0, pi=acos(-1.)


  contains

  !---------------------------------------------------------------------
  ! read in input file/parameters
  !---------------------------------------------------------------------
  subroutine read_input
    integer :: i 

    namelist /input/ tmax, dtwci, restart, MPI_IO_format, & 
      ! simulation domain
      nx, ny, nz, xmax, ymax, zmax, ppcx, ppcy, ppcz, node_conf, periods, &  
      ! xaa, xbb, nax, nbx, yaa, ybb, nay, nby, zaa, zbb, naz, nbz, &
      uniform_load_logical, &
      ! field solver
      n_sub_b, eta_par, mask, mask_zs, mask_r, mask_B0_fac, & 
      dB_B0, wave_cycles, sign_cos, wave_upramp, wave_flat, wave_downramp, &  
      ! wave injection
      inj_waves_b, inj_waves_bv, inj_waves_e, inj_dB_B0, inj_wave_cycles, inj_sign_cos, inj_wave_pol, &
      inj_wave_radius, inj_z_pos, inj_t_upramp, inj_t_flat, inj_t_downramp, &
      ! plasma  
      nspec, n_sort, qspec, wspec, frac, denmin, & 
      wpiwci, beta_spec, beta_elec, &  
      ieta, resis, netax, netay, etamin, etamax, eta_zs, &
      anisot, gamma, smoothing, smooth_pass, &
      ! diagnostics
      n_print, n_diag_mesh, n_diag_energy, n_diag_probe, & 
      n_diag_tracking, n_write_restart, n_diag_particle, &  
      probe_x, &
      tracking_binary, tracking_mpi, &
      xbox_l, xbox_r, ybox_l, ybox_r, zbox_l, zbox_r

    ! print logo
    if (myid == 0) then
      print*
      print*, "********************************************************"
      print*, "           ||     ||    ======        ======            "              
      print*, "           ||     ||           ||    ||      \\         "         
      print*, "           ||=====||     ======||    ||       ||        "         
      print*, "           ||     ||           ||    ||      //         "         
      print*, "           ||     ||    ======        ======            "  
      print*, "********************************************************"        
    endif 

    ! read in input deck
    if (myid == 0) then
      print*
      print*
      print*, "Reading input file"
      print*, "-------------------------------------------------"
      open(5, file='input.f90', form='formatted', status='old')
      read(5, nml=input, iostat=ierr)
      if (ierr == 0) then 
        print*, 'Input deck loaded successfully.' 
      else
        call error_abort("error reading input deck!")
      endif 
    endif

    ! Broadcast input parameters (read in at rank 0) to all other ranks
    ! global sim. info
    call MPI_BCAST(tmax                   ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(dtwci                  ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(restart                ,1     ,MPI_LOGICAL          ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(MPI_IO_format          ,1     ,MPI_LOGICAL          ,0,MPI_COMM_WORLD,IERR)
    ! sim. domain
    call MPI_BCAST(nx                     ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ny                     ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(nz                     ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xmax                   ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ymax                   ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zmax                   ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ppcx                   ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ppcy                   ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ppcz                   ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(node_conf              ,2     ,MPI_INTEGER          ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(periods                ,2     ,MPI_LOGICAL          ,0,MPI_COMM_WORLD,IERR)
    ! call MPI_BCAST(xaa                    ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    ! call MPI_BCAST(xbb                    ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    ! call MPI_BCAST(nax                    ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    ! call MPI_BCAST(nbx                    ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    ! call MPI_BCAST(yaa                    ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    ! call MPI_BCAST(ybb                    ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    ! call MPI_BCAST(nay                    ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    ! call MPI_BCAST(nby                    ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    ! call MPI_BCAST(zaa                    ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    ! call MPI_BCAST(zbb                    ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    ! call MPI_BCAST(naz                    ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    ! call MPI_BCAST(nbz                    ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(uniform_load_logical   ,1     ,MPI_LOGICAL          ,0,MPI_COMM_WORLD,IERR)
    ! field solver
    call MPI_BCAST(n_sub_b                ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(eta_par                ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(mask                   ,1     ,MPI_LOGICAL          ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(mask_zs                ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(mask_r                 ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(mask_B0_fac            ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(dB_B0                  ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(wave_cycles            ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(sign_cos               ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(wave_upramp            ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(wave_flat              ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(wave_downramp          ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    ! wave injection
    call MPI_BCAST(inj_waves_b            ,1     ,MPI_LOGICAL          ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(inj_waves_bv           ,1     ,MPI_LOGICAL          ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(inj_waves_e            ,1     ,MPI_LOGICAL          ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(inj_dB_B0              ,4     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(inj_wave_cycles        ,4     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(inj_sign_cos           ,4     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(inj_wave_pol           ,4     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(inj_wave_radius        ,4     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(inj_z_pos              ,4     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(inj_t_upramp           ,4     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(inj_t_flat             ,4     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(inj_t_downramp         ,4     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    ! plasma setup
    call MPI_BCAST(nspec                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_sort                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(qspec                  ,5     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(wspec                  ,5     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(frac                   ,5     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(denmin                 ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(wpiwci                 ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(beta_spec              ,5     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(beta_elec              ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ieta                   ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(resis                  ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(netax                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(netay                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(etamin                 ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(etamax                 ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(eta_zs                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(anisot                 ,5     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(gamma                  ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(smoothing              ,1     ,MPI_LOGICAL          ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(smooth_pass            ,1     ,MPI_INTEGER          ,0,MPI_COMM_WORLD,IERR)
    ! diagnostic control
    call MPI_BCAST(n_print                ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_diag_mesh            ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_diag_energy          ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_diag_probe           ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(probe_x                ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_diag_tracking        ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_diag_particle        ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(n_write_restart        ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(tracking_binary        ,1     ,MPI_LOGICAL          ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(tracking_mpi           ,1     ,MPI_LOGICAL          ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xbox_l                 ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(xbox_r                 ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ybox_l                 ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(ybox_r                 ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zbox_l                 ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(zbox_r                 ,1     ,MPI_DOUBLE_PRECISION ,0,MPI_COMM_WORLD,IERR)

    ! The unit of dt is 1/wci in input file, 
    ! but now converted to 1/wpi inside the code
    dt = dtwci * wpiwci

    ! boundaries of the uniform region
    xaa = 0.; xbb = xmax; nax = 0; nbx = nx
    yaa = 0.; ybb = ymax; nay = 0; nby = ny
    zaa = 0.; zbb = zmax; naz = 0; nbz = nz

    ! steps of iteration 
    ! (could be modified later in 'init_restart' if restart=.true.)
    it = 0; itrestart = 0; 
    itstart = it; itfinish = tmax/dtwci

    ! set output directories
    data_directory = 'data/'
    restart_directory = 'restart/'
    restart_index_suffix(1) = '.1'
    restart_index_suffix(2) = '.2'

    ! specify decomposition along y, z; no decomposition along x 
    if (nz==1 .and. ny==1) then ! only nx>=1 and 1 rank will be used  
      ndim=0; dims(1)=1; dims(2)=1
    else if (nz == 1) then ! ny>1 and decomposition only occurs in y
      ndim=1; dims(1)=node_conf(1); dims(2)=1
    else ! ny>1, nz>1, and decomposition in both y and z
      ndim=2; dims(1)=node_conf(1); dims(2)=node_conf(2)
    endif

    ! npy(z) now means number of particles in each rank along y(z)
    ! npy=npy/dims(1); npz=npz/dims(2)

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
    ! write(6,'(a,i5,a,i5,i5,a,i5,i5,a,i5,i5)') &
    !     ' myid=', myid, ',  jb, je =', jb, je, &
    !     ',  kb, ke = ',kb, ke, ',  coords =', coords

    ! max number of cells. why adding 2?
    nxmax = nx + 2; nymax = ny + 2; nzmax = nz + 2
    ! local (at each rank) max number of cells
    nylmax = je - jb + 1 ; nzlmax = ke - kb + 1  
    ! local particle number along x, y, z
    nplx = ppcx*nx; nply = ppcy*nylmax; nplz = ppcz*nzlmax

    ! print decomposition information
    if (myid == 0) then
      write(6,*)
      write(6,*) "Total number of processors = ", nprocs
      do i = 1, ndim
        write(6,'(a,i5,a,i5)') " Dimension = ", i, ", Dims = ", dims(i)
      enddo
      ! write(6,*)
      ! write(6,'(a29,5i5)') ' Local paricle number in x = ', npx 
      ! write(6,'(a29,5i5)') ' Local paricle number in y = ', npy 
      ! write(6,'(a29,5i5)') ' Local paricle number in z = ', npz 

      ! write(6,*)
      ! write(6,'(a29,i5)') " Local cell number in x = ", nx
      ! write(6,'(a29,i5)') " Local cell number in y = ", nylmax
      ! write(6,'(a29,i5)') " Local cell number in z = ", nzlmax
    endif

    return

  end subroutine read_input


  !---------------------------------------------------------------------
  ! Set global parameters
  !---------------------------------------------------------------------
  subroutine init_arrays
    integer :: i, j, k

    if (myid==0) then
      print*
      print*
      print*, "Initializing global arrays"
      print*, "-------------------------------------------------"
    endif

    ! INQUIRE (IOLENGTH=iolength) output-items
    ! iolength is a scalar default INTEGER variable having
    ! a value that would result from the use of output-items
    ! in an unformatted output statement.  The value is used
    ! as a RECL=specifier in an OPEN statement that connects
    ! a file for unformatted direct access when there are
    ! input/output statements with the same list of output-items.
    single_prec=0.; double_prec=0.
    inquire(IOLENGTH=recl_for_single) single_prec
    inquire(IOLENGTH=recl_for_double) double_prec

    ! local total particle number
    nptotp = 0 
    do i = 1, nspec
      nptotp = nptotp + nplx(i)*nply(i)*nplz(i)
    enddo
    nplmax = 5* nptotp  ! pad storage requirement by a factor; why?
    ! if (myid==0) then
    !   write(6,"(a35,i10)") " total particle # per rank       = ", nptotp
    !   write(6,"(a35,i10)") " total particle # per rank (x5)  = ", nplmax
    ! endif
    ! print "(2(a,i8))", "myid = ", myid, ',  nptotp = ', nptotp

    ! number of tags used to track particles per species-rank
    ! maxtags was initialized as 100
    maxtags_pe = maxtags/nprocs/nspec
    if (maxtags_pe==0) then
        maxtags_pe = 1 
        maxtags = maxtags_pe * nprocs
    endif
    if (myid==0) then
      write(6,*)
      write(6,"(a35,2i10)") " maxtags_pe, maxtags = ", maxtags_pe, maxtags
    endif 

    ! allocate arrays that depend on 'nprocs', 'nplmax'
    allocate( zbglobal(0:nprocs-1), zeglobal(0:nprocs-1), ybglobal(0:nprocs-1), yeglobal(0:nprocs-1), &
              kbglobal(0:nprocs-1), keglobal(0:nprocs-1), jbglobal(0:nprocs-1), jeglobal(0:nprocs-1), &
              nsendp(0:nprocs-1), nrecvp(0:nprocs-1) )

    allocate( x(nplmax), y(nplmax), z(nplmax), vx(nplmax), vy(nplmax), vz(nplmax), &
              link(nplmax), porder(nplmax), qp(nplmax), ptag(nplmax) )

    allocate( ninj(nspec), ninj_global(nspec), nescape(nspec), nescape_global(nspec), &
              npart(nspec), npart_global(nspec) )
            
    allocate( nescape_xy(nspec), nescape_yx(nspec), &
              nescape_xz(nspec), nescape_zx(nspec), &
              nescape_yz(nspec), nescape_zy(nspec), &
              nescape_xy_global(nspec), nescape_yx_global(nspec), nescape_xz_global(nspec), &
              nescape_zx_global(nspec), nescape_yz_global(nspec), nescape_zy_global(nspec) )

    allocate( x0(nspec), x1(nspec), tx0(nspec), vpar(nspec), vper(nspec) )

    allocate( dfac(nspec),nskip(nspec),ipleft(nspec),iprite(nspec),ipsendleft(nspec),ipsendrite(nspec), &
              iprecv(nspec),ipsendtop(nspec),ipsendbot(nspec),ipsendlefttop(nspec),ipsendleftbot(nspec), &
              ipsendritetop(nspec),ipsendritebot(nspec),ipsend(nspec) )     

    allocate( idmap_yz(0:ny+1,0:nz+1), idmap(0:nzmax), idfft(nzmax), kvec(nzlmax), jvec(nylmax) )

    ! Use CART_SHIFT to determine processor to immediate left (NBRLEFT) 
    ! and right (NBRRITE) of processor MYID
    ! Since code is aperiodic in z, need to manually set the left boundary for processor 0 and right boundary for nprocs-1
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

    ! what's doing here
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

    if (mod(coords(1),2) == 0.and.mod(coords(2),2) == 0) then
      irecvid(1,2)=nbrrite
      irecvid(2,2)=-1
      irecvid(3,2)=nbrleft
      irecvid(4,2)=-1
    else if (mod(coords(1)+1,2) == 0.and.mod(coords(2)+1,2) == 0) then
      irecvid(1,2)=-1
      irecvid(2,2)=nbrtop
      irecvid(3,2)=-1
      irecvid(4,2)=nbrbot
    else if (mod(coords(1),2) == 0.and.mod(coords(2)+1,2) == 0) then
      irecvid(1,2)=nbrritetop
      irecvid(2,2)=nbrlefttop
      irecvid(3,2)=nbrleftbot
      irecvid(4,2)=nbrritebot
    endif
    
    if (mod(coords(1),2) == 0.and.mod(coords(2)+1,2) == 0) then
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

    ! VR: again, this is much simpler than the commented block
    do i = 0, nprocs-1
      do k = kbglobal(i), keglobal(i)
        do j = jbglobal(i), jeglobal(i)
          idmap_yz(j, k) = i
        enddo
      enddo
    enddo

    ! VR: fill in ghost cells in idmap     
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
    call MPI_TYPE_VECTOR(int(nylmax+2,4), int(nx+2,4), int(nx+2,4), MPI_DOUBLE_PRECISION, striderz, IERR)
    call MPI_TYPE_COMMIT(striderz, IERR)

    ! allocate arrays for local ranks
    allocate ( uniform_mesh(nxmax,jb-1:je+1,kb-1:ke+1), &
      ex    (nxmax,jb-1:je+1,kb-1:ke+1), ey     (nxmax,jb-1:je+1,kb-1:ke+1), ez     (nxmax,jb-1:je+1,kb-1:ke+1), &
      bx    (nxmax,jb-1:je+1,kb-1:ke+1), by     (nxmax,jb-1:je+1,kb-1:ke+1), bz     (nxmax,jb-1:je+1,kb-1:ke+1), &
      bx_av (nxmax,jb-1:je+1,kb-1:ke+1), by_av  (nxmax,jb-1:je+1,kb-1:ke+1), bz_av  (nxmax,jb-1:je+1,kb-1:ke+1), &
      fox   (nxmax,jb-1:je+1,kb-1:ke+1), foy    (nxmax,jb-1:je+1,kb-1:ke+1), foz    (nxmax,jb-1:je+1,kb-1:ke+1), &
      curlex(nxmax,jb-1:je+1,kb-1:ke+1), curley (nxmax,jb-1:je+1,kb-1:ke+1), curlez (nxmax,jb-1:je+1,kb-1:ke+1), &
      bxs   (nxmax,jb-1:je+1,kb-1:ke+1), bys    (nxmax,jb-1:je+1,kb-1:ke+1), bzs    (nxmax,jb-1:je+1,kb-1:ke+1), &
      den   (nxmax,jb-1:je+1,kb-1:ke+1), deno   (nxmax,jb-1:je+1,kb-1:ke+1), denh   (nxmax,jb-1:je+1,kb-1:ke+1), &
      dpedx (nxmax,jb-1:je+1,kb-1:ke+1), dpedy  (nxmax,jb-1:je+1,kb-1:ke+1), dpedz  (nxmax,jb-1:je+1,kb-1:ke+1), & 
      vix   (nxmax,jb-1:je+1,kb-1:ke+1), viy    (nxmax,jb-1:je+1,kb-1:ke+1), viz    (nxmax,jb-1:je+1,kb-1:ke+1), &  
      vixo  (nxmax,jb-1:je+1,kb-1:ke+1), viyo   (nxmax,jb-1:je+1,kb-1:ke+1), vizo   (nxmax,jb-1:je+1,kb-1:ke+1), & 
      curlbx(nxmax,jb-1:je+1,kb-1:ke+1), curlby (nxmax,jb-1:je+1,kb-1:ke+1), curlbz (nxmax,jb-1:je+1,kb-1:ke+1), & 
      pe    (nxmax,jb-1:je+1,kb-1:ke+1), eta    (nxmax,jb-1:je+1,kb-1:ke+1), eta_times_b_dot_j(nxmax,jb-1:je+1,kb-1:ke+1) )

    allocate ( dns(nxmax,jb-1:je+1,kb-1:ke+1,nspec), &
               dnsh(nxmax,jb-1:je+1,kb-1:ke+1,nspec), &
               vxs(nxmax,jb-1:je+1,kb-1:ke+1,nspec), &
               vys(nxmax,jb-1:je+1,kb-1:ke+1,nspec), &
               vzs(nxmax,jb-1:je+1,kb-1:ke+1,nspec), &
               tpar(nxmax,jb-1:je+1,kb-1:ke+1,nspec), &
               tperp(nxmax,jb-1:je+1,kb-1:ke+1,nspec), &
               qp_cell(nxmax,jb-1:je+1,kb-1:ke+1,nspec) ) 

    allocate ( p_xx(nxmax,jb-1:je+1,kb-1:ke+1,nspec), &
               p_xy(nxmax,jb-1:je+1,kb-1:ke+1,nspec), &
               p_xz(nxmax,jb-1:je+1,kb-1:ke+1,nspec), &
               p_yy(nxmax,jb-1:je+1,kb-1:ke+1,nspec), &
               p_yz(nxmax,jb-1:je+1,kb-1:ke+1,nspec), &
               p_zz(nxmax,jb-1:je+1,kb-1:ke+1,nspec) )

    allocate ( ainjxz(nxmax,kb-1:ke+1), ainjzx(nxmax,kb-1:ke+1), &
               deavxz(nxmax,kb-1:ke+1), deavzx(nxmax,kb-1:ke+1), &
               vxavxz(nxmax,kb-1:ke+1), vyavxz(nxmax,kb-1:ke+1), vzavxz(nxmax,kb-1:ke+1), &
               vxavzx(nxmax,kb-1:ke+1), vyavzx(nxmax,kb-1:ke+1), vzavzx(nxmax,kb-1:ke+1), &
               vxcaxz(nxmax,kb-1:ke+1), vycaxz(nxmax,kb-1:ke+1), vzcaxz(nxmax,kb-1:ke+1), &
              vxcazx(nxmax,kb-1:ke+1),vycazx(nxmax,kb-1:ke+1), vzcazx(nxmax,kb-1:ke+1) )

    allocate ( ainjyz(jb-1:je+1,kb-1:ke+1), ainjzy(jb-1:je+1,kb-1:ke+1), &
               deavyz(jb-1:je+1,kb-1:ke+1), deavzy(jb-1:je+1,kb-1:ke+1), &
               vxavyz(jb-1:je+1,kb-1:ke+1), vyavyz(jb-1:je+1,kb-1:ke+1), &
               vzavyz(jb-1:je+1,kb-1:ke+1), vxavzy(jb-1:je+1,kb-1:ke+1), &
               vyavzy(jb-1:je+1,kb-1:ke+1), vzavzy(jb-1:je+1,kb-1:ke+1), &
               vxcayz(jb-1:je+1,kb-1:ke+1), vycayz(jb-1:je+1,kb-1:ke+1), &
               vzcayz(jb-1:je+1,kb-1:ke+1), vxcazy(jb-1:je+1,kb-1:ke+1), &
               vycazy(jb-1:je+1,kb-1:ke+1), vzcazy(jb-1:je+1,kb-1:ke+1) )

    allocate ( ainjxy(nxmax,jb-1:je+1),ainjyx(nxmax,jb-1:je+1),deavxy(nxmax,jb-1:je+1), &
               deavyx(nxmax,jb-1:je+1),vxavxy(nxmax,jb-1:je+1),vyavxy(nxmax,jb-1:je+1), &
               vzavxy(nxmax,jb-1:je+1),vxavyx(nxmax,jb-1:je+1),vyavyx(nxmax,jb-1:je+1), &
               vzavyx(nxmax,jb-1:je+1),vxcaxy(nxmax,jb-1:je+1),vycaxy(nxmax,jb-1:je+1), &
               vzcaxy(nxmax,jb-1:je+1),vxcayx(nxmax,jb-1:je+1),vycayx(nxmax,jb-1:je+1), &
               vzcayx(nxmax,jb-1:je+1) )

    allocate ( iphead(nxmax,jb-1:je+1,kb-1:ke+1,nspec), &
               iptemp(nxmax,jb-1:je+1,kb-1:ke+1,nspec) )

    allocate ( xc_uniform(nxmax), yc_uniform(nymax), zc_uniform(nzmax), &
               xv_uniform(nxmax), yv_uniform(nymax), zv_uniform(nzmax) ) 

    allocate ( ixc_2_c_map(nx+1), iyc_2_c_map(ny+1), izc_2_c_map(nz+1) )
    allocate ( ixc_2_v_map(nx+1), iyc_2_v_map(ny+1), izc_2_v_map(nz+1) )
    allocate ( ixv_2_c_map(nx+1), iyv_2_c_map(ny+1), izv_2_c_map(nz+1) )
    allocate ( ixv_2_v_map(nx+1), iyv_2_v_map(ny+1), izv_2_v_map(nz+1) )
    allocate ( buf(nprobes,nbufsteps), buf2(nprobes*nprocs,nbufsteps), buftime(nbufsteps) )
    allocate ( buf_particle(tracking_width,nspec*maxtags,nbufsteps) )
    allocate ( buf_p1(tracking_width,nspec*maxtags) )

    ! make particle list
    if (myid == 0) then
      print*
      print*, 'done allocation of global arrays'
      print*, 'calling make particle list'
    endif 
    call makelist

  end subroutine init_arrays 


  !---------------------------------------------------------------------
  ! make list
  !---------------------------------------------------------------------
  subroutine makelist
    integer*8:: ip

    ipstore = 1 ! integer scalar
    ipleft = 0 ! ipleft(nspec)
    iprite = 0
    iprecv = 0
    iphead = 0 ! iphead(nxmax,jb-1:je+1,kb-1:ke+1,nspec)
    iptemp = 0

    do ip = 1, nplmax-1
        link(ip) = ip+1  
    enddo
    link(nplmax)=0

    return
  end subroutine makelist

end module m_parameter
