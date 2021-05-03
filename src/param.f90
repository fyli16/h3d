! Declare global parameters and global arrays in this module
module m_parameter  
  use mpi
  implicit none

  save

  ! global simulation/MPI
  integer :: it, itstart, itfinish, itrestart, my_short_int, i_source, i_tag, i_length, i_i, &
             time_begin(8,128), time_end(8,128), ierr, n_subcycles
  real*8 :: time_elapsed(128)

  logical :: periods(2), reorder
  integer :: status(mpi_status_size), status1(mpi_status_size), status2(mpi_status_size), status_array(mpi_status_size,8)

  integer*8 :: nxmax, nymax, nzmax, nvar, nylmax, nzlmax, npm
  
  integer :: nprocs, ndim, dims(2), node_conf(2), comm2d, myid, req(8), & 
            nbrtop, nbrbot, nbrritetop, nbrlefttop, nbrritebot, nbrleftbot, &      
            nbrleft, nbrrite, ipe, stridery, striderz, iseed(1), coords(2)

  real*8 :: zb, ze, yb, ye, volume_fraction, cell_volume_ratio, &
            zb_logical, ze_logical, yb_logical, ye_logical, &
            xb_logical, xe_logical, xb, xe

  real*8, dimension(:), allocatable :: zbglobal, zeglobal, ybglobal, yeglobal, &
                                      xc_uniform, yc_uniform, zc_uniform, &
                                      xv_uniform, yv_uniform, zv_uniform

  integer*8, dimension(:), allocatable :: kbglobal, keglobal, jbglobal, jeglobal, &
      nsendp, nrecvp, ixc_2_c_map, iyc_2_c_map, izc_2_c_map, ixc_2_v_map, iyc_2_v_map, izc_2_v_map, &
      ixv_2_c_map, iyv_2_c_map, izv_2_c_map, ixv_2_v_map, iyv_2_v_map, izv_2_v_map  

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

  real*8, dimension(:), allocatable :: x, y, z, vx, vy, vz, qp

  integer*8, dimension(:), allocatable :: ptag, link, porder

  integer*8 :: nplmax, ipstore, np 

  integer*8, dimension(:,:,:,:), allocatable:: iphead, iptemp
  integer*8, dimension(:), allocatable ::  ninj, ninj_global, nescape,nescape_global, npart, npart_global
  integer*8, dimension(:),allocatable :: nescape_yz,nescape_zy,nescape_xy                 &
                                        ,nescape_yx,nescape_xz,nescape_zx                 &
                                        ,nescape_yz_global,nescape_zy_global              &
                                        ,nescape_xy_global,nescape_yx_global              &
                                        ,nescape_xz_global,nescape_zx_global

  real*8, dimension(:), allocatable:: x0,x1,tx0,vpar,vper

  real*8, dimension(5) :: beta_spec, qspec, wspec, frac, anisot

  real*8 :: denmin, resis, wpiwci, beta_e, ave1, ave2, phib, demin2, &
            xmax, ymax, zmax, dt, gamma, dtwci, wall_clock_elapsed, tmax,  &
            xaa, xbb, yaa, ybb, zaa, zbb, t_stopped=0.

  integer*8 :: nax, nbx, nay, nby, naz, nbz
  integer*8, dimension(8) :: wall_clock_begin,wall_clock_end
  integer*8, dimension(5) :: npx, npy, npz
  integer*8 :: iterb, nspec, n_sort, nx, ny, nz, n_print, &
            n_diag_mesh, n_diag_energy, n_diag_probe, n_diag_tracking, &
            n_diag_particle, n_write_restart, nskipx,nskipy,nskipz

  ! resistivity
  integer*8 :: ieta, netax, netay, eta_par, eta_zs
  real*8 :: etamin, etamax

  logical :: restart, uniform_load_logical, MPI_IO_format, smoothing 
  integer :: smooth_pass 

  real*8 ::  hx, hy, hz, hxi, hyi, hzi, efld, bfld, efluid, ethermal, eptcl, time, te0
  integer*8 :: nsteps0, iwt=0, nx1, nx2, ny1, ny2, nz1, nz2, iopen, file_unit(25), &
              file_unit_read(20), nptot, npleaving, npentering, iclock_speed, nptotp
  real*8 :: clock_init, clock_old, clock_now, clock_time1
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
  integer*8 :: restart_index=1
  character(len=2) :: restart_index_suffix(2)
  character(len=160) :: data_directory, restart_directory, cycle_ascii, cycle_ascii_new, &
                        myid_char, cleanup_status
  
  integer*8:: recl_for_single, recl_for_double
  real :: single_prec
  real*8 :: double_prec

  real*8 :: dB_B0, num_wave_cycles ! for initializing waves

  real*8, parameter :: zero=0.0d0, one=1.0d0, two=2.0d0, one_half=0.5d0, pi=acos(-1.)

  contains

  !---------------------------------------------------------------------
  ! Set global parameters
  !---------------------------------------------------------------------
  subroutine init_arrays
    integer :: i, j, k

    if (myid==0) then
      print*, " "
      print*, "Setting up global arrays"
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

    ! estimate on particle storage requirement
    nptotp = 0  ! total local number of particles
    do i = 1, nspec
      nptotp = nptotp + npx(i)*npy(i)*npz(i)
    enddo
    nplmax = 5* nptotp  ! pad storage requirement by a factor; why?
    if (myid==0) then
      print*, "total particle # per rank      = ", nptotp
      print*, "total particle # per rank (x5) = ", nplmax
    endif

    ! number of tags used to track particles per species-rank
    ! maxtags was initialized as 100
    maxtags_pe = maxtags/nprocs/nspec
    if (maxtags_pe==0) then
        maxtags_pe = 1 
        maxtags = maxtags_pe * nprocs
    endif
    if (myid==0) print*, "maxtags_pe, maxtags = ", maxtags_pe, maxtags

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
