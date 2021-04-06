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
            xmax,ymax,zmax,dt,gama,dtwci,quota,wall_clock_elapsed,tmax,buffer_zone,  &
            xaa,xbb,yaa,ybb,zaa,zbb,t_stopped
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
  integer*8 :: nsteps0, itfin, iwt,nx1,nx2,ny1,ny2, nz1, nz2, iopen, file_unit(25), file_unit_time,            &
                file_unit_tmp,file_unit_read(20),nptot,npleaving,npentering,iclock_speed, nptotp
  real*8 :: clock_time_init, clock_time_old, clock_time, clock_time1
  real*8, dimension(:) ,allocatable:: dfac
  integer*8, dimension(:) ,allocatable:: nskip,ipleft,iprite,ipsendleft,ipsendrite,iprecv,ipsendtop,ipsendbot     &
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
  subroutine allocate_global_arrays
    implicit none
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
