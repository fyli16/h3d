! Declare global parameters and global arrays in this module
module parameter_mod  
! implicit double precision (a-h,o-z),integer*8(i-n)
use mpi
implicit none
save
integer :: my_short_int, i_source, i_destination, i_tag, i_length, i_i
integer, dimension(8,128) :: time_begin_array,time_end_array
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
real*8, dimension(:,:,:,:), allocatable :: dns, dnsh, vxs, vys, vzs, tpar, tperp,qp_cell
real*8, dimension(:,:,:,:), allocatable :: p_xx,p_xy,p_xz,p_yy,p_yz,p_zz

real*8, dimension(:,:), allocatable ::ainjxz,ainjzx,deavxz,deavzx,vxavxz,vyavxz,vzavxz,vxavzx,      &
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
real*8, dimension(5) :: btspec,qspec,wspec,frac,anisot
double precision :: denmin, resis, wpiwci, bete, fxsho,ave1,ave2,phib,demin2, &
                    xmax,ymax,zmax,dt,gama,dtwci,quota,wall_clock_elapsed,tmax,buffer_zone,  &
                    xaa,xbb,yaa,ybb,zaa,zbb,t_stopped
integer*8 :: nax,nbx,nay,nby,naz,nbz
integer*8, dimension(8) :: wall_clock_begin,wall_clock_end
integer*8 :: eta_par, nparbuf
integer*8, dimension(5) :: npx, npy, npz
integer*8 :: iterb,norbskip,restrt_write,nxcel,netax,netay,netaz,nspec,nx,ny,nz,nprint, &
          nwrtdata, nwrtparticle, nwrtrestart, nskipx,nskipy,nskipz
real*8 :: etamin,etamax,moat_zone
integer*8 :: ieta, profile_power
logical :: testorbt, restart, uniform_loading_in_logical_grid, MPI_IO_format, smoothing
real*8 ::  hx,hy,hz,hxi,hyi,hzi,efld,bfld,efluidt,ethermt,eptclt,time,te0
logical :: prntinfo, wrtdat
integer :: it,notime
integer*8 :: nsteps0,itfin,iwt,nx1,nx2,ny1,ny2,nz1,nz2,iopen,file_unit(25),file_unit_time,            &
               file_unit_tmp,file_unit_read(20),nptot,npleaving,npentering,iclock_speed, nptotp
real*8 :: clock_time_init,clock_time_old,clock_time,clock_time1
real*8, dimension(:) ,allocatable:: dfac
integer*8, dimension(:) ,allocatable:: nskip,ipleft,iprite,ipsendleft,ipsendrite,iprecv,ipsendtop,ipsendbot     &
                                        ,ipsendlefttop,ipsendleftbot,ipsendritetop,ipsendritebot,ipsend
logical :: global, harris, Yee, post_process
integer*8:: idum
integer*8, dimension(:), allocatable:: idmap
integer*8, dimension(:,:), allocatable:: idmap_yz
integer*8 :: kb, ke, jb, je, nsendtotp, nrecvtotp, nsendtot, nrecvtot
integer*8, dimension(:), allocatable :: idfft, kvec, jvec, myid_stop
integer*8 :: ihstb, ihste, isendid(4), irecvid(4,4)
real*4 :: single_prec
real*8:: double_prec, xtmp1m, xtmp2m,initial_time,final_time    &
                    ,xbox_l,xbox_r,ybox_l,ybox_r,zbox_l,zbox_r,t_begin,t_end
real*8, dimension(:,:), allocatable :: buf, buf2, buf_p1
real*8, dimension(:,:,:), allocatable :: buf_particle
integer, dimension(:), allocatable :: buftime
integer, parameter :: nprobes=6, nbufsteps=100, tracking_width=14
integer :: maxtags=100, maxtags_pe, ntot 
logical :: tracking_binary=.true., tracking_mpi=.true.
integer :: tracking_fh

integer*8 :: recl_for_real, recl_for_double_precision
logical :: periods(2), reorder
integer*8 :: restart_index
character(len=2) :: restart_index_suffix(2)
character(len=160) :: data_directory, restart_directory, cycle_ascii_new, myid_char
character(len=160) :: cycle_ascii,cleanup_status
real*8 :: dB_B0, num_cycles ! for init_wave
real*8, parameter :: zero=0.0d0, one=1.0d0, two=2.0d0, one_half=0.5d0, pi=acos(-1.)

contains
! Set global parameters
subroutine set_parameters()
     implicit none
     integer :: i, j, k
     integer*8 :: nyl, nzl
     double_prec = 0.
     single_prec = 0.
     inquire (IOLENGTH=recl_for_double_precision) double_prec
     inquire (IOLENGTH=recl_for_real) single_prec

     nspecm = nspec  ! nspecm is just a mirror of nspec
     nxmax  = nx + 2  
     nymax  = ny + 2
     nzmax  = nz + 2
     nylmax = je - jb + 1  ! max of local array size in y
     nzlmax = ke - kb + 1  ! max of local array size in z
     nparbuf = nxmax*(nylmax+2)*(nzlmax+2)
     npes = numprocs  ! npes is a copy of numprocs
     npes_over_60 = npes / 512  ! if numprocs > 512

     do i = 1, nspecm
          qleft(i) = 0
          qrite(i) = 0
     enddo
     if (myid == 0) then
          write(6,*) "LOCAL ARRAY SIZE IN Y-DIRECTION = ", nylmax
          write(6,*) "LOCAL ARRAY SIZE IN Z-DIRECTION = ", nzlmax
     endif

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
     nzl = nzlmax
     nyl = nylmax

     ! estimate on particle storage requirement
     nplmax = 0  ! max number of local particles in each process
     do i = 1, nspec
          nplmax = nplmax + npx(i)*npy(i)*npz(i)  ! notice npy, npz are already divided by nodey, nodez respectively
     enddo
     nplmax = 2*nplmax  ! pad storage requirement by a factor of 2 

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

     ! gather jb,je,kb,ke of each rank into *global (where *=jb,je,kb,ke)
     call MPI_ALLGATHER(jb,1,MPI_INTEGER8,jbglobal,1,MPI_INTEGER8,MPI_COMM_WORLD,IERR)
     call MPI_ALLGATHER(je,1,MPI_INTEGER8,jeglobal,1,MPI_INTEGER8,MPI_COMM_WORLD,IERR)
     call MPI_ALLGATHER(kb,1,MPI_INTEGER8,kbglobal,1,MPI_INTEGER8,MPI_COMM_WORLD,IERR)
     call MPI_ALLGATHER(ke,1,MPI_INTEGER8,keglobal,1,MPI_INTEGER8,MPI_COMM_WORLD,IERR)

     !VR: again, this is much simpler than the commented block
     do i = 0, numprocs-1
          do k = kbglobal(i), keglobal(i)
             do j = jbglobal(i), jeglobal(i)
                idmap_yz(j,k) = i
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
 
     !VR: output neigbor info for each process & idmap
     ! write(tmpfname,"(A,I0,A)") "neighbors.",myid,".dat"
     ! open(unit=512,file=TRIM(tmpfname),status='replace',action='write')
     ! write(512,"(A,I6)") "TOP=",NBRTOP
     ! write(512,"(A,I6)") "TOPLEFT=",NBRLEFTTOP
     ! write(512,"(A,I6)") "LEFT=",NBRLEFT
     ! write(512,"(A,I6)") "BOTLEFT=",NBRLEFTBOT
     ! write(512,"(A,I6)") "BOT=",NBRBOT
     ! write(512,"(A,I6)") "BOTRITE=",NBRRITEBOT
     ! write(512,"(A,I6)") "RITE=",NBRRITE
     ! write(512,"(A,I6)") "RITETOP=",NBRRITETOP
     ! close(512)
     
     ! write(tmpfname,"(A,I0,A)") "idmap_yz.",myid,".dat"
     ! open(unit=512,file=TRIM(tmpfname),status='replace',action='write',access='stream',form='unformatted')
     ! write(512) idmap_yz
     ! close(512)

     ! print *, myid, "size of id_map is ",size(idmap_yz)
 
     call MPI_TYPE_VECTOR(int(nzl+2,4), int(nx+2,4), int((nx+2)*(nyl+2),4), MPI_DOUBLE_PRECISION, stridery, IERR)
     call MPI_TYPE_COMMIT(stridery, IERR)
     call MPI_TYPE_VECTOR(int(nyl+2,4), int(nx+2,4), int(nx+2,4), MPI_DOUBLE_PRECISION, STRIDERZ, IERR)
     call MPI_TYPE_COMMIT(STRIDERZ, IERR)
  
     nptotp = 0  ! total number of particles per processor
     do i=1, nspec
          nptotp = nptotp + npx(i)*npy(i)*npz(i)
     enddo
 
     do i = 0, npes-1
          i_i = i
          CALL MPI_BCAST(MYID_STOP(i),1,MPI_INTEGER8,i_i,MPI_COMM_WORLD, IERR)
     enddo
         
     if (.not.testorbt) norbskip=1
 
     call allocate_global_arrays

end subroutine set_parameters


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
end subroutine allocate_global_arrays

end module parameter_mod
