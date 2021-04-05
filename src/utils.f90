!********************************************************
!
!********************************************************
subroutine accumulate_time_difference(time_begin, time_end, time_elapsed)
  implicit none
  integer,dimension(8):: time_begin,time_end
  double precision:: time_elapsed
  
  time_elapsed=time_elapsed &
       +(time_end(3)-time_begin(3))*3600.*24. &
       +(time_end(5)-time_begin(5))*3600. &
       +(time_end(6)-time_begin(6))*60. &
       +(time_end(7)-time_begin(7)) &
       +(time_end(8)-time_begin(8))*0.001
  return
end subroutine accumulate_time_difference


!********************************************************
!XF:  smoothing routine--for periodic B.C.
!3D version of 3-point binomial smoothing
!            y(i)=(x(i-1)+2*x(i)+x(i+1))/4
!i.e. 27 points are involved
!********************************************************
subroutine nsmth (a)
  use parameter_mod
  implicit none
  integer*8 i,j,k
  double precision, dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: temp, a

  ! copy input array "a" to "temp" including ghost cells
  do k=kb-1,ke+1
      do j = jb-1,je+1
        do i=1,nx2
            temp(i,j,k)=a(i,j,k)
        enddo
      enddo
  enddo

  ! smoothing only for inner cells (exclude ghost cells)
  do k=kb,ke
      do j = jb,je
        do i=2,nx1
            a(i,j,k)=temp(i,j,k)/8.&
                +( temp(i-1,j,k)+temp(i+1,j,k)+temp(i,j+1,k)+temp(i,j-1,k)&
                +temp(i,j,k+1)+temp(i,j,k-1))/16.&
                +( temp(i+1,j+1,k)+temp(i+1,j-1,k)+temp(i-1,j+1,k)&
                +temp(i-1,j-1,k)&
                +temp(i,j+1,k+1)+temp(i,j-1,k+1)+temp(i,j+1,k-1)+temp(i,j-1,k-1)&
                +temp(i+1,j,k+1)+temp(i-1,j,k+1)+temp(i+1,j,k-1)&
                +temp(i-1,j,k-1))/32.&
                +( temp(i+1,j+1,k+1)+temp(i-1,j+1,k+1)&
                +temp(i+1,j-1,k+1)+temp(i-1,j-1,k+1)&
                +temp(i+1,j+1,k-1)+temp(i-1,j+1,k-1)&
                +temp(i+1,j-1,k-1)+temp(i-1,j-1,k-1))/64.
        enddo
      enddo
  enddo

  ! apply periodic BCs 
  call XREALBCC(a,0_8,NX,NY,NZ)

  return
end subroutine nsmth


!********************************************************
!>    sort the particles
!********************************************************
subroutine sortit
  use parameter_mod
  use mesh2d
  implicit none
  double precision pstore(nplmax)
  integer pstore2(nplmax)
  integer*8 id, kb1, is, ix, iy, iz, ixe ,iye, ize, l, nttot, nplist
  double precision:: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi
  
  dtxi = 1./meshX%dt
  dtyi = 1./meshY%dt
  dtzi = 1./meshZ%dt

  id = 0
  kb1 = kb-1
  iptemp = 0
  porder = 0
  do is = 1,nspec
    do iz = kb-1,ke
      do iy = jb-1,je
        do ix = 1, nx1
          np = iphead(ix,iy,iz,is)
          DO WHILE (NP.NE.0)

          ! Uniform mesh - Same as in version 5.0
          !  ixe = hxi*x(np)+1.5000000000000001d+00
          !  iye = hyi*y(np)+0.5000000000000001d+00
          !  ize = hzi*z(np)+0.5000000000000001d+00

            ! Nonuniform mesh - using MESH_UNMAP
              rxe=dtxi*MESH_UNMAP(meshX,x(np))+1.50000000000d+00
              rye=dtyi*MESH_UNMAP(meshY,y(np))+1.50000000000d+00
              rze=dtzi*MESH_UNMAP(meshZ,z(np))+1.50000000000d+00
              ixe=rxe
              iye=rye
              ize=rze
              iye=iye-1             ! integer index in y direction starts at 0
              ize=ize-1             ! integer index in z direction starts at 0

            porder(np)=iptemp(ixe,iye,ize,is)
            iptemp(ixe,iye,ize,is)=np
            np = link(np)
          ENDDO
        enddo
      enddo
    enddo
  enddo

  l = 0
  nttot = 0
  do is = 1,nspec
    do iz = kb-1,ke
      do iy = jb-1,je
        do ix = 1, nx1
          np = iptemp(ix,iy,iz,is)
          nplist = 0
          DO WHILE (NP.NE.0)
            nplist = nplist+1
            l = l+1
            link(l) = np
            np = porder(np)
          ENDDO
          nttot = nttot + nplist
          iphead(ix,iy,iz,is) = nplist
        enddo
      enddo
    enddo
  enddo

  id = 0
  kb1 = kb-1
  do l = 1,nttot
    pstore(l) = vx(link(l))
  enddo
  do l = 1,nttot
    vx(l) = pstore(l)
  enddo

  do l = 1,nttot
    pstore(l) = vy(link(l))
  enddo
  do l = 1,nttot
    vy(l) = pstore(l)
  enddo      

  do l = 1,nttot
    pstore(l) = vz(link(l))
  enddo
  do l = 1,nttot
    vz(l) = pstore(l)
  enddo      

  do l = 1,nttot
    pstore(l) = x(link(l))
  enddo
  do l = 1,nttot
    x(l) = pstore(l)
  enddo      

  do l = 1,nttot
    pstore(l) = y(link(l))
  enddo
  do l = 1,nttot
    y(l) = pstore(l)
  enddo      

  do l = 1,nttot
    pstore(l) = z(link(l))
  enddo
  do l = 1,nttot
    z(l) = pstore(l)
  enddo

  do l = 1,nttot
    pstore(l) = qp(link(l))
  enddo
  do l = 1,nttot
    qp(l) = pstore(l)
  enddo

  do l = 1,nttot
    pstore2(l) = ptag(link(l))
  enddo
  do l = 1,nttot
    ptag(l) = pstore2(l)
  enddo

  l=1
  do is = 1,nspec
    do iz = kb-1,ke
      do iy = jb-1,je
        do ix = 1, nx1
    nplist = iphead(ix,iy,iz,is)
    if (nplist.ne.0) then
            iphead(ix,iy,iz,is) = l
            do np = l, l+nplist-1
              link(np) = np+1
            enddo
            link(l+nplist-1) = 0
            l = l + nplist
          else
            iphead(ix,iy,iz,is) = 0
          endif
        enddo
      enddo
    enddo
  enddo

  if (l-1.ne.nttot) then
    print *,'Problem in SORT: l-1 NE NTTOT'
    stop
  endif
  ipstore = nttot + 1
  do l = nttot+1, nplmax-1
    link(l) = l+1
  enddo
  link(nplmax) = 0

  return
end subroutine sortit


!********************************************************
!********************************************************
subroutine clock_write(iunit,message,i2,i1,is,it)
  implicit none
  integer*8 iunit,i2,i1,is,it
  character*10 message
  write(iunit,"(i4,a10,e20.8)") it, real(i2-i1)/real(is)
  return
end subroutine clock_write


!---------------------------------------------------------------------
subroutine makelist
  use parameter_mod
  implicit none
  integer*8:: ip

  ipstore=1
  ipleft =0
  iprite =0
  iprecv =0
  iphead =0
  iptemp =0
  do ip=1,nplmax-1
      link(ip)=ip+1
  enddo
  link(nplmax)=0
  return
end subroutine makelist


!***********************************************************************
!>    domain decompostion util: splits n elements between numprocs processors
!!    @param n number of elements
!!    @param numprocs number of processors
!!    @param s,e : start/end indices
!***********************************************************************
subroutine MPE_DECOMP1D( n, numprocs, myid, s, e )
  implicit none  
  integer   numprocs, myid
  integer*8 n
  integer*8 s, e
  integer*8 nlocal
  integer*8 deficit

  nlocal  = n / numprocs
  s       = myid * nlocal + 1
  deficit = mod(n, int(numprocs,8) )
  s       = s + min( int(myid,8) ,deficit)
  if (myid  <  deficit) then
      nlocal = nlocal + 1
  endif
  e = s + nlocal - 1
  if (e  >  n .or. myid  ==  numprocs-1) e = n
  return
end subroutine MPE_DECOMP1D


!************************************************************************
!************************************************************************
subroutine get_cleanup_status(maxchar)
  use parameter_mod
  implicit none
  integer maxchar
  logical fexists

  if (myid==0) then
      inquire(file=trim(adjustl(data_directory))//'.cleanup_status',exist=fexists)
      if (fexists) then
        open(unit=1,file=trim(adjustl(data_directory))//'.cleanup_status',status='old')
        read(1,*) cleanup_status
        close(unit=1)
      else
        cleanup_status='CLEANUP_STATUS=FALSE'
      endif
  endif
  call MPI_BCAST(cleanup_status,maxchar,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
end subroutine get_cleanup_status


! ========================================
! generate simple KEY based on date & time
! ========================================
subroutine get_sim_id(key)
  implicit none
  character (len=8)  :: date
  character (len=10) :: time
  character (len=5)  :: zone
  character (len=27) :: key
  integer :: values(8)

  call DATE_AND_TIME(DATE, TIME, ZONE, VALUES) 
  key = date//'_'//time//'UTC'//zone
end subroutine get_sim_id


! ========================================
! creat id file 
! ========================================
subroutine create_id_file(key)
  implicit none
  character (len=27) :: key
  integer wunit,ierr
  logical fexist

  inquire(file="sim_id.txt", exist=fexist)

  if (fexist) then
    open(newunit = wunit, file="sim_id.txt",status='replace',iostat=ierr)
  else
    open(newunit = wunit, file="sim_id.txt",status='new',iostat=ierr)
  endif

  if (ierr==0) then
    write(wunit,'(A)') key
    close(wunit)
  else
    write(*,*) "Cannot create sim_id file!"
  endif
end subroutine create_id_file


! ========================================
! debug
! ========================================
subroutine DEBUG(myid,ln)
  implicit none
  integer :: myid,ln
  write(6,*)"myid=",myid,"line number=",ln
end subroutine


!-----------------------------------------------------------------
!    computes velocities?
!    what is the difference between vxs and vix?
!-----------------------------------------------------------------
subroutine trans        
  use parameter_mod
  use MESH2D
  implicit none

  integer*8:: is,i,j,k,jbmin,jbmax,kbmin,kbmax
  integer*4:: time_begin(8),time_end(8)
  double precision:: dttmp,dns_tmp

  call date_and_time(values=time_begin_array(:,20))

  do is=1,nspec
    DO K=KB-1,KE+1
      do j=jb-1,je+1
        do i=1,nx2
          dns(i,j,k,is)=1.e-10
          dnsh(i,j,k,is)=1.e-10
          vxs(i,j,k,is)=0.
          vys(i,j,k,is)=0.
          vzs(i,j,k,is)=0.
          if (is == 1) then
            deno(i,j,k)=den(i,j,k)
            vixo(i,j,k)=vix(i,j,k)
            viyo(i,j,k)=viy(i,j,k)
            vizo(i,j,k)=viz(i,j,k)
            den(i,j,k)=0.
            denh(i,j,k)=0.
            vix(i,j,k)=0.
            viy(i,j,k)=0.
            viz(i,j,k)=0.
          endif
        enddo
      enddo
    enddo
  enddo


  if (notime == 0) then
    call date_and_time(values=time_end)
    clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
    write(file_unit_time,"(i4,' parmov   ',f15.3)") it,real(clock_time-clock_time_init)
  endif
  call date_and_time(values=time_begin_array(:,7))


  if (ndim /= 1) then
     call parmov
   else
     call parmov_2d
   endif

  call date_and_time(values=time_end_array(:,7))
  call accumulate_time_difference(time_begin_array(1,7),time_end_array(1,7),time_elapsed(7))

  if (testorbt) return

  if (.false.) then
  do is=1,nspec
    do k=kb-1,ke+1
      do j=jb-1,je+1
        do i=1,nx2

        ! Nonuniform mesh
          cell_volume_ratio = hx*hy*hz/(meshX%dxc(i)*meshY%dxc(j+1)*meshZ%dxc(k+1))

          dns_tmp=dnsh(i,j,k,is)

        ! Uniform mesh - Same as is in version 5.0
        !  if (dns_tmp <= denmin) dns_tmp=1.d+10

        ! Nonuniform mesh
        !  if (dns_tmp*cell_volume_ratio <= denmin) dns_tmp=1.d+10
          if (dns_tmp*cell_volume_ratio <= denmin) dns_tmp=denmin/cell_volume_ratio ! July 21, 2010

          vxs(i,j,k,is)=vxs(i,j,k,is)/dns_tmp
          vys(i,j,k,is)=vys(i,j,k,is)/dns_tmp
          vzs(i,j,k,is)=vzs(i,j,k,is)/dns_tmp

        ! Uniform mesh - Same as is in version 5.0
        !  dns(i,j,k,is)=dns(i,j,k,is)

        ! Nonuniform mesh
        !  dns(i,j,k,is)=dns(i,j,k,is)*cell_volume_ratio

        enddo
      enddo
    enddo
  enddo
  endif

  do is=1,nspec
    do k=kb-1,ke+1
      do j=jb-1,je+1
        do i=1,nx2 
        den(i,j,k)=den(i,j,k)+dns(i,j,k,is)*qspec(is) 
        denh(i,j,k)=denh(i,j,k)+dnsh(i,j,k,is)*qspec(is) 
        !vix(i,j,k)=vix(i,j,k)+qspec(is)*dnsh(i,j,k,is)*vxs(i,j,k,is) 
        !viy(i,j,k)=viy(i,j,k)+qspec(is)*dnsh(i,j,k,is)*vys(i,j,k,is) 
        !viz(i,j,k)=viz(i,j,k)+qspec(is)*dnsh(i,j,k,is)*vzs(i,j,k,is)
        vix(i,j,k)=vix(i,j,k)+qspec(is)*vxs(i,j,k,is) 
        viy(i,j,k)=viy(i,j,k)+qspec(is)*vys(i,j,k,is) 
        viz(i,j,k)=viz(i,j,k)+qspec(is)*vzs(i,j,k,is)
        enddo
      enddo
    enddo
  enddo

! Apply Boundary Conditions
  if (ndim /= 1) then
  call XREALBCC(DEN,1_8,NX,NY,NZ)
  call XREALBCC(DENH,1_8,NX,NY,NZ)
  call XREALBCC(VIX,1_8,NX,NY,NZ)
  call XREALBCC(VIY,1_8,NX,NY,NZ)
  call XREALBCC(VIZ,1_8,NX,NY,NZ)
  else
  call XREALBCC_2D(DEN,1_8,NX,NY,NZ)
  call XREALBCC_2D(DENH,1_8,NX,NY,NZ)
  call XREALBCC_2D(VIX,1_8,NX,NY,NZ)
  call XREALBCC_2D(VIY,1_8,NX,NY,NZ)
  call XREALBCC_2D(VIZ,1_8,NX,NY,NZ)
  endif

! smooth density and velocity
  if (smoothing) then
    if (ndim /=1) then
      call nsmth(DEN)
      call nsmth(DENH)
      call nsmth(VIX)
      call nsmth(VIY)
      call nsmth(VIZ)
    else
      call nsmth_2d(DEN,NX2,NY2,NZ2)
      call nsmth_2d(DENH,NX2,NY2,NZ2)
      call nsmth_2d(VIX,NX2,NY2,NZ2)
      call nsmth_2d(VIY,NX2,NY2,NZ2)
      call nsmth_2d(VIZ,NX2,NY2,NZ2)
    endif
  endif

  do k=kb-1,ke+1
      do j=jb-1,je+1
      do i=1,nx2
        den(i,j,k)=max(denmin,den(i,j,k))
        !pressure calculation moved to pressgrad
        !pe(i,j,k) =te0*den(i,j,k)**gama
        vix(i,j,k)=vix(i,j,k)/denh(i,j,k)
        viy(i,j,k)=viy(i,j,k)/denh(i,j,k)
        viz(i,j,k)=viz(i,j,k)/denh(i,j,k)
      enddo
    enddo
  enddo

  if (it == 0) then
     deno=den;vixo=vix;viyo=viy;vizo=viz
  endif

  call date_and_time(values=time_begin_array(:,8))

  if (mod(it,10)==0) call energy

  call date_and_time(values=time_end_array(:,8))
  call accumulate_time_difference(time_begin_array(1,8),time_end_array(1,8),time_elapsed(8))

  kbmin = kb-1
  kbmax = ke+1

  jbmin = jb-1
  jbmax = je+1

  call date_and_time(values=time_end_array(:,20))
  call accumulate_time_difference(time_begin_array(1,20),time_end_array(1,20),time_elapsed(20))

  return
end subroutine trans


!-----------------------------------------------------------------
!     compute perp and par temperature and pressure tensor
!-----------------------------------------------------------------
  subroutine caltemp2_global

  use parameter_mod
  use MESH2D
  implicit none

  double precision:: rx,ry,rz,fx,fy,fz,dtxi,dtyi,dtzi,xx,xy,xz,yy,yz,zz
  integer*8 ix,iy,iz,ixp1,iyp1,izp1,iiy,iiye,iiz,iize,is,l,iix,iixe
  double precision vxa,vya,vza,rfrac,vxavg,vxavg1,vxavg2 &
       ,vyavg,vyavg1,vyavg2,vzavg,vzavg1,vzavg2,wperp2,wpar,wmult
  double precision w1,w2,w3,w4,w5,w6,w7,w8,h,hh,dns1,dns2,bxa,bya,bza,btota,dnst

  call date_and_time(values=time_begin_array(:,23))

  dtxi = 1./meshX%dt
  dtyi = 1./meshY%dt
  dtzi = 1./meshZ%dt


  tpar=0.
  tperp=0.

  p_xx=0.;p_xy=0.;p_xz=0.;p_yy=0.;p_yz=0.;p_zz=0.

  if (nspec >= 2) then
     rfrac = frac(2)/frac(1)
  else
     rfrac = 0.
  endif

  call date_and_time(values=time_begin_array(:,26))
  DO IS=1,NSPEC
    wmult=wspec(is)
    h=dt*qspec(is)/wmult
    hh=.5*h
    dpedx = 0.
    DO IIZE = KB-1,KE
      DO IIYE = JB-1,JE
        DO IIXE = 1, NX1
          NP=IPHEAD(IIXE,IIYE,IIZE,IS)

!  begin advance of particle position and velocity
!  If dt=0, skip
!
          DO WHILE (NP.NE.0)
            L=NP

          ! Uniform mesh - Same as is in version 5.0
          !  rx=hxi*x(l)+1.5000000000000001
          !  ry=hyi*y(l)+0.5000000000000001d+00
          !  rz=hzi*z(l)+0.5000000000000001d+00
          !  ix=rx
          !  iy=ry
          !  iz=rz
          !  fx=rx-ix
          !  fy=ry-iy
          !  fz=rz-iz

          ! Nonuniform mesh - using MESH_UNMAP
            rx=dtxi*MESH_UNMAP(meshX,x(l))+1.50000000000d+00
            ry=dtyi*MESH_UNMAP(meshY,y(l))+1.50000000000d+00
            rz=dtzi*MESH_UNMAP(meshZ,z(l))+1.50000000000d+00
            ix=rx
            iy=ry
            iz=rz
            fx=rx-ix
            fy=ry-iy
            fz=rz-iz
            iy=iy-1             ! integer index in y direction starts at 0
            iz=iz-1             ! integer index in z direction starts at 0

            ixp1 = ix+1
            iyp1 = iy+1
            izp1 = iz+1

            w1=(1.-fx)*(1.-fy)*(1.-fz)
            w2=    fx *(1.-fy)*(1.-fz)
            w3=(1.-fx)*    fy *(1.-fz)
            w4=    fx*     fy *(1.-fz)
            w5=(1.-fx)*(1.-fy)*    fz
            w6=    fx *(1.-fy)*    fz
            w7=(1.-fx)*    fy*     fz
            w8=    fx*     fy*     fz

            dnst= dnsh(ix  ,iy  ,iz  ,is)*w1+dnsh(ixp1,iy  ,iz  ,is)*w2  &
            +     dnsh(ix  ,iyp1,iz  ,is)*w3+dnsh(ixp1,iyp1,iz  ,is)*w4  &
            +     dnsh(ix  ,iy  ,izp1,is)*w5+dnsh(ixp1,iy  ,izp1,is)*w6  &
            +     dnsh(ix  ,iyp1,izp1,is)*w7+dnsh(ixp1,iyp1,izp1,is)*w8

            vxavg=vxs(ix  ,iy  ,iz  ,is)*w1+vxs(ixp1,iy  ,iz  ,is)*w2  &
            +      vxs(ix  ,iyp1,iz  ,is)*w3+vxs(ixp1,iyp1,iz  ,is)*w4  &
            +      vxs(ix  ,iy  ,izp1,is)*w5+vxs(ixp1,iy  ,izp1,is)*w6  &
            +      vxs(ix  ,iyp1,izp1,is)*w7+vxs(ixp1,iyp1,izp1,is)*w8
            vxavg=vxavg/dnst
            

            vyavg=vys(ix  ,iy  ,iz  ,is)*w1+vys(ixp1,iy  ,iz  ,is)*w2  &
            +      vys(ix  ,iyp1,iz  ,is)*w3+vys(ixp1,iyp1,iz  ,is)*w4  &
            +      vys(ix  ,iy  ,izp1,is)*w5+vys(ixp1,iy  ,izp1,is)*w6  &
            +      vys(ix  ,iyp1,izp1,is)*w7+vys(ixp1,iyp1,izp1,is)*w8
            vyavg=vyavg/dnst
            

            vzavg=vzs(ix  ,iy  ,iz  ,is)*w1+vzs(ixp1,iy  ,iz  ,is)*w2  &
            +      vzs(ix  ,iyp1,iz  ,is)*w3+vzs(ixp1,iyp1,iz  ,is)*w4  &
            +      vzs(ix  ,iy  ,izp1,is)*w5+vzs(ixp1,iy  ,izp1,is)*w6  &
            +      vzs(ix  ,iyp1,izp1,is)*w7+vzs(ixp1,iyp1,izp1,is)*w8
            vzavg=vzavg/dnst


            vxa=vx(l)-vxavg
            vya=vy(l)-vyavg
            vza=vz(l)-vzavg

            bxa  =bx(ix  ,iy  ,iz  )*w1+bx(ixp1,iy  ,iz  )*w2  &
            +     bx(ix  ,iyp1,iz  )*w3+bx(ixp1,iyp1,iz  )*w4  &
            +     bx(ix  ,iy  ,izp1)*w5+bx(ixp1,iy  ,izp1)*w6  &
            +     bx(ix  ,iyp1,izp1)*w7+bx(ixp1,iyp1,izp1)*w8

            bya  =by(ix  ,iy  ,iz  )*w1+by(ixp1,iy  ,iz  )*w2  &
            +     by(ix  ,iyp1,iz  )*w3+by(ixp1,iyp1,iz  )*w4  &
            +     by(ix  ,iy  ,izp1)*w5+by(ixp1,iy  ,izp1)*w6  &
            +     by(ix  ,iyp1,izp1)*w7+by(ixp1,iyp1,izp1)*w8

            bza  =bz(ix  ,iy  ,iz  )*w1+bz(ixp1,iy  ,iz  )*w2  &
            +     bz(ix  ,iyp1,iz  )*w3+bz(ixp1,iyp1,iz  )*w4  &
            +     bz(ix  ,iy  ,izp1)*w5+bz(ixp1,iy  ,izp1)*w6  &
            +     bz(ix  ,iyp1,izp1)*w7+bz(ixp1,iyp1,izp1)*w8

            btota=sqrt(bxa**2+bya**2+bza**2)
            if (btota < 1.e-20) btota=1.e-20
            wpar=(vxa*bxa+vya*bya+vza*bza)/btota
            wperp2=vxa**2+vya**2+vza**2-wpar**2
            xx=vxa*vxa
            xy=vxa*vya
            xz=vxa*vza
            yy=vya*vya
            yz=vya*vza
            zz=vza*vza

            tpar (ix  ,iy  ,iz  ,is)=tpar (ix  ,iy  ,iz  ,is)+qp(np)*w1*wpar*wpar
            tpar (ixp1,iy  ,iz  ,is)=tpar (ixp1,iy  ,iz  ,is)+qp(np)*w2*wpar*wpar 
            tpar (ix  ,iyp1,iz  ,is)=tpar (ix  ,iyp1,iz  ,is)+qp(np)*w3*wpar*wpar 
            tpar (ixp1,iyp1,iz  ,is)=tpar (ixp1,iyp1,iz  ,is)+qp(np)*w4*wpar*wpar 
            tpar (ix  ,iy  ,izp1,is)=tpar (ix  ,iy  ,izp1,is)+qp(np)*w5*wpar*wpar 
            tpar (ixp1,iy  ,izp1,is)=tpar (ixp1,iy  ,izp1,is)+qp(np)*w6*wpar*wpar 
            tpar (ix  ,iyp1,izp1,is)=tpar (ix  ,iyp1,izp1,is)+qp(np)*w7*wpar*wpar 
            tpar (ixp1,iyp1,izp1,is)=tpar (ixp1,iyp1,izp1,is)+qp(np)*w8*wpar*wpar 
            tperp(ix  ,iy  ,iz  ,is)=tperp(ix  ,iy  ,iz  ,is)+qp(np)*w1*wperp2 
            tperp(ixp1,iy  ,iz  ,is)=tperp(ixp1,iy  ,iz  ,is)+qp(np)*w2*wperp2 
            tperp(ix  ,iyp1,iz  ,is)=tperp(ix  ,iyp1,iz  ,is)+qp(np)*w3*wperp2 
            tperp(ixp1,iyp1,iz  ,is)=tperp(ixp1,iyp1,iz  ,is)+qp(np)*w4*wperp2
            tperp(ix  ,iy  ,izp1,is)=tperp(ix  ,iy  ,izp1,is)+qp(np)*w5*wperp2 
            tperp(ixp1,iy  ,izp1,is)=tperp(ixp1,iy  ,izp1,is)+qp(np)*w6*wperp2
            tperp(ix  ,iyp1,izp1,is)=tperp(ix  ,iyp1,izp1,is)+qp(np)*w7*wperp2 
            tperp(ixp1,iyp1,izp1,is)=tperp(ixp1,iyp1,izp1,is)+qp(np)*w8*wperp2
            dpedx(ix  ,iy  ,iz  )=dpedx(ix  ,iy  ,iz  )+qp(np)*w1
            dpedx(ixp1,iy  ,iz  )=dpedx(ixp1,iy  ,iz  )+qp(np)*w2 
            dpedx(ix  ,iyp1,iz  )=dpedx(ix  ,iyp1,iz  )+qp(np)*w3 
            dpedx(ixp1,iyp1,iz  )=dpedx(ixp1,iyp1,iz  )+qp(np)*w4 
            dpedx(ix  ,iy  ,izp1)=dpedx(ix  ,iy  ,izp1)+qp(np)*w5 
            dpedx(ixp1,iy  ,izp1)=dpedx(ixp1,iy  ,izp1)+qp(np)*w6 
            dpedx(ix  ,iyp1,izp1)=dpedx(ix  ,iyp1,izp1)+qp(np)*w7 
            dpedx(ixp1,iyp1,izp1)=dpedx(ixp1,iyp1,izp1)+qp(np)*w8 

            p_xx (ix  ,iy  ,iz  ,is)=p_xx (ix  ,iy  ,iz  ,is)+qp(np)*w1*xx
            p_xx (ixp1,iy  ,iz  ,is)=p_xx (ixp1,iy  ,iz  ,is)+qp(np)*w2*xx 
            p_xx (ix  ,iyp1,iz  ,is)=p_xx (ix  ,iyp1,iz  ,is)+qp(np)*w3*xx 
            p_xx (ixp1,iyp1,iz  ,is)=p_xx (ixp1,iyp1,iz  ,is)+qp(np)*w4*xx 
            p_xx (ix  ,iy  ,izp1,is)=p_xx (ix  ,iy  ,izp1,is)+qp(np)*w5*xx 
            p_xx (ixp1,iy  ,izp1,is)=p_xx (ixp1,iy  ,izp1,is)+qp(np)*w6*xx 
            p_xx (ix  ,iyp1,izp1,is)=p_xx (ix  ,iyp1,izp1,is)+qp(np)*w7*xx 
            p_xx (ixp1,iyp1,izp1,is)=p_xx (ixp1,iyp1,izp1,is)+qp(np)*w8*xx 

            p_xy (ix  ,iy  ,iz  ,is)=p_xy (ix  ,iy  ,iz  ,is)+qp(np)*w1*xy
            p_xy (ixp1,iy  ,iz  ,is)=p_xy (ixp1,iy  ,iz  ,is)+qp(np)*w2*xy 
            p_xy (ix  ,iyp1,iz  ,is)=p_xy (ix  ,iyp1,iz  ,is)+qp(np)*w3*xy 
            p_xy (ixp1,iyp1,iz  ,is)=p_xy (ixp1,iyp1,iz  ,is)+qp(np)*w4*xy 
            p_xy (ix  ,iy  ,izp1,is)=p_xy (ix  ,iy  ,izp1,is)+qp(np)*w5*xy 
            p_xy (ixp1,iy  ,izp1,is)=p_xy (ixp1,iy  ,izp1,is)+qp(np)*w6*xy 
            p_xy (ix  ,iyp1,izp1,is)=p_xy (ix  ,iyp1,izp1,is)+qp(np)*w7*xy 
            p_xy (ixp1,iyp1,izp1,is)=p_xy (ixp1,iyp1,izp1,is)+qp(np)*w8*xy 

            p_xz (ix  ,iy  ,iz  ,is)=p_xz (ix  ,iy  ,iz  ,is)+qp(np)*w1*xz
            p_xz (ixp1,iy  ,iz  ,is)=p_xz (ixp1,iy  ,iz  ,is)+qp(np)*w2*xz 
            p_xz (ix  ,iyp1,iz  ,is)=p_xz (ix  ,iyp1,iz  ,is)+qp(np)*w3*xz 
            p_xz (ixp1,iyp1,iz  ,is)=p_xz (ixp1,iyp1,iz  ,is)+qp(np)*w4*xz 
            p_xz (ix  ,iy  ,izp1,is)=p_xz (ix  ,iy  ,izp1,is)+qp(np)*w5*xz 
            p_xz (ixp1,iy  ,izp1,is)=p_xz (ixp1,iy  ,izp1,is)+qp(np)*w6*xz 
            p_xz (ix  ,iyp1,izp1,is)=p_xz (ix  ,iyp1,izp1,is)+qp(np)*w7*xz 
            p_xz (ixp1,iyp1,izp1,is)=p_xz (ixp1,iyp1,izp1,is)+qp(np)*w8*xz 

            p_yy (ix  ,iy  ,iz  ,is)=p_yy (ix  ,iy  ,iz  ,is)+qp(np)*w1*yy
            p_yy (ixp1,iy  ,iz  ,is)=p_yy (ixp1,iy  ,iz  ,is)+qp(np)*w2*yy 
            p_yy (ix  ,iyp1,iz  ,is)=p_yy (ix  ,iyp1,iz  ,is)+qp(np)*w3*yy 
            p_yy (ixp1,iyp1,iz  ,is)=p_yy (ixp1,iyp1,iz  ,is)+qp(np)*w4*yy 
            p_yy (ix  ,iy  ,izp1,is)=p_yy (ix  ,iy  ,izp1,is)+qp(np)*w5*yy 
            p_yy (ixp1,iy  ,izp1,is)=p_yy (ixp1,iy  ,izp1,is)+qp(np)*w6*yy 
            p_yy (ix  ,iyp1,izp1,is)=p_yy (ix  ,iyp1,izp1,is)+qp(np)*w7*yy 
            p_yy (ixp1,iyp1,izp1,is)=p_yy (ixp1,iyp1,izp1,is)+qp(np)*w8*yy 

            p_yz (ix  ,iy  ,iz  ,is)=p_yz (ix  ,iy  ,iz  ,is)+qp(np)*w1*yz
            p_yz (ixp1,iy  ,iz  ,is)=p_yz (ixp1,iy  ,iz  ,is)+qp(np)*w2*yz 
            p_yz (ix  ,iyp1,iz  ,is)=p_yz (ix  ,iyp1,iz  ,is)+qp(np)*w3*yz 
            p_yz (ixp1,iyp1,iz  ,is)=p_yz (ixp1,iyp1,iz  ,is)+qp(np)*w4*yz 
            p_yz (ix  ,iy  ,izp1,is)=p_yz (ix  ,iy  ,izp1,is)+qp(np)*w5*yz 
            p_yz (ixp1,iy  ,izp1,is)=p_yz (ixp1,iy  ,izp1,is)+qp(np)*w6*yz 
            p_yz (ix  ,iyp1,izp1,is)=p_yz (ix  ,iyp1,izp1,is)+qp(np)*w7*yz 
            p_yz (ixp1,iyp1,izp1,is)=p_yz (ixp1,iyp1,izp1,is)+qp(np)*w8*yz 

            p_zz (ix  ,iy  ,iz  ,is)=p_zz (ix  ,iy  ,iz  ,is)+qp(np)*w1*zz
            p_zz (ixp1,iy  ,iz  ,is)=p_zz (ixp1,iy  ,iz  ,is)+qp(np)*w2*zz 
            p_zz (ix  ,iyp1,iz  ,is)=p_zz (ix  ,iyp1,iz  ,is)+qp(np)*w3*zz 
            p_zz (ixp1,iyp1,iz  ,is)=p_zz (ixp1,iyp1,iz  ,is)+qp(np)*w4*zz 
            p_zz (ix  ,iy  ,izp1,is)=p_zz (ix  ,iy  ,izp1,is)+qp(np)*w5*zz 
            p_zz (ixp1,iy  ,izp1,is)=p_zz (ixp1,iy  ,izp1,is)+qp(np)*w6*zz 
            p_zz (ix  ,iyp1,izp1,is)=p_zz (ix  ,iyp1,izp1,is)+qp(np)*w7*zz 
            p_zz (ixp1,iyp1,izp1,is)=p_zz (ixp1,iyp1,izp1,is)+qp(np)*w8*zz 

            np=link(np)
          ENDDO
        ENDDO
      ENDDO
    ENDDO


    call XREAL(tpar (1,jb-1,kb-1,is),NX,NY,NZ)
    call XREAL(tperp(1,jb-1,kb-1,is),NX,NY,NZ)
    call XREAL(dpedx(1,jb-1,kb-1   ),NX,NY,NZ)

    call XREAL(p_xx (1,jb-1,kb-1,is),NX,NY,NZ)
    call XREAL(p_xy (1,jb-1,kb-1,is),NX,NY,NZ)
    call XREAL(p_xz (1,jb-1,kb-1,is),NX,NY,NZ)
    call XREAL(p_yy (1,jb-1,kb-1,is),NX,NY,NZ)
    call XREAL(p_yz (1,jb-1,kb-1,is),NX,NY,NZ)
    call XREAL(p_zz (1,jb-1,kb-1,is),NX,NY,NZ)

    DO IZ = KB-1,KE
      DO IY = JB-1,JE
        DO IX = 1, NX1
          if (dpedx(ix,iy,iz) /= 0.) then
            tpar (ix,iy,iz,is) = tpar (ix,iy,iz,is)/(   tx0(is)*dpedx(ix,iy,iz))
            tperp(ix,iy,iz,is) = tperp(ix,iy,iz,is)/(2.*tx0(is)*dpedx(ix,iy,iz))
          endif
        ENDDO
      ENDDO
    ENDDO

    DO IIZ=KB-1,KE+1
      DO IIY=JB-1,JE+1
        DO IIX=1,NX2
          p_xx(iix,iiy,iiz,is) = p_xx(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_xy(iix,iiy,iiz,is) = p_xy(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_xz(iix,iiy,iiz,is) = p_xz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_yy(iix,iiy,iiz,is) = p_yy(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_yz(iix,iiy,iiz,is) = p_yz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_zz(iix,iiy,iiz,is) = p_zz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
        ENDDO
      ENDDO
    ENDDO

  !  p_xx(:,:,:,is)=p_xx(:,:,:,is)/(tx0(is)*frac(is))
  !  p_xy(:,:,:,is)=p_xy(:,:,:,is)/(tx0(is)*frac(is))
  !  p_xz(:,:,:,is)=p_xz(:,:,:,is)/(tx0(is)*frac(is))
  !  p_yy(:,:,:,is)=p_yy(:,:,:,is)/(tx0(is)*frac(is))
  !  p_yz(:,:,:,is)=p_yz(:,:,:,is)/(tx0(is)*frac(is))
  !  p_zz(:,:,:,is)=p_zz(:,:,:,is)/(tx0(is)*frac(is))

  ENDDO
  call date_and_time(values=time_end_array(:,26))
  call accumulate_time_difference(time_begin_array(1,26) &
 &                               ,time_end_array(1,26) &
 &                                ,time_elapsed(26))


!  do is=1,nspec
!    call date_and_time(values=time_begin_array(:,24))
!    call XREAL(tpar (1,jb-1,kb-1,is),NX,NY,NZ)
!    call XREAL(tperp(1,jb-1,kb-1,is),NX,NY,NZ)
!    call date_and_time(values=time_end_array(:,24))
!    call accumulate_time_difference(time_begin_array(1,24) &
! &                                 ,time_end_array(1,24) &
! &                                 ,time_elapsed(24))

!    call date_and_time(values=time_begin_array(:,25))
!    call XREALBCC(tpar (1,jb-1,kb-1,is),1,NX,NY,NZ)
!    call XREALBCC(tperp(1,jb-1,kb-1,is),1,NX,NY,NZ)
!    call date_and_time(values=time_end_array(:,25))
!    call accumulate_time_difference(time_begin_array(1,25) &
! &                                 ,time_end_array(1,25) &
! &                                 ,time_elapsed(25))

!  enddo


!  call date_and_time(values=time_begin_array(:,26))
!  do is=1,nspec
!    do k=kb-1,ke+1
!      do j = jb-1,je+1
!        do i=1,nx2
!          if (is == 1) then
!            dns1=dns(i,j,k,1)/(dfac(1)*frac(1))
!            dns2=0.
!            denum=dns1+rfrac*dns2
!          else
!            denum=dns(i,j,k,is)/(dfac(is)*frac(is))
!          endif
!          if (denum < denmin)  then
!           tpar(i,j,k,is)=1.e-5
!           tperp(i,j,k,is)=1.e-5
!          else
!           denum=denum*tx0(is)
!           tpar(i,j,k,is)=tpar(i,j,k,is)*wspec(is)/denum
!           tperp(i,j,k,is)=0.5*tperp(i,j,k,is)*wspec(is)/denum
!          endif
!        enddo
!      enddo
!    enddo
!  enddo
!  call date_and_time(values=time_end_array(:,26))
!  call accumulate_time_difference(time_begin_array(1,26) &
! &                               ,time_end_array(1,26) &
! &                               ,time_elapsed(26))

   call date_and_time(values=time_end_array(:,23))
   call accumulate_time_difference(time_begin_array(1,23) &
 &                                ,time_end_array(1,23) &
 &                                ,time_elapsed(23))

  return
end subroutine caltemp2_global


!********************************************************
!>    computes field energy ex^2+ey^2+ez^2 and bx^2+by^2+bz^2
!!    and particle energies
!********************************************************
subroutine energy
  use parameter_mod
  use MESH2D
  implicit none
  double precision:: rx,ry,rz,fx,fy,fz,dtxi,dtyi,dtzi,xx,xy,xz,yy,yz,zz
  integer*8 ix,iy,iz,ixp1,iyp1,izp1,iiy,iiye,iiz,iize,is,l,iix,iixe
  double precision vxa,vya,vza,rfrac,vxavg,vxavg1,vxavg2 &
        ,vyavg,vyavg1,vyavg2,vzavg,vzavg1,vzavg2,wperp2,wpar,wmult
  double precision w1,w2,w3,w4,w5,w6,w7,w8,h,hh,dns1,dns2,bxa,bya,bza,btota,dnst
  double precision bfldp,efldp,e_fluid,e_thermal,v2
  integer*8 i,j,k
  
  efldp=0.
  bfldp=0.
  e_fluid=0.
  e_thermal=0.
  v2=0.

  ! particle energy calculation -- works for 3D only !!!
  dtxi = 1./meshX%dt
  dtyi = 1./meshY%dt
  dtzi = 1./meshZ%dt
  p_xx=0.; p_yy=0.; p_zz=0.

  DO IS=1,1
    DO IIZE = KB-1,KE
      DO IIYE = JB-1,JE
        DO IIXE = 1, NX1
          NP=IPHEAD(IIXE,IIYE,IIZE,IS)

          DO WHILE (NP.NE.0)
            L=NP

          ! Uniform mesh - Same as is in version 5.0
          !  rx=hxi*x(l)+1.5000000000000001
          !  ry=hyi*y(l)+0.5000000000000001d+00
          !  rz=hzi*z(l)+0.5000000000000001d+00
          !  ix=rx
          !  iy=ry
          !  iz=rz
          !  fx=rx-ix
          !  fy=ry-iy
          !  fz=rz-iz

          ! Nonuniform mesh - using MESH_UNMAP
            rx=dtxi*MESH_UNMAP(meshX,x(l))+1.50000000000d+00
            ry=dtyi*MESH_UNMAP(meshY,y(l))+1.50000000000d+00
            rz=dtzi*MESH_UNMAP(meshZ,z(l))+1.50000000000d+00
            ix=rx
            iy=ry
            iz=rz
            fx=rx-ix
            fy=ry-iy
            fz=rz-iz
            iy=iy-1             ! integer index in y direction starts at 0
            iz=iz-1             ! integer index in z direction starts at 0

            ixp1 = ix+1
            iyp1 = iy+1
            izp1 = iz+1

            w1=(1.-fx)*(1.-fy)*(1.-fz)
            w2=    fx *(1.-fy)*(1.-fz)
            w3=(1.-fx)*    fy *(1.-fz)
            w4=    fx*     fy *(1.-fz)
            w5=(1.-fx)*(1.-fy)*    fz
            w6=    fx *(1.-fy)*    fz
            w7=(1.-fx)*    fy*     fz
            w8=    fx*     fy*     fz

            dns1= dns(ix  ,iy  ,iz  ,1)*w1+dns(ixp1,iy  ,iz  ,1)*w2  &
            +     dns(ix  ,iyp1,iz  ,1)*w3+dns(ixp1,iyp1,iz  ,1)*w4  &
            +     dns(ix  ,iy  ,izp1,1)*w5+dns(ixp1,iy  ,izp1,1)*w6  &
            +     dns(ix  ,iyp1,izp1,1)*w7+dns(ixp1,iyp1,izp1,1)*w8

            dns2= 0.

            dnst = dns1 + dns2
  
            vxavg1=vxs(ix  ,iy  ,iz  ,1)*w1+vxs(ixp1,iy  ,iz  ,1)*w2  &
            +      vxs(ix  ,iyp1,iz  ,1)*w3+vxs(ixp1,iyp1,iz  ,1)*w4  &
            +      vxs(ix  ,iy  ,izp1,1)*w5+vxs(ixp1,iy  ,izp1,1)*w6  &
            +      vxs(ix  ,iyp1,izp1,1)*w7+vxs(ixp1,iyp1,izp1,1)*w8

            vxavg2= 0.

            vxavg = (dns1*vxavg1 + dns2*vxavg2)/dnst

            vyavg1=vys(ix  ,iy  ,iz  ,1)*w1+vys(ixp1,iy  ,iz  ,1)*w2  &
            +      vys(ix  ,iyp1,iz  ,1)*w3+vys(ixp1,iyp1,iz  ,1)*w4  &
            +      vys(ix  ,iy  ,izp1,1)*w5+vys(ixp1,iy  ,izp1,1)*w6  &
            +      vys(ix  ,iyp1,izp1,1)*w7+vys(ixp1,iyp1,izp1,1)*w8

            vyavg2=0.

            vyavg = (dns1*vyavg1 + dns2*vyavg2)/dnst

            vzavg1=vzs(ix  ,iy  ,iz  ,1)*w1+vzs(ixp1,iy  ,iz  ,1)*w2  &
            +      vzs(ix  ,iyp1,iz  ,1)*w3+vzs(ixp1,iyp1,iz  ,1)*w4  &
            +      vzs(ix  ,iy  ,izp1,1)*w5+vzs(ixp1,iy  ,izp1,1)*w6  &
            +      vzs(ix  ,iyp1,izp1,1)*w7+vzs(ixp1,iyp1,izp1,1)*w8

            vzavg2=0.

            vzavg = (dns1*vzavg1 + dns2*vzavg2)/dnst

            vxa=vx(l)-vxavg
            vya=vy(l)-vyavg
            vza=vz(l)-vzavg

            xx=vxa*vxa
            yy=vya*vya
            zz=vza*vza

            p_xx (ix  ,iy  ,iz  ,is)=p_xx (ix  ,iy  ,iz  ,is)+qp(np)*w1*xx
            p_xx (ixp1,iy  ,iz  ,is)=p_xx (ixp1,iy  ,iz  ,is)+qp(np)*w2*xx 
            p_xx (ix  ,iyp1,iz  ,is)=p_xx (ix  ,iyp1,iz  ,is)+qp(np)*w3*xx 
            p_xx (ixp1,iyp1,iz  ,is)=p_xx (ixp1,iyp1,iz  ,is)+qp(np)*w4*xx 
            p_xx (ix  ,iy  ,izp1,is)=p_xx (ix  ,iy  ,izp1,is)+qp(np)*w5*xx 
            p_xx (ixp1,iy  ,izp1,is)=p_xx (ixp1,iy  ,izp1,is)+qp(np)*w6*xx 
            p_xx (ix  ,iyp1,izp1,is)=p_xx (ix  ,iyp1,izp1,is)+qp(np)*w7*xx 
            p_xx (ixp1,iyp1,izp1,is)=p_xx (ixp1,iyp1,izp1,is)+qp(np)*w8*xx 

            p_yy (ix  ,iy  ,iz  ,is)=p_yy (ix  ,iy  ,iz  ,is)+qp(np)*w1*yy
            p_yy (ixp1,iy  ,iz  ,is)=p_yy (ixp1,iy  ,iz  ,is)+qp(np)*w2*yy 
            p_yy (ix  ,iyp1,iz  ,is)=p_yy (ix  ,iyp1,iz  ,is)+qp(np)*w3*yy 
            p_yy (ixp1,iyp1,iz  ,is)=p_yy (ixp1,iyp1,iz  ,is)+qp(np)*w4*yy 
            p_yy (ix  ,iy  ,izp1,is)=p_yy (ix  ,iy  ,izp1,is)+qp(np)*w5*yy 
            p_yy (ixp1,iy  ,izp1,is)=p_yy (ixp1,iy  ,izp1,is)+qp(np)*w6*yy 
            p_yy (ix  ,iyp1,izp1,is)=p_yy (ix  ,iyp1,izp1,is)+qp(np)*w7*yy 
            p_yy (ixp1,iyp1,izp1,is)=p_yy (ixp1,iyp1,izp1,is)+qp(np)*w8*yy 

            p_zz (ix  ,iy  ,iz  ,is)=p_zz (ix  ,iy  ,iz  ,is)+qp(np)*w1*zz
            p_zz (ixp1,iy  ,iz  ,is)=p_zz (ixp1,iy  ,iz  ,is)+qp(np)*w2*zz 
            p_zz (ix  ,iyp1,iz  ,is)=p_zz (ix  ,iyp1,iz  ,is)+qp(np)*w3*zz 
            p_zz (ixp1,iyp1,iz  ,is)=p_zz (ixp1,iyp1,iz  ,is)+qp(np)*w4*zz 
            p_zz (ix  ,iy  ,izp1,is)=p_zz (ix  ,iy  ,izp1,is)+qp(np)*w5*zz 
            p_zz (ixp1,iy  ,izp1,is)=p_zz (ixp1,iy  ,izp1,is)+qp(np)*w6*zz 
            p_zz (ix  ,iyp1,izp1,is)=p_zz (ix  ,iyp1,izp1,is)+qp(np)*w7*zz 
            p_zz (ixp1,iyp1,izp1,is)=p_zz (ixp1,iyp1,izp1,is)+qp(np)*w8*zz 

            v2=v2+qp(l)*(vx(l)**2+vy(l)**2+vz(l)**2)

            np=link(np)
          ENDDO
        ENDDO
      ENDDO
    ENDDO


    call XREAL(p_xx (1,jb-1,kb-1,is),NX,NY,NZ)
    call XREAL(p_yy (1,jb-1,kb-1,is),NX,NY,NZ)
    call XREAL(p_zz (1,jb-1,kb-1,is),NX,NY,NZ)

    DO IIZ=KB-1,KE+1
      DO IIY=JB-1,JE+1
        DO IIX=1,NX2
          p_xx(iix,iiy,iiz,is) = p_xx(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_yy(iix,iiy,iiz,is) = p_yy(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
          p_zz(iix,iiy,iiz,is) = p_zz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
        ENDDO
      ENDDO
    ENDDO

  ENDDO

  ! field energy calculation
  do k=kb,ke
      do j=jb,je
        do i=2,nx1
            efldp=efldp+ex(i,j,k)**2+ey(i,j,k)**2+ez(i,j,k)**2
            bfldp=bfldp+bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2
            ! particle energy for species 1 only
            e_fluid=e_fluid+dns(i,j,k,1)*(vxs(i,j,k,1)**2+vys(i,j,k,1)**2+vzs(i,j,k,1)**2)
            e_thermal=e_thermal+p_xx(i,j,k,1)+p_yy(i,j,k,1)+p_yy(i,j,k,1)
        enddo
      enddo
  enddo
  efldp=efldp*hx*hy*hz*0.5 ! assuming uniform grids
  bfldp=bfldp*hx*hy*hz*0.5
  e_fluid=e_fluid*hx*hy*hz*0.5
  e_thermal=e_thermal*hx*hy*hz*0.5
  v2=v2*0.5

  ! collect energies
  call MPI_ALLREDUCE(efldp,efld,1,MPI_DOUBLE_PRECISION,&
        MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(bfldp,bfld,1,MPI_DOUBLE_PRECISION,&
        MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(e_fluid,efluidt,1,MPI_DOUBLE_PRECISION,&
        MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(e_thermal,ethermt,1,MPI_DOUBLE_PRECISION,&
        MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(v2,eptclt,1,MPI_DOUBLE_PRECISION,&
        MPI_SUM,MPI_COMM_WORLD,IERR)
  return
end subroutine energy