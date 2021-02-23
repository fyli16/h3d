!********************************************************
!
!********************************************************
subroutine accumulate_time_difference(time_begin,time_end &
     &                                     ,time_elapsed)
  
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
      use MESH2D
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

!               Uniform mesh - Same as in version 5.0
!                ixe = hxi*x(np)+1.5000000000000001d+00
!                iye = hyi*y(np)+0.5000000000000001d+00
!                ize = hzi*z(np)+0.5000000000000001d+00

!                 Nonuniform mesh - using MESH_UNMAP
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
!********************************************************
!>    computes field energy ex^2+ey^2+ez^2 and
!!    bx^2+by^2+bz^2
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
!
              DO WHILE (NP.NE.0)
                L=NP

!               Uniform mesh - Same as is in version 5.0
!                rx=hxi*x(l)+1.5000000000000001
!                ry=hyi*y(l)+0.5000000000000001d+00
!                rz=hzi*z(l)+0.5000000000000001d+00
!                ix=rx
!                iy=ry
!                iz=rz
!                fx=rx-ix
!                fy=ry-iy
!                fz=rz-iz

!               Nonuniform mesh - using MESH_UNMAP
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
!********************************************************
!***********************************************************************
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
!
!***********************************************************************
!>    domain decompostion util: splits n elements between numprocs processors
!!    @param n number of elements
!!    @param numprocs number of processors
!!    @param s,e : start/end indices
!***********************************************************************
      subroutine MPE_DECOMP1D( n, numprocs, myid, s, e )
      implicit none  
      integer   numprocs,myid
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
!        open(unit=1,file='.cleanup_status',status='old')
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

subroutine DEBUG(myid,ln)
  implicit none
  integer :: myid,ln
  write(6,*)"myid=",myid,"line number=",ln
end subroutine
