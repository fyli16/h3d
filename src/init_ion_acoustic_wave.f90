!=======================================================================
!=======================================================================
      subroutine init_IA_wave
      use parameter_mod
      use MESH2D
      implicit none

      integer*8:: ibp1,ibp2,nptot_max,i,remake,field_subcycle
      double precision:: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi     &
                        ,x_p,y_p,z_p,x_p_logical,y_p_logical        &
                        ,z_p_logical,r_c,q_p,dtsav
 
      integer*8 ip,ipb1,ipb2,is,ixe,iye,ize,j,k
      double precision vxa,vya,vza,vmag,th,ranval(4)

      integer :: seed_size
      integer,allocatable :: seed(:)

      real(kind=8) dB_B0, x_pos,y_pos,z_pos, B0, VA, dB0, dn_n0, dx, f_tot

      real(kind=8) :: kx,ky,kz,kxmin,kymin,kzmin,dvx_,dvy_,dvz_,sin_factor
      real(kind=8) :: loaded_percentage, print_percentage

      double precision:: w1e,w2e,w3e,w4e,w5e,w6e,w7e,w8e
      double precision:: vix1,viy1,viz1
      double precision:: vix2,viy2,viz2
      double precision:: vix3,viy3,viz3
      double precision:: vix4,viy4,viz4
      double precision:: vix5,viy5,viz5
      double precision:: vix6,viy6,viz6
      double precision:: vix7,viy7,viz7
      double precision:: vix8,viy8,viz8
      integer tag,tag0
      integer, parameter :: nbins=65536
      real(kind=8), dimension(nbins) :: f,xx
      real(kind=8), dimension(nbins+1) :: p


      remake = 0

      dtxi = one/meshX%dt
      dtyi = one/meshY%dt
      dtzi = one/meshZ%dt

  
!  initialize the random number generator with a different
!  seed for different processors
 
      call random_seed(size=seed_size) ! find out size of seed
      allocate(seed(seed_size))
      !VR arbirary, but reproducable  numbers for phases
      if (myid==0) then
         call random_seed() ! set default seed
         call random_seed(get=seed) ! get default seed
      endif
      call MPI_BCAST(seed,seed_size,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)  ! propagate master process's version
      call random_seed(put=myid*seed) ! set current seed

      pi=acos(-one)
      it=0
      itfin = 0
      nx1=nx+1
      nx2=nx+2
      ny1=ny+1
      ny2=ny+2
      nz1=nz+1
      nz2=nz+2
      hx=xmax/nx
      hy=ymax/ny
      hz=zmax/nz
      hxi=one/hx
      hyi=one/hy
      hzi=one/hz

!     Nonuniform mesh
        zb=meshZ%xn(kb+1)
        ze=meshZ%xn(ke+2)
        do ipe=0,npes-1
          zbglobal(ipe)=meshZ%xn(kbglobal(ipe)+1)
          zeglobal(ipe)=meshZ%xn(keglobal(ipe)+2)
        enddo
        yb=meshY%xn(jb+1)
        ye=meshY%xn(je+2)
        do ipe=0,npes-1
          ybglobal(ipe)=meshY%xn(jbglobal(ipe)+1)
          yeglobal(ipe)=meshY%xn(jeglobal(ipe)+2)
        enddo
        volume_fraction = (ye-yb)*(ze-zb)/(ymax*zmax)
        !write(*,*)'volume_fraction=',volume_fraction, myid

        xb        = zero
        xe        = xmax
        xb_logical=MESH_UNMAP(meshX,xb)
        xe_logical=MESH_UNMAP(meshX,xe)
        yb_logical=MESH_UNMAP(meshY,yb)
        ye_logical=MESH_UNMAP(meshY,ye)
        zb_logical=MESH_UNMAP(meshZ,zb)
        ze_logical=MESH_UNMAP(meshZ,ze)


        do is=1,nspec
          npm=npx(is)*npy(is)*npz(is)*npes
          dfac(is)=real(ny*nz*nx)/real(npm)
          do ixe=1,nx2 
            do iye=jb-1,je+1
              do ize=kb-1,ke+1
!                 qp_cell(ixe,iye,ize,is) = (meshX%dxc(ixe)*meshY%dxc(iye+1)*meshZ%dxc(ize+1)/(hx*hy*hz))*dfac(is)*frac(is)
                 qp_cell(ixe,iye,ize,is) = meshX%dxc(ixe)*meshY%dxc(iye+1)*meshZ%dxc(ize+1)*dfac(is)*frac(is)
              enddo
            enddo
          enddo
        enddo

      ! allocate the necessary vectors
 
!!!! density fluctuation of an ion acoustic wave
#define DENSITY(k,phi) (1+dn_n0*cos((k)*kxmin*x_pos+ (phi)))
#define DEN DENSITY(3,0)

      !VR: initialize wave parameters

      dn_n0=0.00

      B0 = one/wpiwci
      kxmin = two*pi/xmax
      kymin = two*pi/ymax
      kzmin = two*pi/zmax

      bx=zero;  by=zero; bz = B0;
      ex=zero;ey=zero;ez=zero;
      
    
!  LOAD PARTCILES
     !calculate number of particles in a subdomain
     f_tot=0d0
     dx=xmax/nbins
     do k=1, nbins
       x_pos=(k-1)*dx
       f(k)=DEN
       xx(k)=x_pos
       f_tot=f_tot+f(k)
     enddo
     p(1)=0d0
     do k=1,nbins
       p(k+1)=p(k)+f(k)
     enddo
     p=p/f_tot
     !write(*,*)p(nbins:nbins+1)
     f_tot=f_tot*dx/xmax
     write(*,*)'f_tot=',f_tot

      if (uniform_loading_in_logical_grid) then
        write(*,*)'NOT IMPLEMENTED!'
        stop
      else
        nptotp=0
        do is=1,nspec
          nptotp = nptotp + npx(is)*npy(is)*npz(is)*npes*volume_fraction*f_tot
        enddo
      endif

      call MPI_ALLREDUCE(nptotp,nptot_max,1,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,IERR)
      if (nptot_max > nplmax) then

          remake = 1
          nplmax = 2*nptot_max

          deallocate (x)
          allocate   (x(nplmax))

          deallocate (y)
          allocate   (y(nplmax))

          deallocate (z)
          allocate   (z(nplmax))

          deallocate (vx)
          allocate   (vx(nplmax))

          deallocate (vy)
          allocate   (vy(nplmax))

          deallocate (vz)
          allocate   (vz(nplmax))

          deallocate (qp)
          allocate   (qp(nplmax))

          deallocate (ptag)
          allocate   (ptag(nplmax))

          deallocate (link)
          allocate   (link(nplmax))

          deallocate (porder)
          allocate   (porder(nplmax))

      endif

      if (remake /= 0) call makelist

 
      do is=1,nspec
        ninj(is)=0
        ninj_global(is)=0
        npart(is)=0
        tx0(is)=btspec(is)/(two*wpiwci**2)/wspec(is)
        x0(is)=zero
        x1(is)=xmax
        call MPI_ALLREDUCE(npart(is),npart_global(is),1,MPI_INTEGER8,&
             MPI_SUM,MPI_COMM_WORLD,IERR)
      enddo
      te0=bete/(two*wpiwci**2)
      vbal=one

      if (myid==0) print *, "Initializing particles."
      do is = 1, nspec
        if (myid==0) print *, "species", is
        !tag0=(is-1)*maxtags+maxtags_pe*myid
        tag0=maxtags_pe*nspec*myid+(is-1)*maxtags_pe
        tag=1 

        ipb1 = 1

!       Nonuniform mesh
        if (uniform_loading_in_logical_grid) then
          ipb2 = npx(is)*npy(is)*npz(is)
        else
          ipb2 = npx(is)*npy(is)*npz(is)*npes*volume_fraction*f_tot
        endif

        npm=npx(is)*npy(is)*npz(is)*npes
        dfac(is)=real(ny*nz)*(x1(is)-x0(is))/(hx*real(npm))
        vpar(is)=sqrt(btspec(is)/(wspec(is)*wpiwci*wpiwci))
        vper(is)=vpar(is)*sqrt(anisot(is))

        print_percentage = zero

        !write(*,*)'is, q_p=',is,hx*hy*hz*dfac(is)*frac(is)
        !write(*,*)'myid, is, ipb1, ipb2',myid,is,ipb1,ipb2
        !write(*,*)'myid, is, ipstore',myid,is,ipstore
        do ip=ipb1,ipb2
            
          if (ipb2 == ipb1) write(6,*) "myid = , # particles = ",myid,ipb1,ipb2

          call random_number(harvest=ranval)
          if (uniform_loading_in_logical_grid) then
             X_P_LOGICAL  = XB_LOGICAL+(XE_LOGICAL-XB_LOGICAL)*ranval(1)
             Y_P_LOGICAL  = YB_LOGICAL+(YE_LOGICAL-YB_LOGICAL)*ranval(2)
             Z_P_LOGICAL  = ZB_LOGICAL+(ZE_LOGICAL-ZB_LOGICAL)*ranval(3)
             X_P          = MESH_MAP(meshX,X_P_LOGICAL)
             Y_P          = MESH_MAP(meshY,Y_P_LOGICAL)
             Z_P          = MESH_MAP(meshZ,Z_P_LOGICAL)
             IXE          = dtxi*X_P_LOGICAL+1.50000000000d+00
             IYE          = dtyi*Y_P_LOGICAL+1.50000000000d+00
             IZE          = dtzi*Z_P_LOGICAL+1.50000000000d+00
!            q_p          = (meshX%dxc(ixe)*meshY%dxc(iye)*meshZ%dxc(ize)/(hx*hy*hz))*dfac(is)*frac(is)
             q_p          = meshX%dxc(ixe)*meshY%dxc(iye)*meshZ%dxc(ize)*dfac(is)*frac(is)
          else
            !X_P  = X0(IS)+(X1(IS)-X0(IS))*ranval(1)
            j=1
            do
              if (ranval(1)>p(j+1)) then
                j=j+1
              else
                X_P=xx(j)
                exit
              endif
            enddo
            Y_P  = YB+(YE-YB)*ranval(2)
            Z_P  = ZB+(ZE-ZB)*ranval(3)
              
            q_p  = hx*hy*hz*dfac(is)*frac(is)
          endif
         
          if (q_p == 0) write(6,*) "q_p = ",q_p,ixe,iye,ize
  
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
          tag=tag+1


!         Nonuniform mesh - using MESH_UNMAP
          rxe=dtxi*MESH_UNMAP(meshX,x(np))+1.50000000000d+00
          rye=dtyi*MESH_UNMAP(meshY,y(np))+1.50000000000d+00
          rze=dtzi*MESH_UNMAP(meshZ,z(np))+1.50000000000d+00
          ixe=rxe
          iye=rye
          ize=rze
          iye=iye-1             ! integer index in y direction starts at 0
          ize=ize-1             ! integer index in z direction starts at 0
!
          call random_number(harvest=ranval)
          vmag=sqrt(-log(one-ranval(1)))
          th=two*pi*ranval(2)
          vza=vpar(is)*vmag*cos(th)
          
          vmag=sqrt(-log(one-ranval(3)))
          th=two*pi*ranval(4)
          vxa=vper(is)*vmag*cos(th)
          vya=vper(is)*vmag*sin(th)

          vx(np)=vxa
          vy(np)=vya
          vz(np)=vza

          ipstore=link(np)
          link(np)=iphead(ixe,iye,ize,is)
          iphead(ixe,iye,ize,is)=np
 
 10       CONTINUE
 
          loaded_percentage =  100.0*real(ip-ipb1)/(ipb2-ipb1)
          
          if ((myid==0).and.(loaded_percentage>=print_percentage)) then
             print '(A,F5.1,A)', "loaded ", loaded_percentage," % of particles"
             print_percentage = print_percentage + 5.0d0
          endif
          
       enddo
    enddo
     
      
     if (ndim /= 1) then
        call xrealbcc(ex,1_8,nx,ny,nz)
        call xrealbcc(ey,1_8,nx,ny,nz)
        call xrealbcc(ez,1_8,nx,ny,nz)
     else
        call xrealbcc_pack_e_2d(ex,ey,ez,1_8,nx,ny,nz)
     endif
     
 
 
!***********************
!   set the friction force and resitivity
!   friction force can be zero at t=0
!**********************
 
      do k=kb-1,ke+1
        do j=jb-1,je+1
          do i=1,nx2
            fox(i,j,k)=zero
            foy(i,j,k)=zero
            foz(i,j,k)=zero
            eta(i,j,k)=resis
          enddo
        enddo
      enddo
      dtsav=dt
      dt=zero
      call trans
 
      dt=dtsav
      if (myid == 0) write(6,*) " BEFORE FIELD (CALLED BY INIT)"
      if (.not.testorbt) then
          do field_subcycle=1,n_subcycles
          if (ndim /= 1) then
            call field
          else
            call field_2d
          endif
          enddo
      endif

      
      if (myid == 0) write(6,*) " AFTER  FIELD (CALLED BY INIT)"
      if (ndim /= 1) then
         call etacalc      ! Dietmar's resistivity
      else
         call etacalc_2d   ! Dietmar's resistivity
      endif
      
999   CONTINUE
      
      deallocate(seed) 
      return
    end subroutine
!
!***********************************************************************
