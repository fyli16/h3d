    !---------------------------------------------------------------------
    ! init waves
    !---------------------------------------------------------------------
    subroutine init_wave
      use parameter_mod
      use mesh2d
      implicit none

      integer*8 :: ibp1,ibp2,nptot_max,i,remake,field_subcycle
      real*8 :: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi     &
                        ,x_p,y_p,z_p,x_p_logical,y_p_logical        &
                        ,z_p_logical,r_c,q_p,dtsav
      integer*8 :: ip,ipb1,ipb2,is,ixe,iye,ize,j,k
      real*8 :: vxa,vya,vza,vmag,th,ranval(4)
      integer :: seed_size
      integer, allocatable :: seed(:)
      real*8 :: x_pos, y_pos, z_pos, B0, VA, dB0, mi
      real*8 :: bx_, by_, bz_, ex_, ey_, ez_
      real*8 :: kx,ky,kz,kxmin,kymin,kzmin,dvx_,dvy_,dvz_,sin_factor
      real*8 :: loaded_percentage, print_percentage
      real*8 :: w1e,w2e,w3e,w4e,w5e,w6e,w7e,w8e
      real*8 :: vix1,viy1,viz1
      real*8 :: vix2,viy2,viz2
      real*8 :: vix3,viy3,viz3
      real*8 :: vix4,viy4,viz4
      real*8 :: vix5,viy5,viz5
      real*8 :: vix6,viy6,viz6
      real*8 :: vix7,viy7,viz7
      real*8 :: vix8,viy8,viz8
      integer :: ixep1,iyep1,izep1
      integer :: tag,tag0

      if (myid==0) then
        write(6,*) " "
        write(6,*) "Initializing wave ..."
      endif 

      remake = 0  !??

      dtxi = one/meshX%dt
      dtyi = one/meshY%dt
      dtzi = one/meshZ%dt

      ! initialize the random number generator with a different
      ! seed for different processors
      call random_seed(size=seed_size) ! find out size of seed
      allocate(seed(seed_size))
      ! VR arbirary, but reproducable numbers for phases
      if (myid==0) then
        call random_seed() ! set default seed
        call random_seed(get=seed) ! get default seed
      endif
      ! propagate master process's version
      call MPI_BCAST(seed,seed_size,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)  
      call random_seed(put=myid*seed) ! set current seed

      it=0; itfin = 0
      nx1 = nx+1; nx2 = nx+2
      ny1 = ny+1; ny2 = ny+2
      nz1 = nz+1; nz2 = nz+2
      hx = xmax/nx; hy = ymax/ny; hz = zmax/nz
      hxi = one/hx; hyi = one/hy; hzi = one/hz

      ! Nonuniform mesh
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

      xb        = zero
      xe        = xmax
      xb_logical = MESH_UNMAP(meshX,xb)
      xe_logical = MESH_UNMAP(meshX,xe)
      yb_logical = MESH_UNMAP(meshY,yb)
      ye_logical = MESH_UNMAP(meshY,ye)
      zb_logical = MESH_UNMAP(meshZ,zb)
      ze_logical = MESH_UNMAP(meshZ,ze)

      do is = 1, nspec
        npm = npx(is)*npy(is)*npz(is)*npes
        dfac(is)=real(ny*nz*nx)/real(npm)
        do ixe=1,nx2 
          do iye=jb-1,je+1
            do ize=kb-1,ke+1
              qp_cell(ixe,iye,ize,is) = meshX%dxc(ixe)*meshY%dxc(iye+1)*meshZ%dxc(ize+1)*dfac(is)*frac(is)
            enddo
          enddo
        enddo
      enddo

! allocate the necessary vectors

! Alfvenic perturbation with deltaB in the x direction
! i,j,k are wave numbers in x,y,z
#define DBX_1(k,j,phi) (dB_B0*B0*cos((k)*kzmin*z_pos + (j)*kymin*y_pos + (phi)))
#define DEY_1(k,j,phi) (-dB_B0*(k/abs(k))*VA*B0*cos((k)*kzmin*z_pos + (j)*kymin*y_pos + (phi)))
! These give velocity & current consistent with Alfven wave
#define DUX_1(k,j,phi) (-dB_B0*(k/abs(k))*VA*cos((k)*kzmin*z_pos + (j)*kymin*y_pos + (phi)))
#define DJY_1(k,j,phi) (-dB_B0*B0*(k)*kzmin*sin((k)*kzmin*z_pos + (j)*kymin*y_pos + (phi)))
#define DJZ_1(k,j,phi) (dB_B0*B0*(j)*kymin*sin((k)*kzmin*z_pos + (j)*kymin*y_pos + (phi)))
#define BX_PERT_1 DBX_1(1,1,0) + DBX_1(1,2,1.5) + DBX_1(-2,3,3.9)   
#define EY_PERT_1 DEY_1(1,1,0) + DEY_1(1,2,1.5) + DEY_1(-2,3,3.9)  
#define UX_PERT_1 DUX_1(1,1,0) + DUX_1(1,2,1.5) + DUX_1(-2,3,3.9)
#define JY_PERT_1 DJY_1(1,1,0) + DJY_1(1,2,1.5) + DJY_1(-2,3,3.9)
#define JZ_PERT_1 DJZ_1(1,1,0) + DJZ_1(1,2,1.5) + DJZ_1(-2,3,3.9)

! Alfvenic perturbation with deltaB in the y direction
! works only for a pair plasma
#define DBY_2(k,i,phi) (dB_B0*B0*cos((k)*kzmin*z_pos + (i)*kxmin*x_pos + (phi)))
#define DEX_2(k,i,phi) (dB_B0*(k/abs(k))*VA*B0*cos((k)*kzmin*z_pos + (i)*kxmin*x_pos + (phi)))
! These give velocity & current consistent with Alfven wave
#define DUY_2(k,i,phi) (-dB_B0*(k/abs(k))*VA*cos((k)*kzmin*z_pos + (i)*kxmin*x_pos + (phi)))
#define DJX_2(k,i,phi) (dB_B0*B0*(k)*kzmin*sin((k)*kzmin*z_pos + (i)*kxmin*x_pos + (phi)))
#define DJZ_2(k,i,phi) (-dB_B0*B0*(i)*kxmin*sin((k)*kzmin*z_pos + (i)*kxmin*x_pos + (phi)))
#define BY_PERT_2 DBY_2(-1,1,0.4) + DBY_2(-1,-2,2.56) + DBY_2(2,-3,4.19)   
#define EX_PERT_2 DEX_2(-1,1,0.4) + DEX_2(-1,-2,2.56) + DEX_2(2,-3,4.19) 
#define UY_PERT_2 DUY_2(-1,1,0.4) + DUY_2(-1,-2,2.56) + DUY_2(2,-3,4.19)
#define JX_PERT_2 DJX_2(-1,1,0.4) + DJX_2(-1,-2,2.56) + DJX_2(2,-3,4.19)
#define JZ_PERT_2 DJZ_2(-1,1,0.4) + DJZ_2(-1,-2,2.56) + DJZ_2(2,-3,4.19)

      !VR: initialize wave parameters
      ! RMS amplitude of the pertubation [B0=RMS(B)]
      ! dB_B0=0.1  ! now get it from input      
      B0 = one/wpiwci
      !VR Alfven speed
      mi=0.
      do j =1, nspec
        mi = mi + frac(j)*wspec(j)
      enddo
      VA = one/wpiwci/sqrt(mi)      

      kxmin = two*pi/xmax
      kymin = two*pi/ymax
      kzmin = two*pi/zmax
      kx = zero
      ky = zero
      kz = num_cycles*kzmin  ! num_cycles is taken from input

      bx = zero; by = zero; bz = zero
      ex = zero; ey = zero; ez = zero
     
      ! initialie perturbation on the grid 
      if (myid==0) then
        write(6,*) "  Initializing waves on the mesh grid ..."
      endif 

      do k=kb-1,ke+1
        z_pos = meshZ%xc(k+1)
        do j=jb-1,je+1  
          y_pos = meshY%xc(j+1)
          do i=1,nx2
            x_pos = meshX%xc(i)   !VR this is not a typo. For some reason, x has different indexing compared to y and z (!!!)               
            ! no wave
            ! bx_ = zero
            ! by_ = zero
            ! bz_ = B0
            ! ex_ = zero
            ! ey_ = zero
            ! ez_ = zero
            ! dvx_ = zero
            ! dvy_ = zero
            ! dvz_ = zero

            ! single Alfven wave
            bx_ =  dB_B0*B0*sin(kz*z_pos)
            by_ = -dB_B0*B0*cos(kz*z_pos)
            bz_ = B0
            ex_ = zero
            ey_ = zero
            ez_ = zero
            dvx_ = -VA*bx_/B0 
            dvy_ = -VA*by_/B0 
            dvz_ = zero

            ! multiple waves
            ! bx_ = BX_PERT_1
            ! by_ = BY_PERT_2
            ! bz_ = B0
            ! ex_ = EX_PERT_2
            ! ey_ = EY_PERT_1
            ! ez_ = zero
            ! dvx_ = UX_PERT_1
            ! dvy_ = UY_PERT_2
            ! dvz_ = zero

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
     
      ! Load particles
      if (uniform_loading_in_logical_grid) then
        nptotp=0
        do is=1,nspec
          nptotp = nptotp + npx(is)*npy(is)*npz(is)
        enddo
      else
        nptotp=0
        do is=1,nspec
          nptotp = nptotp + npx(is)*npy(is)*npz(is)*npes*volume_fraction
        enddo
      endif

      call MPI_ALLREDUCE(nptotp,nptot_max,1,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,IERR)
      if (nptot_max > nplmax) then
        remake = 1
        nplmax = 2*nptot_max
        deallocate(x); allocate(x(nplmax))
        deallocate(y); allocate(y(nplmax))
        deallocate(z); allocate(z(nplmax))
        deallocate(vx); allocate(vx(nplmax))
        deallocate(vy); allocate(vy(nplmax))
        deallocate(vz); allocate(vz(nplmax))
        deallocate(qp); allocate(qp(nplmax))
        deallocate(ptag); allocate(ptag(nplmax))
        deallocate(link); allocate(link(nplmax))
        deallocate(porder); allocate(porder(nplmax))
      endif

      if (remake /= 0) call makelist

      do is = 1, nspec
        ninj(is)=0
        ninj_global(is)=0
        npart(is)=0
        tx0(is)=btspec(is)/(two*wpiwci**2)/wspec(is)
        x0(is)=zero
        x1(is)=xmax
        call MPI_ALLREDUCE(npart(is),npart_global(is),1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
      enddo
      te0=bete/(two*wpiwci**2)
      vbal=one

      if (myid==0) write(6,*) "  Initializing particles ..." 
      do is = 1, nspec
        tag0 = maxtags_pe*nspec*myid + (is-1)*maxtags_pe
        tag=1 
        ipb1 = 1

        ! Nonuniform mesh
        if (uniform_loading_in_logical_grid) then
          ipb2 = npx(is)*npy(is)*npz(is)
        else
          ipb2 = npx(is)*npy(is)*npz(is)*npes*volume_fraction
        endif

        npm = npx(is)*npy(is)*npz(is)*npes
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
          write(6,*) "  nptot_max = ", nptot_max
          write(6,*) "  q_p =", hx*hy*hz*dfac(is)*frac(is)
          write(6,*) " "
        endif

        print_percentage = zero

        do ip=ipb1,ipb2
          if (ipb2 == ipb1) write(6,*) "myid = , # particles = ", myid,ipb1,ipb2
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
            q_p          = meshX%dxc(ixe)*meshY%dxc(iye)*meshZ%dxc(ize)*dfac(is)*frac(is)
          else
            X_P  = X0(IS)+(X1(IS)-X0(IS))*ranval(1)
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
          tag = tag+1

          ! Nonuniform mesh - using MESH_UNMAP
          rxe=dtxi*MESH_UNMAP(meshX,x(np))+1.50000000000d+00
          rye=dtyi*MESH_UNMAP(meshY,y(np))+1.50000000000d+00
          rze=dtzi*MESH_UNMAP(meshZ,z(np))+1.50000000000d+00
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
          
          dvx_=w1e*vix1+w2e*vix2+w3e*vix3+w4e*vix4   &
               +w5e*vix5+w6e*vix6+w7e*vix7+w8e*vix8  
          dvy_=w1e*viy1+w2e*viy2+w3e*viy3+w4e*viy4   &
               +w5e*viy5+w6e*viy6+w7e*viy7+w8e*viy8  
          dvz_=w1e*viz1+w2e*viz2+w3e*viz3+w4e*viz4   &
               +w5e*viz5+w6e*viz6+w7e*viz7+w8e*viz8  
                    
          !interpolate V at the particle position from pre-computed values at the grid

          ! dvx_ = -dB_B0*VA*sin(kz*z_p)
            
          call random_number(harvest=ranval)
          vmag=sqrt(-log(one-ranval(1)))
          th=two*pi*ranval(2)
          vza=vpar(is)*vmag*cos(th) + dvz_
          
          vmag=sqrt(-log(one-ranval(3)))
          th=two*pi*ranval(4)
          vxa=vper(is)*vmag*cos(th) + dvx_
          vya=vper(is)*vmag*sin(th) + dvy_

          vx(np)=vxa
          vy(np)=vya
          vz(np)=vza

          ipstore=link(np)
          link(np)=iphead(ixe,iye,ize,is)
          iphead(ixe,iye,ize,is)=np
 
 10       CONTINUE
 
          loaded_percentage =  100.0*real(ip-ipb1)/(ipb2-ipb1)
          
          if ((myid==0).and.(loaded_percentage>=print_percentage)) then
             write(6,"(A,F5.1,A)") "  loaded ", loaded_percentage," % of particles"
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
     
      !   set the friction force and resitivity
      !   friction force can be zero at t=0
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

      ! why do we need to set dt=0, and push particles once?
      dtsav=dt
      dt=zero
      call trans
      dt=dtsav

      if (.not.testorbt) then
        do field_subcycle = 1, n_subcycles ! n_subcycles=0, so this block is not executed
          if (ndim /= 1) then
            call field
          else
            call field_2d
          endif
        enddo
      endif

      if (ndim /= 1) then
         call eta_calc      ! Dietmar's resistivity
      else
         call eta_calc_2d   ! Dietmar's resistivity
      endif
      
999   CONTINUE
      
      deallocate(seed) 
      return
    end subroutine init_wave