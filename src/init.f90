module m_init
  use m_parameter
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
    call init_arrays
    ! initialize mesh 
    call init_mesh

    ! restart or a fresh start
    call makelist
    if (restart) then 
      call init_restart 
    else
      call init_wavepart 
    endif 

    return 
  end subroutine init_sim


  !---------------------------------------------------------------------
  ! init waves and particles
  !---------------------------------------------------------------------
  subroutine init_wavepart
    use m_eta
    use m_field
    use m_particle

    integer*8 :: ibp1, ibp2, i, remake, field_subcycle
    real*8 :: rxe, rye, rze, fxe, fye, fze, dtxi, dtyi, dtzi, &
              x_p, y_p, z_p, x_p_logical, y_p_logical, z_p_logical, &
              r_c, q_p, dt_save
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
    if (myid==0) then
      print*, "dtxi, dtyi, dtzi = ", dtxi, dtyi, dtzi
    endif 

    ! initialize random number generator for each rank
    call random_seed(size=seed_size) ! obtain seed size
    allocate(seed(seed_size))
    if (myid==0) then
      ! call random_seed() ! method 1: with seed generated by system
      ! call random_seed(get=seed)
      seed = 31415  ! method 2: with a given seed
    endif
    call MPI_BCAST(seed, seed_size, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)  
    call random_seed(put=myid*seed) ! set current seed

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

    ! fraction of local mesh size to global mesh
    volume_fraction = (ye-yb)*(ze-zb)/(ymax*zmax)

    ! what are they for?
    xb_logical = mesh_unmap(meshX,xb)
    xe_logical = mesh_unmap(meshX,xe)
    yb_logical = mesh_unmap(meshY,yb)
    ye_logical = mesh_unmap(meshY,ye)
    zb_logical = mesh_unmap(meshZ,zb)
    ze_logical = mesh_unmap(meshZ,ze)
    ! notice here 'xb' is different from 'meshX%xb'
    ! the former refers to the beginning of x
    ! the latter refers to 'xbb' which is actually the end of x
    if (myid==0) then
      print*, "xb, meshX%xa, meshX%xb, meshX%ta = ", xb, meshX%xa, meshX%xb, meshX%ta
      print*, "xb_logical, xe_logical = ", xb_logical, xe_logical
      print*, "yb_logical, ye_logical = ", yb_logical, ye_logical
      print*, "zb_logical, ze_logical = ", zb_logical, ze_logical
    endif 

    ! what is qp_cell
    do is = 1, nspec
      npm = npx(is)*npy(is)*npz(is)*nprocs
      dfac(is)=real(ny*nz*nx)/real(npm)
      do ixe = 1, nx2 
        do iye = jb-1, je+1
          do ize = kb-1, ke+1
            qp_cell(ixe,iye,ize,is) = meshX%dxc(ixe) * &
                meshY%dxc(iye+1) * meshZ%dxc(ize+1) * dfac(is)*frac(is)
          enddo
        enddo
      enddo
    enddo

    !---------------------------------------------------------------------
    ! initialize wave perturbation on the mesh 
    !---------------------------------------------------------------------
    if (myid==0) then
      print*, " "
      print*, "Initializing wave on the mesh"
      print*, "-------------------------------------------------"
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
          x_pos = meshX%xc(i) ! x has a different indexing              

          ! single Alfven wave
          if (ieta==6) then
            if (k>eta_zs .and. k<nz-eta_zs) then
              bx_ =  dB_B0*B0*sin(kz*z_pos)
              by_ = -dB_B0*B0*cos(kz*z_pos)
            else
              bx_ = 0.
              by_ = 0.
            endif 
          else 
            bx_ =  dB_B0*B0*sin(kz*z_pos)
            by_ = -dB_B0*B0*cos(kz*z_pos)
          endif 
          bz_ = B0
          ex_ = zero  ! why e-fld is zero?
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

          ! use vix to temporarily store values of V on the grid
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
      print*, "Initializing particles"
      print*, "-------------------------------------------------" 
    endif

    do is = 1, nspec
      ninj(is) = 0
      ninj_global(is) = 0
      npart(is) = 0
      tx0(is) = beta_spec(is)/(two*wpiwci**2)/wspec(is) ! plasma temp.
      x0(is) = zero
      x1(is) = xmax
      call MPI_ALLREDUCE(npart(is),npart_global(is),1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
    enddo
    te0 = beta_e/(two*wpiwci**2) ! electron temp.
    if (myid==0) then
      print*, "naprt(is) = ", npart
      print*, "npart_global(is) = ", npart_global
    endif 

    do is = 1, nspec
      tag0 = maxtags_pe*nspec*myid + (is-1)*maxtags_pe
      tag = 1 
      ipb1 = 1

      ! Nonuniform mesh
      if (uniform_load_logical) then
        ipb2 = npx(is)*npy(is)*npz(is)
      else
        ipb2 = npx(is)*npy(is)*npz(is)*nprocs*volume_fraction ! local particle #
      endif
      if (myid==0) print*, "local particle number: ipb2 = ", ipb2

      npm = npx(is)*npy(is)*npz(is)*nprocs
      dfac(is) = real(ny*nz)*(x1(is)-x0(is))/(hx*real(npm))
      vpar(is) = sqrt(beta_spec(is)/(wspec(is)*wpiwci*wpiwci))
      vper(is) = vpar(is)*sqrt(anisot(is))

      if (myid==0) then
        print*
        print*, "species #", is
        print*, "frac     = ", frac(is)
        print*, "npx      = ", npx(is)
        print*, "npy      = ", npy(is)
        print*, "npz      = ", npz(is)
        print*, "dfrac    = ", dfac(is)
        print*, 'npm      = ', npm 
        print*, 'vpar(is) = ', vpar
        print*, 'vper(is) = ', vper
        print*, " "
      endif

      ! actual particle loading
      print_percentage = zero
      do ip = ipb1, ipb2
        call random_number(harvest=ranval)
        if (uniform_load_logical) then
          x_p_logical  = xb_logical+(xe_logical-xb_logical)*ranval(1)
          y_p_logical  = yb_logical+(ye_logical-yb_logical)*ranval(2)
          z_p_logical  = zb_logical+(ze_logical-zb_logical)*ranval(3)
          x_p          = mesh_map(meshX,x_p_logical)
          y_p          = mesh_map(meshY,y_p_logical)
          z_p          = mesh_map(meshZ,z_p_logical)
          ixe          = dtxi*x_p_logical+1.50000000000d+00
          iye          = dtyi*y_p_logical+1.50000000000d+00
          ize          = dtzi*z_p_logical+1.50000000000d+00
          q_p          = meshX%dxc(ixe) * meshY%dxc(iye) * meshZ%dxc(ize) * dfac(is)*frac(is)
        else
          x_p  = x0(is)+(x1(is)-x0(is))*ranval(1)
          y_p  = yb+(ye-yb)*ranval(2)
          z_p  = zb+(ze-zb)*ranval(3)
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
            write(6,"(A,F5.1,A)") "loaded ", loaded_percentage," % of particles"
            print_percentage = print_percentage + 20.0d0
        endif

      enddo
    enddo

    ! updates ghost cells for e-fields
    if (ndim /= 1) then
      call xrealbcc(ex,1_8,nx,ny,nz)
      call xrealbcc(ey,1_8,nx,ny,nz)
      call xrealbcc(ez,1_8,nx,ny,nz)
    else
      call xrealbcc_pack_e_2d(ex,ey,ez,1_8,nx,ny,nz)
    endif
    
    ! initialize friction force and resitivity (temporal)
    ! friction force is zero, and eta is just 'resis'
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

    ! what's done here (no actual particle push)
    dt_save = dt
    dt    = zero ! temporarily set dt=0
    call trans  ! because dt=0, no actual push is done (see 'parmov')
    dt    = dt_save

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

    integer*8 :: ixe, iye, ize, i, j, k, is

    ! read in restart data set and corresponding step of iteration
    if (myid == 0) then
      open(unit=222, file=trim(restart_directory)//'restart_index.dat', status='old')
      read(222,*) restart_index, itrestart 
      print*, " "
      print*, "Restart from set # ", restart_index
      print*, "Restart from iteration # ", itrestart
      print*, " "
      close(222)
    endif
    call MPI_BCAST(restart_index,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
    call MPI_BCAST(itrestart    ,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)

    ! read in restart data
    call write_read_restart_files(-1.0)  
    call MPI_BCAST(itrestart,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
    ! modify start and finish steps of iteration
    itstart = itrestart; it = itstart
    itfinish = (tmax-t_stopped)/dtwci + itstart
    ! swap index for next restart dump
    if (restart_index == 1) then
      restart_index=2
    else
      restart_index=1
    endif
      
    ! Nonuniform mesh
    xb = 0.; xe = xmax

    zb = meshZ%xn(kb+1); ze = meshZ%xn(ke+2)
    do ipe = 0,nprocs-1
      zbglobal(ipe) = meshZ%xn(kbglobal(ipe)+1)
      zeglobal(ipe) = meshZ%xn(keglobal(ipe)+2)
    enddo

    yb = meshY%xn(jb+1); ye = meshY%xn(je+2)
    do ipe = 0,nprocs-1
      ybglobal(ipe) = meshY%xn(jbglobal(ipe)+1)
      yeglobal(ipe) = meshY%xn(jeglobal(ipe)+2)
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
      dfac(is) = real(ny*nz*nx)/real(npm)
      do ixe=1,nx2
          do iye=jb-1,je+1
            do ize=kb-1,ke+1
                qp_cell(ixe,iye,ize,is) = meshX%dxc(ixe)*meshY%dxc(iye+1)*meshZ%dxc(ize+1)*dfac(is)*frac(is)
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