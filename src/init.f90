module m_init
  use m_parameter
  use m_io
  use m_mesh
  implicit none

  contains 

  !---------------------------------------------------------------------
  ! init waves
  !---------------------------------------------------------------------
  subroutine init_waves 
    real*8 :: B0, VA, mi
    real*8 :: kx, ky, kz, kxmin, kymin, kzmin, x_pos, y_pos, z_pos
    real*8 :: bx_, by_, bz_, ex_, ey_, ez_
    real*8 :: dvx_, dvy_, dvz_
    integer :: i, j, k

    if (myid == 0) then
      print*
      print*
      print*, "Initializing waves on the mesh"
      print*, "-------------------------------------------------"
    endif 

    B0 = one/wpiwci  ! RMS amplitude of background B field  
    mi = 0.
    do j = 1, nspec
      mi = mi + frac(j)*wspec(j)
    enddo
    VA = one/wpiwci/sqrt(mi) ! Alfven speed     

    kxmin = two*pi/xmax
    kymin = two*pi/ymax
    kzmin = two*pi/zmax
    kx = zero
    ky = zero
    kz = n_wave_cycles*kzmin

    if (myid == 0) then
      write(6,*) 'loading a single Alfven wave'
      write(6,'(a20,ES15.4)') ' dB_B0 = ', dB_B0
      write(6,'(a20,ES15.4)') ' n_wave_cycles = ', n_wave_cycles
      write(6,'(a20,ES15.4)') ' VA = ', VA
      write(6,*)
      
    endif 

    bx = zero; by = zero; bz = zero
    ex = zero; ey = zero; ez = zero

    do k = kb-1, ke+1
      z_pos = meshZ%xc(k+1)
      do j = jb-1, je+1  
        y_pos = meshY%xc(j+1)
        do i = 1, nx2
          x_pos = meshX%xc(i) ! x has a different indexing              

          ! single Alfven wave
          if (ieta == 6) then
            if (k >= eta_zs .and. k <= nz-eta_zs) then
              bx_ =  dB_B0*B0*sin(kz*z_pos)
              by_ = - dB_B0*B0*cos(kz*z_pos)
            else
              bx_ = 0.
              by_ = 0.
            endif 
          else if (mask .eqv. .true.) then 
            if (k <= nz-mask_zs) then
              bx_ =  dB_B0*B0*sin(kz*z_pos)
              by_ = - dB_B0*B0*cos(kz*z_pos)
            else
              bx_ = 0.
              by_ = 0.
            endif 
          else 
            bx_ =  dB_B0*B0*sin(kz*z_pos)
            by_ = - dB_B0*B0*cos(kz*z_pos)
          endif 

          bz_ = B0
          ex_ = zero  ! e-fields will be determined in field solver
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
  end subroutine init_waves


  !---------------------------------------------------------------------
  ! init particles
  !---------------------------------------------------------------------
  subroutine init_particles
    use m_particle

    integer*8 :: ip, ipb1, ipb2, is, ixe, iye, ize, i, j, k, ibp1, ibp2
    real*8 :: rxe, rye, rze, fxe, fye, fze, r_c, q_p, dt_tmp, &
              x_p, y_p, z_p, x_p_logical, y_p_logical, z_p_logical
    real*8 :: vxa, vya, vza, vmag, th, ranval(4)
    real*8 :: dvx_, dvy_, dvz_
    real*8 :: w1e, w2e, w3e, w4e, w5e, w6e, w7e, w8e
    real*8 :: vix1, viy1, viz1, vix2, viy2, viz2, vix3, viy3, viz3, vix4, viy4, viz4, &
              vix5, viy5, viz5, vix6, viy6, viz6, vix7, viy7, viz7, vix8, viy8, viz8
    integer :: ixep1, iyep1, izep1, tag, tag0
    real*8 :: load_percent, print_percent

    if (myid==0) then
      print*
      print*
      print*, "Initializing particles"
      print*, "-------------------------------------------------" 
    endif

    ! obtain ion/electron temperatures
    do is = 1, nspec
      ninj(is) = 0
      ninj_global(is) = 0
      npart(is) = 0
      tx0(is) = beta_spec(is)/(two*wpiwci**2)/wspec(is) ! ion temp.
      x0(is) = zero; x1(is) = xmax
      call MPI_ALLREDUCE(npart(is),npart_global(is),1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
    enddo
    te0 = beta_elec/(two*wpiwci**2) ! electron temp.
    if (myid==0) then
      write(6,*)
      write(6,'(a20,ES14.6)') 'te0 = ', te0
      do is = 1, nspec
        write(6,'(a15,i1,a,ES14.6)') 'tx0(', is, ') = ', tx0(is)
        write(6,'(a15,i1,a,i6)') 'npart_global(', is, ') = ', npart_global(is)
      enddo 
      print*
    endif 

    ! init particles per species
    call init_seed
    do is = 1, nspec
      tag0 = maxtags_pe*nspec*myid + (is-1)*maxtags_pe
      tag = 1 

      ! determine local start/end indexes of particles
      ipb1 = 1
      ipb2 = nplx(is)*nply(is)*nplz(is)
      print*, 'myid = ', myid, '   ipb2 = ', ipb2
      ! if (uniform_load_logical) then
      !   ipb2 = npx(is)*npy(is)*npz(is)
      ! else
      !   ipb2 = npx(is)*npy(is)*npz(is)*nprocs*volume_fraction
      ! endif

      ! parallel/perpendicular thermal speeds
      ! npm = npx(is)*npy(is)*npz(is)*nprocs ! total particle #
      npm = ppcx(is)*ppcy(is)*ppcz(is)*nx*ny*nz ! total particle #
      dfac(is) = real(ny*nz)*(x1(is)-x0(is))/(hx*real(npm))
      vpar(is) = sqrt( beta_spec(is)/(wspec(is)*wpiwci*wpiwci) )
      vper(is) = vpar(is)*sqrt(anisot(is))

      ! actual loading for this species
      print_percent = zero
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
          x_p  = x0(is) + (x1(is)-x0(is))*ranval(1)
          y_p  = yb + (ye-yb)*ranval(2)
          z_p  = zb + (ze-zb)*ranval(3)
          q_p  = hx*hy*hz*dfac(is)*frac(is)
        endif

        np = ipstore ! np points to the 1st particle
        x(np) = x_p; y(np) = y_p; z(np) = z_p; qp(np) = q_p

        ! tag particles for tracking
        if (tag <= maxtags_pe .and. tag0+tag <= maxtags) then ! only track first maxtags particles
          ptag(np) = tag0 + tag
        else
          ptag(np) = 0 ! do not track
        endif
        tag = tag+1

        ! dtxi=nx; mehs_unmap=(x-xaa)/(xbb-xaa)
        rxe = dtxi*mesh_unmap(meshX,x(np)) + 1.50000000000d+00
        rye = dtyi*mesh_unmap(meshY,y(np)) + 1.50000000000d+00
        rze = dtzi*mesh_unmap(meshZ,z(np)) + 1.50000000000d+00
        ixe = rxe ! take integer part
        iye = rye
        ize = rze
        iye = iye - 1  ! index in y direction starts at 0
        ize = ize - 1  ! index in z direction starts at 0

        fxe = rxe - ixe ! get float part
        fye = rye - iye
        fze = rze - ize

        ixep1 = ixe + 1
        iyep1 = iye + 1
        izep1 = ize + 1

        ! some coefficients in manipulating dvx_
        w1e = (1.-fxe)*(1.-fye)*(1.-fze)
        w2e = fxe*(1.-fye)*(1.-fze)
        w3e = (1.-fxe)*fye*(1.-fze)
        w4e = fxe*fye*(1.-fze)
        w5e = (1.-fxe)*(1.-fye)*fze
        w6e = fxe*(1.-fye)*fze
        w7e = (1.-fxe)*fye*fze
        w8e = fxe*fye*fze
        
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
        
        dvx_ = w1e*vix1+w2e*vix2+w3e*vix3+w4e*vix4+w5e*vix5+w6e*vix6+w7e*vix7+w8e*vix8  
        dvy_ = w1e*viy1+w2e*viy2+w3e*viy3+w4e*viy4+w5e*viy5+w6e*viy6+w7e*viy7+w8e*viy8  
        dvz_ = w1e*viz1+w2e*viz2+w3e*viz3+w4e*viz4+w5e*viz5+w6e*viz6+w7e*viz7+w8e*viz8
          
        ! init thermal speed
        call random_number(harvest=ranval)
        vmag = sqrt(-log(one-ranval(1)))
        th = two*pi*ranval(2)
        vza = vpar(is)*vmag*cos(th) + dvz_
        
        vmag = sqrt(-log(one-ranval(3)))
        th = two*pi*ranval(4)
        vxa = vper(is)*vmag*cos(th) + dvx_
        vya = vper(is)*vmag*sin(th) + dvy_

        ! final ion velocity
        vx(np) = vxa; vy(np) = vya; vz(np) = vza

        ! point to next particle
        ipstore = link(np) ! link(1)-->2-->ipstore
        link(np) = iphead(ixe,iye,ize,is) ! link(1)-->0
        iphead(ixe,iye,ize,is) = np ! iphead=1

        ! print progress of loading
        load_percent = 100.0*real(ip-ipb1)/(ipb2-ipb1)  
        if (myid==0 .and. load_percent>=print_percent) then
            write(6,"(A,F5.1,A)") " loaded ", load_percent," % of particles"
            print_percent = print_percent + 20.0d0
        endif

      enddo
    enddo
    deallocate(seed)

    ! updates ghost cells for e-fields
    if (ndim /= 1) then
      call xrealbcc(ex,1_8,nx,ny,nz)
      call xrealbcc(ey,1_8,nx,ny,nz)
      call xrealbcc(ez,1_8,nx,ny,nz)
    else
      call xrealbcc_pack_e_2d(ex,ey,ez,1_8,nx,ny,nz)
    endif
    
    ! init friction force and resitivity arrays
    do k = kb-1, ke+1
      do j = jb-1, je+1
        do i = 1, nx2
          fox(i,j,k) = zero
          foy(i,j,k) = zero
          foz(i,j,k) = zero
          eta(i,j,k) = resis
        enddo
      enddo
    enddo

    ! 'move' particles with dt=0 (no actual push), just to collect 
    ! moments which will be used in initializing fields at t=0
    dt_tmp = dt
    dt = zero
    call update_particles 
    dt = dt_tmp

  end subroutine init_particles


  !---------------------------------------------------------------------
  ! initialize random number generator for each rank
  !---------------------------------------------------------------------
  subroutine init_seed
    call random_seed(size=seed_size) ! obtain seed size
    allocate(seed(seed_size))
    if (myid==0) then
      ! call random_seed() ! method 1: with seed generated by system
      ! call random_seed(get=seed)
      seed = 31415         ! method 2: with a given seed
    endif
    call MPI_BCAST(seed, seed_size, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)  
    call random_seed(put=myid*seed) ! set current seed
  end subroutine init_seed

end module m_init 