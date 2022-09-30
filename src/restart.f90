module m_restart
  use m_parameter
  use m_mesh
  implicit none 

  contains 

  !---------------------------------------------------------------------
  ! init restart from a previous run
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
    call read_restart 
    call MPI_BCAST(itrestart,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
    ! modify start and finish steps of iteration
    itstart = itrestart; it = itstart
    itfinish = (tmax-t_stopped)/dtwci + itstart
    ! swap index for next restart dump
    if (restart_index == 1) then
      restart_index = 2
    else
      restart_index = 1
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
      ! npm = npx(is)*npy(is)*npz(is)*nprocs
      npm = ppcx(is)*ppcy(is)*ppcz(is)*nx*ny*nz ! total particle #
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

    return 
  end subroutine init_restart


  !---------------------------------------------------------------------
  ! write or read restart files
  ! rw = +1.0: write 
  ! rw = -1.0: read
  !---------------------------------------------------------------------
  subroutine rw_restart(rw)
    integer*8 :: f_unit, np_count, is, ixe, iye, ize, noresete
    real :: rw
    real*8, dimension(:), allocatable :: particle_tmp_array
    integer, dimension(:), allocatable :: particle_tmp_array2 
  
    if (rw == +1.0) then  ! write restart data 
      t_stopped = t_stopped + (it - itstart + 1) * dtwci
      f_unit = 215 + myid
      open(unit=f_unit, file=trim(restart_directory)//'restfld_'//trim(adjustl(myid_char))  &
              //'.bin'//restart_index_suffix(restart_index), form='unformatted', status='unknown')
    
      ! Old restart format - dump all allocated particle memory
      ! write(f_unit) x,y,z,vx,vy,vz,qp,link,porder

      ! New restart format - dump only non-trivial particles
      write(f_unit) nspec

      do is = 1, nspec
        nptotp=0
        do ize = kb-1, ke
          do iye = jb-1, je
            do ixe = 1, nx1   
                np = IPHEAD(ixe, iye, ize, is)
                do while (np.ne.0)
                  nptotp = nptotp + 1
                  np = link(np)
                enddo
            enddo
          enddo
        enddo
        write(f_unit) nptotp

        allocate (particle_tmp_array(nptotp))
        allocate (particle_tmp_array2(nptotp))

        ! x
        nptotp=0
        do ize=kb-1,ke
          do iye=jb-1,je
            do ixe=1,nx1   
                np=IPHEAD(ixe,iye,ize,is)
                do while (np.ne.0)
                  nptotp=nptotp+1
                  particle_tmp_array(nptotp) = x(np)
                  np=link(np)
                enddo
            enddo
          enddo
        enddo
        write(f_unit) particle_tmp_array

        ! y
        nptotp=0
        do ize=kb-1,ke
          do iye=jb-1,je
            do ixe=1,nx1   
                np=IPHEAD(ixe,iye,ize,is)
                do while (np.ne.0)
                  nptotp=nptotp+1
                  particle_tmp_array(nptotp) = y(np)
                  np=link(np)
                enddo
            enddo
          enddo
        enddo
        write(f_unit) particle_tmp_array

        ! z
        nptotp=0
        do ize=kb-1,ke
          do iye=jb-1,je
            do ixe=1,nx1   
                np=IPHEAD(ixe,iye,ize,is)
                do while (np.ne.0)
                  nptotp=nptotp+1
                  particle_tmp_array(nptotp) = z(np)
                  np=link(np)
                enddo
            enddo
          enddo
        enddo
        write(f_unit) particle_tmp_array

        ! vx
        nptotp=0
        do ize=kb-1,ke
          do iye=jb-1,je
            do ixe=1,nx1   
                np=IPHEAD(ixe,iye,ize,is)
                do while (np.ne.0)
                  nptotp=nptotp+1
                  particle_tmp_array(nptotp) = vx(np)
                  np=link(np)
                enddo
            enddo
          enddo
        enddo
        write(f_unit) particle_tmp_array

        ! vy
        nptotp=0
        do ize=kb-1,ke
          do iye=jb-1,je
            do ixe=1,nx1   
                np=IPHEAD(ixe,iye,ize,is)
                do while (np.ne.0)
                  nptotp=nptotp+1
                  particle_tmp_array(nptotp) = vy(np)
                  np=link(np)
                enddo
            enddo
          enddo
        enddo
        write(f_unit) particle_tmp_array

        ! vz
        nptotp=0
        do ize=kb-1,ke
          do iye=jb-1,je
            do ixe=1,nx1   
                np=IPHEAD(ixe,iye,ize,is)
                do while (np.ne.0)
                  nptotp=nptotp+1
                  particle_tmp_array(nptotp) = vz(np)
                  np=link(np)
                enddo
            enddo
          enddo
        enddo
        write(f_unit) particle_tmp_array

        ! qp
        nptotp=0
        do ize=kb-1,ke
          do iye=jb-1,je
            do ixe=1,nx1   
                np=IPHEAD(ixe,iye,ize,is)
                do while (np.ne.0)
                  nptotp=nptotp+1
                  particle_tmp_array(nptotp) = qp(np)
                  np=link(np)
                enddo
            enddo
          enddo
        enddo
        write(f_unit) particle_tmp_array

        ! ptag
        nptotp=0
        do ize=kb-1,ke
          do iye=jb-1,je
            do ixe=1,nx1   
                np=IPHEAD(ixe,iye,ize,is)
                do while (np.ne.0)
                  nptotp=nptotp+1
                  particle_tmp_array2(nptotp) = ptag(np)
                  np=link(np)
                enddo
            enddo
          enddo
        enddo
        write(f_unit) particle_tmp_array2

        deallocate (particle_tmp_array)
        deallocate (particle_tmp_array2)

      enddo

      write(f_unit) ninj,ninj_global,nescape,nescape_global,npart, &
      npart_global,t_stopped

      write(f_unit) x0,x1,tx0,vpar,vper

      write(f_unit) beta_spec, qspec, wspec, frac,                    &
      anisot, denmin, resis, wpiwci, beta_elec,                &
      xmax,ymax,zmax,gamma,                              &
      nplx, nply, nplz,                                               &
      n_sub_b, eta_par, netax, netay, nspec,   &
      nx, ny, nz, etamin, etamax, ieta

      write(f_unit) hx,hy,hz,hxi,hyi,hzi                           &
      ,pi,efld,bfld,efluid,ethermal,eptcl,time,te0                &
      ,itrestart,iwt                                   &
      ,nx1,nx2,ny1,ny2,nz1,nz2,it                                  &
      ! ,ipstore,nptot,npleaving,npentering                &
      ,nptot,npleaving,npentering                        &
      ,iclock_speed,iopen,iseed, file_unit,file_unit_read          &
      ,clock_init,clock_old,clock_now                   &
      ,clock_time1

      write(f_unit) dfac,nskip,ipleft,iprite,ipsendleft,ipsendrite &
      ,iprecv,ipsendtop,ipsendbot,ipsendlefttop,ipsendleftbot      &
      ,ipsendritetop,ipsendritebot,ipsend

      write(f_unit) bx,by,bz,den,pe,eta,ex,ey,ez,fox,foy,foz       &
      ,eta_times_b_dot_j

      write(f_unit) vix, viy, viz, vixo, viyo, vizo

      ! write(f_unit) tpar,tperp,dns,vxs,vys,vzs,iphead,iptemp
      write(f_unit) tpar,tperp,dns,vxs,vys,vzs

      ! save user data
      ! call user_data_write_restart(f_unit)
      
      close(unit=f_unit)

    else if (rw == -1.0) then ! read restart data

      f_unit = 215 + myid
      open(unit=f_unit, file=trim(restart_directory)//'restfld_'//trim(adjustl(myid_char))  &
              //'.bin'//restart_index_suffix(restart_index), form='unformatted', status='unknown')
        
      ! Old restart format - dump all allocated particle memory
      ! read(f_unit) x,y,z,vx,vy,vz,qp,link,porder

      ! New restart format - dump only non-trivial particles
      read(f_unit) nspec

      do is = 1, nspec
        read(f_unit) nptotp
        allocate (particle_tmp_array(nptotp))
        allocate (particle_tmp_array2(nptotp))

        ! x
        read(f_unit) particle_tmp_array
        ixe=2
        iye=jb
        ize=kb
        do np_count = 1, nptotp
          np = ipstore
          x(np)=particle_tmp_array(np_count)
          ipstore=link(np)
          link(np)=iphead(ixe,iye,ize,is)
          iphead(ixe,iye,ize,is)=np
        enddo
        iptemp(ixe,iye,ize,is)=0
        np=iphead(ixe,iye,ize,is)
        do while (np.ne.0)
          iphead(ixe,iye,ize,is)=link(np)
          link(np)=iptemp(ixe,iye,ize,is)
          iptemp(ixe,iye,ize,is)=np
          np=iphead(ixe,iye,ize,is)
        enddo
        iphead(ixe,iye,ize,is)=iptemp(ixe,iye,ize,is)
        iptemp(ixe,iye,ize,is)=0

        ! y
        read(f_unit) particle_tmp_array
        np_count=0
        np=iphead(ixe,iye,ize,is)
        do while (np.ne.0)
          np_count=np_count+1
          y(np)=particle_tmp_array(np_count)
          np=link(np)
        enddo

        ! z
        read(f_unit) particle_tmp_array
        np_count=0
        np=iphead(ixe,iye,ize,is)
        do while (np.ne.0)
          np_count=np_count+1
          z(np)=particle_tmp_array(np_count)
          np=link(np)
        enddo

        ! vx
        read(f_unit) particle_tmp_array
        np_count=0
        np=iphead(ixe,iye,ize,is)
        do while (np.ne.0)
          np_count=np_count+1
          vx(np)=particle_tmp_array(np_count)
          np=link(np)
        enddo

        ! vy
        read(f_unit) particle_tmp_array
        np_count=0
        np=iphead(ixe,iye,ize,is)
        do while (np.ne.0)
          np_count=np_count+1
          vy(np)=particle_tmp_array(np_count)
          np=link(np)
        enddo
    
        ! vz
        read(f_unit) particle_tmp_array
        np_count=0
        np=iphead(ixe,iye,ize,is)
        do while (np.ne.0)
          np_count=np_count+1
          vz(np)=particle_tmp_array(np_count)
          np=link(np)
        enddo

        ! qp
        read(f_unit) particle_tmp_array
        np_count=0
        np=iphead(ixe,iye,ize,is)
        do while (np.ne.0)
          np_count=np_count+1
          qp(np)=particle_tmp_array(np_count)
          np=link(np)
        enddo

        ! ptag
        read(f_unit) particle_tmp_array2
        np_count=0
        np=iphead(ixe,iye,ize,is)
        do while (np.ne.0)
          np_count=np_count+1
          ptag(np)=particle_tmp_array2(np_count)
          np=link(np)
        enddo

        deallocate (particle_tmp_array)
        deallocate (particle_tmp_array2)

      enddo

      read(f_unit) ninj,ninj_global,nescape,nescape_global,npart,  &
      npart_global,t_stopped

      read(f_unit) x0,x1,tx0,vpar,vper

      read(f_unit) beta_spec, qspec, wspec, frac,       &
      anisot, denmin, resis, wpiwci, beta_elec,       &
      xmax,ymax,zmax,gamma,                  &
      nplx, nply, nplz,                                               &
      n_sub_b, eta_par, netax, netay, nspec,   &
      nx, ny, nz, etamin, etamax, ieta

      read(f_unit) hx,hy,hz,hxi,hyi,hzi                            &
      ,efld,bfld,efluid,ethermal,eptcl,time,te0                                        &
      ,itrestart,iwt                     &
      ,nx1,nx2,ny1,ny2,nz1,nz2,it                                   &
      ! ,ipstore,nptot,npleaving,npentering                 &
      ,nptot,npleaving,npentering                 &
      ,iclock_speed,iopen,iseed, file_unit,file_unit_read           &
      ,clock_init,clock_old,clock_now                    &
      ,clock_time1

      read(f_unit) dfac,nskip,ipleft,iprite,ipsendleft,ipsendrite  &
      ,iprecv,ipsendtop,ipsendbot,ipsendlefttop,ipsendleftbot      &
      ,ipsendritetop,ipsendritebot,ipsend

      read(f_unit) bx,by,bz,den,pe,eta,ex,ey,ez,fox,foy,foz        &
      ,eta_times_b_dot_j

      read(f_unit) vix, viy, viz, vixo, viyo, vizo

      ! read(f_unit) tpar,tperp,dns,vxs,vys,vzs,iphead,iptemp
      read(f_unit) tpar,tperp,dns,vxs,vys,vzs

      ! save user data
      ! call user_diagnostics_restart(f_unit)        

      close(unit=f_unit)

      ! Reset electic field and fox,y,z 
      noresete = 1
      if (noresete == 0) then
        deno=den; vixo=vix; viyo=viy; vizo=viz
        call ecalc( 1 )
        call focalc
      endif

    endif

    return
  end subroutine rw_restart


  !---------------------------------------------------------------------
  ! write restart
  !---------------------------------------------------------------------
  subroutine write_restart
    call rw_restart(1.0)
  end subroutine write_restart


  !---------------------------------------------------------------------
  ! read restart
  !---------------------------------------------------------------------
  subroutine read_restart
    call rw_restart(-1.0)
  end subroutine read_restart

end module m_restart
