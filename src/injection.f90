module m_injection
  use m_parameter
  ! use m_io
  ! use m_mesh
  implicit none

  contains 

  !---------------------------------------------------------------------
  ! inject 3D RMF antenna waves via B field
  !---------------------------------------------------------------------
  subroutine inject_waves_b_rmf
    real*8 :: B0, VA, mi
    real*8 :: kx, ky, kz, kxmin, kymin, kzmin
    real*8 :: inj_time
    real*8 :: bx_, by_, bz_
    integer :: i, j, k
    integer :: iw  ! index of wave
    real*8 :: time_env, radial_env  ! wave envelope

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
    inj_time = it * dtwci  ! in 1/wci

    do iw = 1, 4 
      bx_ = 0.0; by_ = 0.0

      if ( inj_dB_B0(iw)>0.0 .and. (kb-1)<=inj_z_pos(iw) .and. inj_z_pos(iw)<=ke+1 ) then 

        if (inj_time <=inj_t_upramp(iw)) then
          time_env = (sin(0.5*pi*inj_time/inj_t_upramp(iw)))**2.
        else if (inj_time <= inj_t_upramp(iw)+inj_t_flat(iw)) then
          time_env = 1.0
        else if (inj_time <= inj_t_upramp(iw)+inj_t_flat(iw)+inj_t_downramp(iw)) then
          time_env = (cos(0.5*pi*(inj_time-inj_t_upramp(iw)-inj_t_flat(iw))/inj_t_downramp(iw)))**2.
        else 
          time_env = 0.0
        endif 

        ! bx_ = inj_dB_B0(iw)*B0*time_env*inj_rmf_ampl_corr 
        ! by_ = inj_dB_B0(iw)*B0*time_env*inj_rmf_ampl_corr
        bx_ = inj_dB_B0(iw)*B0*time_env
        by_ = inj_dB_B0(iw)*B0*time_env

        kz = inj_wave_cycles(iw) * kzmin
        do j = jb, je
          do i = 1, nx
            if (inj_wave_pol(iw)==0) then  ! x-pol
              bx(i,j,inj_z_pos(iw)) = bx_*rmf_bx1(i,j)*sin(kz*inj_time)
              by(i,j,inj_z_pos(iw)) = by_*rmf_by1(i,j)*sin(kz*inj_time)
            else if (inj_wave_pol(iw)==1) then  ! left-hand
              bx(i,j,inj_z_pos(iw)) = bx_*rmf_bx1(i,j)*sin(kz*inj_time) + bx_*rmf_bx2(i,j)*cos(kz*inj_time)
              by(i,j,inj_z_pos(iw)) = by_*rmf_by1(i,j)*sin(kz*inj_time) + by_*rmf_by2(i,j)*cos(kz*inj_time)
            else if (inj_wave_pol(iw)==-1) then  ! right-hand
              bx(i,j,inj_z_pos(iw)) = bx_*rmf_bx1(i,j)*cos(kz*inj_time) + bx_*rmf_bx2(i,j)*sin(kz*inj_time)
              by(i,j,inj_z_pos(iw)) = by_*rmf_by1(i,j)*cos(kz*inj_time) + by_*rmf_by2(i,j)*sin(kz*inj_time)
            endif 
          enddo
        enddo
        ! do j = jb-1, je+1
        !   do i = 1, nx2
        !     if (inj_wave_pol(iw)==0) then  ! x-pol
        !       bx(i,j,inj_z_pos(iw)) = bx_*rmf_bx1(i,j)*sin(kz*inj_time)
        !       by(i,j,inj_z_pos(iw)) = by_*rmf_by1(i,j)*sin(kz*inj_time)
        !     else if (inj_wave_pol(iw)==1) then  ! left-hand
        !       bx(i,j,inj_z_pos(iw)) = bx_*rmf_bx1(i,j)*sin(kz*inj_time) + bx_*rmf_bx2(i,j)*cos(kz*inj_time)
        !       by(i,j,inj_z_pos(iw)) = by_*rmf_by1(i,j)*sin(kz*inj_time) + by_*rmf_by2(i,j)*cos(kz*inj_time)
        !     else if (inj_wave_pol(iw)==-1) then  ! right-hand
        !       bx(i,j,inj_z_pos(iw)) = bx_*rmf_bx1(i,j)*cos(kz*inj_time) + bx_*rmf_bx2(i,j)*sin(kz*inj_time)
        !       by(i,j,inj_z_pos(iw)) = by_*rmf_by1(i,j)*cos(kz*inj_time) + by_*rmf_by2(i,j)*sin(kz*inj_time)
        !     endif 
        !   enddo
        ! enddo

      endif ! end inj_dB_B0(iw)
    enddo ! end wave indexing (iw)

  end subroutine inject_waves_b_rmf


  !---------------------------------------------------------------------
  ! inject waves via B field
  !---------------------------------------------------------------------
  subroutine inject_waves_b
    real*8 :: B0, VA, mi
    real*8 :: kx, ky, kz, kxmin, kymin, kzmin
    real*8 :: inj_time
    real*8 :: bx_, by_, bz_
    integer :: i, j, k
    integer :: iw  ! index of wave
    real*8 :: time_env, radial_env  ! wave envelope

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
    ! kz = wave_cycles*kzmin
    inj_time = it * dtwci  ! in 1/wci

    do iw = 1, 4 
      bx_ = 0.0; by_ = 0.0

      if ( inj_dB_B0(iw)>0.0 .and. (kb-1)<=inj_z_pos(iw) .and. inj_z_pos(iw)<=ke+1 ) then 

        if (inj_time <=inj_t_upramp(iw)) then
          time_env = (sin(0.5*pi*inj_time/inj_t_upramp(iw)))**2.
        else if (inj_time <= inj_t_upramp(iw)+inj_t_flat(iw)) then
          time_env = 1.0
        else if (inj_time <= inj_t_upramp(iw)+inj_t_flat(iw)+inj_t_downramp(iw)) then
          time_env = (cos(0.5*pi*(inj_time-inj_t_upramp(iw)-inj_t_flat(iw))/inj_t_downramp(iw)))**2.
        else 
          time_env = 0.0
        endif 

        kz = inj_wave_cycles(iw) * kzmin
        if (inj_wave_pol(iw)==0) then  ! x-pol
          bx_ = inj_dB_B0(iw)*B0*time_env*sin(kz*inj_time)
        else if (inj_wave_pol(iw)==1) then ! y-pol, left-hand pol
          by_ = inj_dB_B0(iw)*B0*time_env*cos(kz*inj_time)
        else if (inj_wave_pol(iw)==-1) then ! y-pol, right-hand pol
          by_ = -inj_dB_B0(iw)*B0*time_env*cos(kz*inj_time)
        endif

        if (inj_wave_radius(iw)==0) then ! inject at all x, y
          do j = jb-1, je+1
            do i = 1, nx2
              ! add injection value to previous wave if they have the same injection position
              if ( iw>1 .and. inj_z_pos(iw)==inj_z_pos(iw-1) ) then
                bx(i,j,inj_z_pos(iw)) = bx(i,j,inj_z_pos(iw)) + bx_
                by(i,j,inj_z_pos(iw)) = by(i,j,inj_z_pos(iw)) + by_
              else ! injection at a new position, simply replace with the injection value
                bx(i,j,inj_z_pos(iw)) = bx_
                by(i,j,inj_z_pos(iw)) = by_
              endif 
            enddo
          enddo
        else ! inj_wave_radius(iw)>0, and selectively apply to x, y
          do j = jb-1, je+1
            do i = 1, nx2
              if ( sqrt((i-nxmax/2.0)**2.0+(j-nymax/2.0)**2.0)<=inj_wave_radius(iw) ) then
                ! add radial envelope
                radial_env = cos(0.5*pi*(i-nxmax/2)/inj_wave_radius(iw))*cos(0.5*pi*(j-nymax/2)/inj_wave_radius(iw))
                ! add injection value to previous wave if they have the same injection position
                if ( iw>1 .and. inj_z_pos(iw)==inj_z_pos(iw-1) ) then
                  bx(i,j,inj_z_pos(iw)) = bx(i,j,inj_z_pos(iw)) + bx_*radial_env
                  by(i,j,inj_z_pos(iw)) = by(i,j,inj_z_pos(iw)) + by_*radial_env
                else ! injection at a new position, simply replace with the injection value
                  bx(i,j,inj_z_pos(iw)) = bx_*radial_env
                  by(i,j,inj_z_pos(iw)) = by_*radial_env
                endif 
              endif 
            enddo
          enddo
        endif ! end inj_wave_radius

      endif ! end inj_dB_B0(iw)

    enddo ! end wave indexing (iw)

  end subroutine inject_waves_b


  !---------------------------------------------------------------------
  ! inject waves via B & V field
  !---------------------------------------------------------------------
  subroutine inject_waves_bv
    real*8 :: B0, VA, mi
    real*8 :: kx, ky, kz, kxmin, kymin, kzmin
    real*8 :: inj_time
    real*8 :: bx_, by_, bz_
    real*8 :: dvx_, dvy_, dvz_
    integer :: i, j, k
    integer :: iw  ! index of wave
    real*8 :: time_env, radial_env  ! wave envelope

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
    ! kz = wave_cycles*kzmin
    inj_time = it * dtwci  ! in 1/wci

    do iw = 1, 4 
      bx_ = 0.0; by_ = 0.0

      if ( inj_dB_B0(iw)>0.0 .and. (kb-1)<=inj_z_pos(iw) .and. inj_z_pos(iw)<=ke+1 ) then 

        if (inj_time <=inj_t_upramp(iw)) then
          time_env = (sin(0.5*pi*inj_time/inj_t_upramp(iw)))**2.
        else if (inj_time <= inj_t_upramp(iw)+inj_t_flat(iw)) then
          time_env = 1.0
        else if (inj_time <= inj_t_upramp(iw)+inj_t_flat(iw)+inj_t_downramp(iw)) then
          time_env = (cos(0.5*pi*(inj_time-inj_t_upramp(iw)-inj_t_flat(iw))/inj_t_downramp(iw)))**2.
        else 
          time_env = 0.0
        endif

        kz = inj_wave_cycles(iw) * kzmin
        if (inj_wave_pol(iw)==0) then  ! x-pol
          bx_ = inj_dB_B0(iw)*B0*time_env*sin(kz*inj_time)
        else if (inj_wave_pol(iw)==1) then ! y-pol, left-hand pol
          by_ = inj_dB_B0(iw)*B0*time_env*cos(kz*inj_time)
        else if (inj_wave_pol(iw)==-1) then ! y-pol, right-hand pol
          by_ = -inj_dB_B0(iw)*B0*time_env*cos(kz*inj_time)
        endif

        dvx_ = -VA*bx_/B0 * inj_sign_cos(iw) 
        dvy_ = -VA*by_/B0 * inj_sign_cos(iw)

        if (inj_wave_radius(iw)==0) then ! inject at all x, y
          do j = jb-1, je+1
            do i = 1, nx2
              ! add injection value to previous wave if they have the same injection position
              if ( iw>1 .and. inj_z_pos(iw)==inj_z_pos(iw-1) ) then
                bx(i,j,inj_z_pos(iw)) = bx(i,j,inj_z_pos(iw)) + bx_
                by(i,j,inj_z_pos(iw)) = by(i,j,inj_z_pos(iw)) + by_
                ! use vix to temporarily store values of V on the grid
                vix(i,j,inj_z_pos(iw)) = vix(i,j,inj_z_pos(iw)) + dvx_
                viy(i,j,inj_z_pos(iw)) = viy(i,j,inj_z_pos(iw)) + dvy_
              else ! injection at a new position, simply replace with the injection value
                bx(i,j,inj_z_pos(iw)) = bx_
                by(i,j,inj_z_pos(iw)) = by_
                ! use vix to temporarily store values of V on the grid
                vix(i,j,inj_z_pos(iw)) = dvx_
                viy(i,j,inj_z_pos(iw)) = dvy_
              endif 
            enddo
          enddo
        else ! inj_wave_radius(iw)>0, and selectively apply to x, y
          do j = jb-1, je+1
            do i = 1, nx2
              if ( sqrt((i-nxmax/2.0)**2.0+(j-nymax/2.0)**2.0)<=inj_wave_radius(iw) ) then
                ! add radial envelope
                radial_env = cos(0.5*pi*(i-nxmax/2)/inj_wave_radius(iw))*cos(0.5*pi*(j-nymax/2)/inj_wave_radius(iw))
                ! add injection value to previous wave if they have the same injection position
                if ( iw>1 .and. inj_z_pos(iw)==inj_z_pos(iw-1) ) then
                  bx(i,j,inj_z_pos(iw)) = bx(i,j,inj_z_pos(iw)) + bx_*radial_env
                  by(i,j,inj_z_pos(iw)) = by(i,j,inj_z_pos(iw)) + by_*radial_env
                  ! use vix to temporarily store values of V on the grid
                  vix(i,j,inj_z_pos(iw)) = vix(i,j,inj_z_pos(iw)) + dvx_*radial_env
                  viy(i,j,inj_z_pos(iw)) = viy(i,j,inj_z_pos(iw)) + dvy_*radial_env
                else ! injection at a new position, simply replace with the injection value
                  bx(i,j,inj_z_pos(iw)) = bx_*radial_env
                  by(i,j,inj_z_pos(iw)) = by_*radial_env
                  ! use vix to temporarily store values of V on the grid
                  vix(i,j,inj_z_pos(iw)) = dvx_*radial_env
                  viy(i,j,inj_z_pos(iw)) = dvy_*radial_env
                endif 
              endif 
            enddo
          enddo
        endif ! end inj_wave_radius

      endif ! end inj_dB_B0(iw)

    enddo ! end wave indexing (iw)

  end subroutine inject_waves_bv


  !---------------------------------------------------------------------
  ! inject waves via E field
  !---------------------------------------------------------------------
  subroutine inject_waves_e
    real*8 :: B0, VA, mi
    real*8 :: kx, ky, kz, kxmin, kymin, kzmin
    real*8 :: inj_time
    real*8 :: ex_, ey_, ez_
    integer :: i, j, k
    integer :: iw  ! index of wave
    real*8 :: time_env, radial_env  ! wave envelope

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
    inj_time = it * dtwci  ! in 1/wci

    do iw = 1, 4 
      ex_ = 0.0; ey_ = 0.0

      if ( inj_dB_B0(iw)>0.0 .and. (kb-1)<=inj_z_pos(iw) .and. inj_z_pos(iw)<=ke+1 ) then 

        if (inj_time <=inj_t_upramp(iw)) then
          time_env = (sin(0.5*pi*inj_time/inj_t_upramp(iw)))**2.
        else if (inj_time <= inj_t_upramp(iw)+inj_t_flat(iw)) then
          time_env = 1.0
        else if (inj_time <= inj_t_upramp(iw)+inj_t_flat(iw)+inj_t_downramp(iw)) then
          time_env = (cos(0.5*pi*(inj_time-inj_t_upramp(iw)-inj_t_flat(iw))/inj_t_downramp(iw)))**2.
        else 
          time_env = 0.0
        endif

        kz = inj_wave_cycles(iw) * kzmin
        if (inj_wave_pol(iw)==0) then  ! x-pol
          ex_ = inj_dB_B0(iw)*B0/wpiwci*time_env*sin(kz*inj_time)
        else if (inj_wave_pol(iw)==1) then ! y-pol, left-hand pol
          ey_ = inj_dB_B0(iw)*B0/wpiwci*time_env*cos(kz*inj_time)
        else if (inj_wave_pol(iw)==-1) then ! y-pol, right-hand pol
          ey_ = -inj_dB_B0(iw)*B0/wpiwci*time_env*cos(kz*inj_time)
        endif

        if (inj_wave_radius(iw)==0) then ! inject at all x, y
          do j = jb-1, je+1
            do i = 1, nx2
              ! add injection value to previous wave if they have the same injection position
              if ( iw>1 .and. inj_z_pos(iw)==inj_z_pos(iw-1) ) then
                ex(i,j,inj_z_pos(iw)) = ex(i,j,inj_z_pos(iw)) + ex_
                ey(i,j,inj_z_pos(iw)) = ey(i,j,inj_z_pos(iw)) + ey_     
              else ! injection at a new position, simply replace with the injection value
                ex(i,j,inj_z_pos(iw)) = ex_
                ey(i,j,inj_z_pos(iw)) = ey_
              endif 
            enddo
          enddo
        else ! inj_wave_radius(iw)>0, and selectively apply to x, y
          do j = jb-1, je+1
            do i = 1, nx2
              if ( sqrt((i-nxmax/2.0)**2.0+(j-nymax/2.0)**2.0)<=inj_wave_radius(iw) ) then
                ! add radial envelope
                radial_env = cos(0.5*pi*(i-nxmax/2)/inj_wave_radius(iw))*cos(0.5*pi*(j-nymax/2)/inj_wave_radius(iw))
                ! add injection value to previous wave if they have the same injection position
                if ( iw>1 .and. inj_z_pos(iw)==inj_z_pos(iw-1) ) then
                  ex(i,j,inj_z_pos(iw)) = ex(i,j,inj_z_pos(iw)) + ex_*radial_env
                  ey(i,j,inj_z_pos(iw)) = ey(i,j,inj_z_pos(iw)) + ey_*radial_env
                else ! injection at a new position, simply replace with the injection value
                  ex(i,j,inj_z_pos(iw)) = ex_*radial_env
                  ey(i,j,inj_z_pos(iw)) = ey_*radial_env
                endif 
              endif 
            enddo
          enddo
        endif ! end inj_wave_radius
      endif ! end inj_dB_B0(iw)
    enddo ! end wave indexing (iw)

  end subroutine inject_waves_e

end module m_injection