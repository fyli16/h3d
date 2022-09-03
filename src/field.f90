module m_field
  use m_parameter
  use m_mesh
  use m_injection
  implicit none 

  contains 

  !---------------------------------------------------------------------
  ! advances electromagnetic field in time
  !---------------------------------------------------------------------
  subroutine update_fields
    if (ndim /=1) then 
      call field_3d
    else
      call field_2d
    endif 
  end subroutine update_fields


  !---------------------------------------------------------------------
  ! advances electromagnetic field in time (3D)
  !---------------------------------------------------------------------
  subroutine field_3d
    call date_and_time(values=time_begin(:,51))
    call pressgrad(1) ! full step at N
    call date_and_time(values=time_end(:,51))
    call add_time(time_begin(1,51), time_end(1,51), time_elapsed(51))

    ! wave injection via B field
    if (inj_waves_b .eqv. .true.) then
      call inject_waves_b  
    else if (inj_waves_bv .eqv. .true.) then
      call inject_waves_bv
    endif 

    call date_and_time(values=time_begin(:,52))
    call bcalc
    call date_and_time(values=time_end(:,52))
    call add_time(time_begin(1,52), time_end(1,52), time_elapsed(52))

    call date_and_time(values=time_begin(:,51))
    call pressgrad(0) ! half step N+1/2
    call date_and_time(values=time_end(:,51))
    call add_time(time_begin(1,51), time_end(1,51), time_elapsed(51))

    call date_and_time(values=time_begin(:,53))
    call ecalc(0)
    call date_and_time(values=time_end(:,53))
    call add_time(time_begin(1,53), time_end(1,53), time_elapsed(53))

    call date_and_time(values=time_begin(:,54))
    call focalc
    call date_and_time(values=time_end(:,54))
    call add_time(time_begin(1,54), time_end(1,54), time_elapsed(54))

    return
  end subroutine field_3d


  !---------------------------------------------------------------------
  ! computes electron pressure gradient
  !---------------------------------------------------------------------
  subroutine pressgrad(iflag)
    integer, intent(in) :: iflag
    integer*8 :: i, j, k
    real*8 :: dena, dxa, dya, dza, a
    
    if (iflag==0) then ! at half-step
      pe = te0*denh**gamma
    else  ! at full step
      pe = te0*den**gamma
    endif

    do k = kb, ke
      do j = jb, je
        do i = 2, nx1
          dena = iflag*den(i,j,k) + (1-iflag)*denh(i,j,k)
          a = one/dena
          dxa=a/(4.*(meshX%dxn(i  )+meshX%dxn(i+1)))
          dya=a/(4.*(meshY%dxn(j+1)+meshY%dxn(j+2))) ! integer index in y direction starts at 0
          dza=a/(4.*(meshZ%dxn(k+1)+meshZ%dxn(k+2))) ! integer index in z direction starts at 0

          dpedx(i,j,k) = ( (pe(i+1,j-1,k+1)+2.*pe(i+1,j,k+1) + pe(i+1,j+1,k+1))/4. &
                + 2.*(pe(i+1,j-1,k  )+2.*pe(i+1,j,k  )+pe(i+1,j+1,k  ))/4. &
                + (pe(i+1,j-1,k-1)+2.*pe(i+1,j,k-1)+pe(i+1,j+1,k-1))/4. &
                - (pe(i-1,j-1,k+1)+2.*pe(i-1,j,k+1)+pe(i-1,j+1,k+1))/4. &
                - 2.*(pe(i-1,j-1,k  )+2.*pe(i-1,j,k  )+pe(i-1,j+1,k))/4. &
                - (pe(i-1,j-1,k-1)+2.*pe(i-1,j,k-1)+pe(i-1,j+1,k-1))/4. ) * dxa

          dpedy(i,j,k)=( (pe(i-1,j+1,k+1)+2.*pe(i,j+1,k+1) + pe(i+1,j+1,k+1))/4.&
                + 2.*(pe(i-1,j+1,k  )+2.*pe(i,j+1,k  )+pe(i+1,j+1,k))/4. &
                +   (pe(i-1,j+1,k-1)+2.*pe(i,j+1,k-1)+pe(i+1,j+1,k-1))/4. &
                -   (pe(i-1,j-1,k+1)+2.*pe(i,j-1,k+1)+pe(i+1,j-1,k+1))/4. &
                -2.*(pe(i-1,j-1,k  )+2.*pe(i,j-1,k  )+pe(i+1,j-1,k))/4. &
                -   (pe(i-1,j-1,k-1)+2.*pe(i,j-1,k-1)+pe(i+1,j-1,k-1))/4. ) * dya

          dpedz(i,j,k)=( (pe(i+1,j-1,k+1)+2.*pe(i+1,j,k+1) + pe(i+1,j+1,k+1))/4. & 
                + 2.*(pe(i  ,j-1,k+1)+2.*pe(i,j,k+1)+pe(i,j+1,k+1))/4. &
                + (pe(i-1,j-1,k+1)+2.*pe(i-1,j,k+1)+pe(i-1,j+1,k+1))/4. &
                - (pe(i+1,j-1,k-1)+2.*pe(i+1,j,k-1)+pe(i+1,j+1,k-1))/4. &
                - 2.*(pe(i  ,j-1,k-1)+2.*pe(i,j,k-1)+pe(i,j+1,k-1))/4. &
                - (pe(i-1,j-1,k-1)+2.*pe(i-1,j,k-1)+pe(i-1,j+1,k-1))/4. ) * dza                    
        enddo
      enddo
    enddo

    return
  end subroutine pressgrad


  !---------------------------------------------------------------------
  ! computes electric field and curl(E)
  !---------------------------------------------------------------------
  subroutine ecalc(iflag)
    integer, intent(in) :: iflag
    integer*8 :: i, j, k
    real*8 :: term1, term2, term3, term4, fm
    real*8 :: tenx,teny,tenz,xj,yj,zj,bxx,byy,bzz,btot,tjdotb,curr_tot
    real*8 :: vixa, viya, viza, dena, a, dxa, dya, dza 
    real*8 :: dbxdy, dbxdz, dbydx, dbydz, dbzdx, dbzdy 
    real*8 :: curlbx_scalar,curlby_scalar,curlbz_scalar
    real*8 :: bxav, byav, bzav  
    real*8 :: dexdy, dexdz, deydx, deydz, dezdx, dezdy  

    do k = kb, ke
      ! fm = masking_func(k, iflag)
      do j = jb, je
        do i = 2, nx1
          vixa = (1.-iflag)*(1.5*vix(i,j,k)-0.5*vixo(i,j,k)) + iflag*vix(i,j,k)
          viya = (1.-iflag)*(1.5*viy(i,j,k)-0.5*viyo(i,j,k)) + iflag*viy(i,j,k)
          viza = (1.-iflag)*(1.5*viz(i,j,k)-0.5*vizo(i,j,k)) + iflag*viz(i,j,k)

          ! iflag=1 (using density at full step) used by ecal in bcal
          ! iflag=0 (using density at half step) used by final ecal 
          !                                 [advancing E^n to E^(n+1)]
          dena = iflag*den(i,j,k) + (1-iflag)*denh(i,j,k)
          ! dena = iflag*0.5*(den(i,j,k)+deno(i,j,k)) + (1.-iflag)*den(i,j,k)
          a = one/dena

          dxa = a/(4.*meshX%dxc(i))
          dya = a/(4.*meshY%dxc(j+1)) ! index in y direction starts at 0
          dza = a/(4.*meshZ%dxc(k+1)) ! index in z direction starts at 0

          dbxdy = bx(i+1,j+1,k+1) + bx(i,j+1,k+1) + bx(i,j+1,k) + bx(i+1,j+1,k) &
                - bx(i+1,j,  k+1) - bx(i,j,  k+1) - bx(i,j,  k) - bx(i+1,j,  k)
          dbxdz = bx(i+1,j+1,k+1) + bx(i,j+1,k+1) + bx(i,j,k+1) + bx(i+1,j,k+1) &
                - bx(i+1,j+1,k  ) - bx(i,j+1,k  ) - bx(i,j,k  ) - bx(i+1,j,k  )
          dbydx = by(i+1,j+1,k+1)+by(i+1,j  ,k+1)+by(i+1,j  ,k  )+by(i+1,j+1,k  )&
                -by(i  ,j+1,k+1)-by(i  ,j  ,k+1)-by(i  ,j  ,k  )-by(i  ,j+1,k  )
          dbydz= by(i+1,j+1,k+1)+by(i  ,j+1,k+1)+by(i  ,j  ,k+1)+by(i+1,j  ,k+1)&
                -by(i+1,j+1,k  )-by(i  ,j+1,k  )-by(i  ,j  ,k  )-by(i+1,j  ,k  )
          dbzdx= bz(i+1,j+1,k+1)+bz(i+1,j  ,k+1)+bz(i+1,j  ,k  )+bz(i+1,j+1,k  )&
                -bz(i  ,j+1,k+1)-bz(i  ,j  ,k+1)-bz(i  ,j  ,k  )-bz(i  ,j+1,k  )
          dbzdy= bz(i+1,j+1,k+1)+bz(i  ,j+1,k+1)+bz(i  ,j+1,k  )+bz(i+1,j+1,k  )&
                -bz(i+1,j  ,k+1)-bz(i  ,j  ,k+1)-bz(i  ,j  ,k  )-bz(i+1,j  ,k  )

          curlbx_scalar = dya*dbzdy - dza*dbydz
          curlby_scalar = dza*dbxdz - dxa*dbzdx
          curlbz_scalar = dxa*dbydx - dya*dbxdy

          bxav = 0.125*( bx(i+1,j+1,k) + bx(i,j+1,k) + bx(i,j,k) + bx(i+1,j ,k) &
                  + bx(i+1,j+1,k+1) + bx(i,j+1,k+1) + bx(i,j,k+1) + bx(i+1,j,k+1) )
          byav = 0.125*( by(i+1,j+1,k) + by(i,j+1,k) + by(i,j,k) + by(i+1,j ,k) &
                  + by(i+1,j+1,k+1) + by(i,j+1,k+1) + by(i,j,k+1) + by(i+1,j,k+1) )
          bzav = 0.125*( bz(i+1,j+1,k) + bz(i,j+1,k) + bz(i,j,k) + bz(i+1,j ,k) &
                  + bz(i+1,j+1,k+1) + bz(i,j+1,k+1) + bz(i,j,k+1) + bz(i+1,j,k+1) )

          xj = curlbx_scalar
          yj = curlby_scalar
          zj = curlbz_scalar

          if (eta_par==0) then
            tenx = eta(i,j,k)*xj
            teny = eta(i,j,k)*yj
            tenz = eta(i,j,k)*zj
          else 
            bxx = bxav
            byy = byav
            bzz = bzav
            btot = sqrt(bxx**2 + byy**2 + bzz**2)
            if (eta_par==1) then     
              tjdotb = eta(i,j,k)*(bxx*xj + byy*yj + bzz*zj)/btot
              tenx = tjdotb*bxx/btot
              teny = tjdotb*byy/btot
              tenz = tjdotb*bzz/btot
            else if (eta_par==2) then
              curr_tot = max(1.d-12,sqrt(xj**2 + yj**2 + zj**2))
              tenx = abs(eta(i,j,k)*bxx*xj/(btot*curr_tot))
              tenx = min(resis,tenx)
              tenx = tenx*xj
              teny = abs(eta(i,j,k)*byy*yj/(btot*curr_tot))
              teny = min(resis,teny)
              teny = teny*yj
              tenz = abs(eta(i,j,k)*bzz*zj/(btot*curr_tot))
              tenz = min(resis,tenz)
              tenz = tenz*zj
            endif 
          endif 

          ex(i,j,k) = (viza*byav-viya*bzav) + (curlby_scalar*bzav-curlbz_scalar*byav) - dpedx(i,j,k) + tenx/a
          ey(i,j,k) = (vixa*bzav-viza*bxav) + (curlbz_scalar*bxav-curlbx_scalar*bzav) - dpedy(i,j,k) + teny/a
          ez(i,j,k) = (viya*bxav-vixa*byav) + (curlbx_scalar*byav-curlby_scalar*bxav) - dpedz(i,j,k) + tenz/a 

          ! ex(i,j,k) = (viza*byav-viya*bzav)*fm + (curlby_scalar*bzav-curlbz_scalar*byav)*fm - dpedx(i,j,k) + tenx/a
          ! ey(i,j,k) = (vixa*bzav-viza*bxav)*fm + (curlbz_scalar*bxav-curlbx_scalar*bzav)*fm - dpedy(i,j,k) + teny/a
          ! ez(i,j,k) = (viya*bxav-vixa*byav)*fm + (curlbx_scalar*byav-curlby_scalar*bxav)*fm - dpedz(i,j,k) + tenz/a 
          
          ! if (myid == 0) then
          !   print*, 'viya*bxav-vixa*byav = ', viya*bxav-vixa*byav
          !   print*, 'urlbx_scalar*byav-curlby_scalar*bxav = ', curlbx_scalar*byav-curlby_scalar*bxav
          !   print*, 'dpedz(i,j,k) =', dpedz(i,j,k)
          !   print*, 'tenz/a = ', tenz/a
          ! endif 

          ! if (myid==0) then
          !   if (viya*bxav-vixa*byav>1e-4/wpiwci**2 .or. curlbx_scalar*byav-curlby_scalar*bxav>1e-4/wpiwci**2 &
          !   .or. dpedz(i,j,k)>1e-4/wpiwci**2 .or. tenz/a>1e-4/wpiwci**2) then
          !     print*, 'it,i,j,k = ', it,i,j,k
          !     print*, 'viya*bxav-vixa*byav = ', viya*bxav-vixa*byav
          !     print*, 'curlbx_scalar*byav-curlby_scalar*bxav = ', curlbx_scalar*byav-curlby_scalar*bxav
          !     print*, 'dpedz(i,j,k) =', dpedz(i,j,k)
          !     print*, 'tenz/a = ', tenz/a
          !     print*, 'ez(i,j,k)*(wpiwci**2) =', ez(i,j,k)*wpiwci**2
          !     print*, " "
          !   endif 
          ! endif 

          ! if (i==2 .and. j==jb .and. k==kb) then
          !   term1 = viya*bxav-vixa*byav
          !   term2 = curlbx_scalar*byav-curlby_scalar*bxav
          !   term3 = - dpedz(i,j,k)
          !   term4 = tenz/a
          ! endif 

        enddo
      enddo
    enddo

    ! debug ez component (only when iflag==0, i.e., outside 'bcal' call)
    ! if (n_debug_ez > 0 .and. mod(it, n_debug_ez) ==0 .and. iflag==0) then
    !   write(int(100), '(I6,1x,4E14.6,1x)') it, &
    !       term1, term2, term3, term4
    ! endif 

    ! wave injection via E field
    ! only inject E field after updating B field (i.e., iflag=0)
    if (iflag>0 .and. inj_waves_e .eqv. .true.) then
      call inject_waves_e
    endif 

    ! boundary conditions
    call date_and_time(values=time_begin(:,18))
    call xrealbcc_pack_e(ex,ey,ez,1_8,nx,ny,nz)
    call date_and_time(values=time_end(:,18))
    call add_time(time_begin(1,18),time_end(1,18),time_elapsed(18))

    ! calculate curl E
    do k = kb, ke+1
      do j = jb, je+1
        do i = 2, nx2
          dexdy=  ex(i  ,j  ,k  ) + ex(i-1,j  ,k  ) + ex(i-1,j  ,k-1) + ex(i  ,j  ,k-1)   &
                - ex(i  ,j-1,k  ) - ex(i-1,j-1,k  ) - ex(i-1,j-1,k-1) - ex(i  ,j-1,k-1)
          dexdz=  ex(i  ,j  ,k  ) + ex(i-1,j  ,k  ) + ex(i-1,j-1,k  ) + ex(i  ,j-1,k  )   &
                - ex(i  ,j  ,k-1) - ex(i-1,j  ,k-1) - ex(i-1,j-1,k-1) - ex(i  ,j-1,k-1)
          deydx=  ey(i  ,j  ,k  ) + ey(i  ,j-1,k  ) + ey(i  ,j-1,k-1) + ey(i  ,j  ,k-1)   &
                - ey(i-1,j  ,k  ) - ey(i-1,j-1,k  ) - ey(i-1,j-1,k-1) - ey(i-1,j  ,k-1)
          deydz=  ey(i  ,j  ,k  ) + ey(i-1,j  ,k  ) + ey(i-1,j-1,k  ) + ey(i  ,j-1,k  )   &
                - ey(i  ,j  ,k-1) - ey(i-1,j  ,k-1) - ey(i-1,j-1,k-1) - ey(i  ,j-1,k-1)
          dezdx=  ez(i  ,j  ,k  ) + ez(i  ,j-1,k  ) + ez(i  ,j-1,k-1) + ez(i  ,j  ,k-1)   &
                - ez(i-1,j  ,k  ) - ez(i-1,j-1,k  ) - ez(i-1,j-1,k-1) - ez(i-1,j  ,k-1)
          dezdy=  ez(i  ,j  ,k  ) + ez(i-1,j  ,k  ) + ez(i-1,j  ,k-1) + ez(i  ,j  ,k-1)   &
                - ez(i  ,j-1,k  ) - ez(i-1,j-1,k  ) - ez(i-1,j-1,k-1) - ez(i  ,j-1,k-1)

          curlex(i,j,k) = dezdy/(4.*meshY%dxn(j+1)) - deydz/(4.*meshZ%dxn(k+1))  ! index in y, z start  at 0
          curley(i,j,k) = dexdz/(4.*meshZ%dxn(k+1)) - dezdx/(4.*meshX%dxn(i  ))  ! index in z starts at 0
          curlez(i,j,k) = deydx/(4.*meshX%dxn(i  )) - dexdy/(4.*meshY%dxn(j+1))  ! index in y starts at 0
          
        enddo
      enddo
    enddo

    return
  end subroutine ecalc


  !---------------------------------------------------------------------
  ! advance magnetic field
  !---------------------------------------------------------------------
  subroutine bcalc
    integer*8 :: i, j, k
    integer :: ii
    real*8 :: dts, dts2, dts6
    real*8 :: fm
    real*8 :: tempx1(nxmax,jb-1:je+1,kb-1:ke+1) &
             ,tempy1(nxmax,jb-1:je+1,kb-1:ke+1) &
             ,tempz1(nxmax,jb-1:je+1,kb-1:ke+1)
    
    dts  = dt/real(n_sub_b)
    dts2 = dts/2.
    dts6 = dts/6.

    ! subcycle into n_sub_b interations
    !VR : E is synchronized between processors at each step of RK.
    !VR : but is it enough to make sure that B is consistent?
    do ii = 1, n_sub_b
      bxs=bx; bys=by; bzs=bz ! save B at the start of each subcycle

      ! R-K part 1
      call ecalc(1)  
      endif 
      do k = kb, ke+1
        do j = jb, je+1
          do i = 2, nx2
            ! B = B(n)+dt*K1/2
            bx(i,j,k) = bxs(i,j,k) - dts2*curlex(i,j,k)
            by(i,j,k) = bys(i,j,k) - dts2*curley(i,j,k)
            bz(i,j,k) = bzs(i,j,k) - dts2*curlez(i,j,k)
            ! temp1 = -K1
            tempx1(i,j,k) = curlex(i,j,k)
            tempy1(i,j,k) = curley(i,j,k)
            tempz1(i,j,k) = curlez(i,j,k)
          enddo
        enddo
      enddo   
        
      ! R-K part 2
      call ecalc(1)
      do k = kb, ke+1
        do j = jb, je+1
          do i = 2, nx2
            ! B = B(n)+dt*K2/2
            bx(i,j,k) = bxs(i,j,k) - dts2*curlex(i,j,k)
            by(i,j,k) = bys(i,j,k) - dts2*curley(i,j,k)
            bz(i,j,k) = bzs(i,j,k) - dts2*curlez(i,j,k)
            ! temp1 = -K1 - 2*K2
            tempx1(i,j,k) = tempx1(i,j,k) + 2.*curlex(i,j,k)
            tempy1(i,j,k) = tempy1(i,j,k) + 2.*curley(i,j,k)
            tempz1(i,j,k) = tempz1(i,j,k) + 2.*curlez(i,j,k)
          enddo
        enddo
      enddo
        
      ! R-K part 3
      call ecalc(1)
      do k = kb, ke+1
        do j = jb, je+1
          do i = 2, nx2
            ! B = B(n)+dt*K3
            bx(i,j,k) = bxs(i,j,k) - dts*curlex(i,j,k)
            by(i,j,k) = bys(i,j,k) - dts*curley(i,j,k)
            bz(i,j,k) = bzs(i,j,k) - dts*curlez(i,j,k)
            ! temp1 = -K1 - 2*K2 - 2*K3
            tempx1(i,j,k) = tempx1(i,j,k) + 2.*curlex(i,j,k)
            tempy1(i,j,k) = tempy1(i,j,k) + 2.*curley(i,j,k)
            tempz1(i,j,k) = tempz1(i,j,k) + 2.*curlez(i,j,k)
          enddo
        enddo
      enddo
        
      ! R-K part 4
      call ecalc(1)
      do k = kb, ke+1
        fm = masking_func(k, ii)
        do j = jb, je+1
          do i = 2, nx2
            ! B = B(n) - (dt/6)*(-K1 - 2*K2 - 2*K3 - K4)
            ! bx(i,j,k) = bxs(i,j,k) - dts6*(tempx1(i,j,k)+curlex(i,j,k))
            ! by(i,j,k) = bys(i,j,k) - dts6*(tempy1(i,j,k)+curley(i,j,k))
            ! bz(i,j,k) = bzs(i,j,k) - dts6*(tempz1(i,j,k)+curlez(i,j,k))
            bx(i,j,k) = ( bxs(i,j,k) - dts6*(tempx1(i,j,k)+curlex(i,j,k)) )*fm
            by(i,j,k) = ( bys(i,j,k) - dts6*(tempy1(i,j,k)+curley(i,j,k)) )*fm
            bz(i,j,k) = bzs(i,j,k) - dts6*(tempz1(i,j,k)+curlez(i,j,k)) 
          enddo
        enddo
      enddo

      ! end of iteration loop
      CALL xrealbcc_pack_b(bx,by,bz,1_8,nx,ny,nz) 

    end do

    return
  end subroutine bcalc


  !---------------------------------------------------------------------
  ! compute friction force due to resistivity?
  !---------------------------------------------------------------------
  subroutine focalc
    integer*8 :: i,j,k
    real*8 :: tenx,teny,tenz,xj,yj,zj,bxx,byy,bzz,btot,tjdotb,curr_tot
    real*8 :: dbxdy,dbydx,dbzdx,dbxdz,dbzdy,dbydz
    real*8 :: curlbx_scalar,curlby_scalar,curlbz_scalar,bxav,byav,bzav

    do k = kb, ke 
      do j = jb, je
        do i = 2, nx1
          dbxdy= bx(i,j,k)+bx(i-1,j,k)+bx(i-1,j,k-1)+bx(i,j,k-1)&
                -bx(i,j-1,k)-bx(i-1,j-1,k)-bx(i-1,j-1,k-1)-bx(i,j-1,k-1)
          dbxdz= bx(i,j,k  )+bx(i-1,j,k  )+bx(i-1,j-1,k  )+bx(i,j-1,k  )&
                -bx(i,j,k-1)-bx(i-1,j,k-1)-bx(i-1,j-1,k-1)-bx(i,j-1,k-1)
          dbydx= by(i  ,j,k)+by(i  ,j-1,k)+by(i  ,j-1,k-1)+by(i  ,j,k-1)&
                -by(i-1,j,k)-by(i-1,j-1,k)-by(i-1,j-1,k-1)-by(i-1,j,k-1)
          dbydz= by(i,j,k  )+by(i-1,j,k  )+by(i-1,j-1,k  )+by(i,j-1,k  )&
                -by(i,j,k-1)-by(i-1,j,k-1)-by(i-1,j-1,k-1)-by(i,j-1,k-1)
          dbzdx= bz(i  ,j,k)+bz(i  ,j-1,k)+bz(i  ,j-1,k-1)+bz(i  ,j,k-1)&
                -bz(i-1,j,k)-bz(i-1,j-1,k)-bz(i-1,j-1,k-1)-bz(i-1,j,k-1)
          dbzdy= bz(i,j  ,k)+bz(i-1,j  ,k)+bz(i-1,j  ,k-1)+bz(i,j  ,k-1)&
                -bz(i,j-1,k)-bz(i-1,j-1,k)-bz(i-1,j-1,k-1)-bz(i,j-1,k-1)

          curlbx_scalar=dbzdy/(4.*meshY%dxc(j+1))-dbydz/(4.*meshZ%dxc(k+1))
          curlby_scalar=dbxdz/(4.*meshZ%dxc(k+1))-dbzdx/(4.*meshX%dxc(i  ))
          curlbz_scalar=dbydx/(4.*meshX%dxc(i  ))-dbxdy/(4.*meshY%dxc(j+1))

          ! 6/25/2006 New eta_par option: tensor eta
          bxav = 0.125*( bx(i+1,j+1,k) + bx(i,j+1,k) + bx(i,j,k) + bx(i+1,j ,k) &
                  + bx(i+1,j+1,k+1) + bx(i,j+1,k+1) + bx(i,j,k+1) + bx(i+1,j,k+1) )
          byav = 0.125*( by(i+1,j+1,k) + by(i,j+1,k) + by(i,j,k) + by(i+1,j ,k) &
                  + by(i+1,j+1,k+1) + by(i,j+1,k+1) + by(i,j,k+1) + by(i+1,j,k+1) )
          bzav = 0.125*( bz(i+1,j+1,k) + bz(i,j+1,k) + bz(i,j,k) + bz(i+1,j ,k) &
                  + bz(i+1,j+1,k+1) + bz(i,j+1,k+1) + bz(i,j,k+1) + bz(i+1,j,k+1) )
          
          xj = curlbx_scalar
          yj = curlby_scalar
          zj = curlbz_scalar

          tenx = eta(i,j,k)*xj
          teny = eta(i,j,k)*yj
          tenz = eta(i,j,k)*zj

          fox(i,j,k) = -tenx
          foy(i,j,k) = -teny
          foz(i,j,k) = -tenz
        enddo
      enddo
    enddo

    ! boundary conditions
    ! first update internal ghost cells so that all processors
    ! have latest information. Care must be exercised so that
    ! BCs are set wrt proecessor that corresponds to that location. 
    ! To that end, keep z loops set to limits of kb and ke.
    call xrealbcc_pack_e(fox,foy,foz,1_8,nx,ny,nz)

    return
  end subroutine focalc


  !---------------------------------------------------------------------
  ! advances electromagnetic field in time (2D)
  !---------------------------------------------------------------------
  subroutine field_2d
    call date_and_time(values=time_begin(:,51))
    call pressgrad_2d(1)
    call date_and_time(values=time_end(:,51))
    call add_time(time_begin(1,51),time_end(1,51),time_elapsed(51))

    call date_and_time(values=time_begin(:,52))
    call bcalc_2d
    call date_and_time(values=time_end(:,52))
    call add_time(time_begin(1,52),time_end(1,52),time_elapsed(52))

    call date_and_time(values=time_begin(:,51))
    call pressgrad_2d(0)
    call date_and_time(values=time_end(:,51))
    call add_time(time_begin(1,51),time_end(1,51),time_elapsed(51))

    call date_and_time(values=time_begin(:,53))
    call ecalc_2d( 0)
    call date_and_time(values=time_end(:,53))
    call add_time(time_begin(1,53),time_end(1,53),time_elapsed(53))

    call date_and_time(values=time_begin(:,54))
    call focalc_2d
    call date_and_time(values=time_end(:,54))
    call add_time(time_begin(1,54),time_end(1,54),time_elapsed(54))

    return
  end subroutine field_2d


  !---------------------------------------------------------------------
  subroutine pressgrad_2d(iflag)
    integer :: iflag
    integer*8 :: i,j,k
    real*8 :: dena,dxa,dya,dza,a

    do k=kb,ke
      do j = jb,je
        do i=2,nx1
          dena = iflag*0.5*(den(i,j,k)+deno(i,j,k)) + (1.-iflag)*den(i,j,k)
          a=1/dena

          dxa=a/((meshX%dxn(i  )+meshX%dxn(i+1)))
          dya=a/((meshY%dxn(j+1)+meshY%dxn(j+2)))  ! integer index in y direction starts at 0
          dza=a/((meshZ%dxn(k+1)+meshZ%dxn(k+2)))  ! integer index in z direction starts at 0

          dpedx(i,j,k)=(  (pe(i+1,j-1,k  )+2.*pe(i+1,j,k  )+pe(i+1,j+1,k  ))/4.    &
                        - (pe(i-1,j-1,k  )+2.*pe(i-1,j,k  )+pe(i-1,j+1,k  ))/4. )*dxa
          dpedy(i,j,k)=(  (pe(i-1,j+1,k  )+2.*pe(i,j+1,k  )+pe(i+1,j+1,k  ))/4.    &
                        - (pe(i-1,j-1,k  )+2.*pe(i,j-1,k  )+pe(i+1,j-1,k  ))/4. )*dya
          dpedz(i,j,k)=0.0
        enddo
      enddo
    enddo
    return
  end subroutine pressgrad_2d


  !---------------------------------------------------------------------
  subroutine ecalc_2d(iflag)
    integer :: iflag
    integer*8 :: i,j,k

    real*8 :: tenx,teny,tenz,xj,yj,zj,bxx,byy,bzz,btot,tjdotb &
                      ,curr_tot
    real*8 :: bx1,bx2,bx3,bx4,bx5,bx6,bx7,bx8 
    real*8 :: by1,by2,by3,by4,by5,by6,by7,by8 
    real*8 :: bz1,bz2,bz3,bz4,bz5,bz6,bz7,bz8  
    real*8 :: vixa, viya, viza, dena, a, dxa, dya, dza 
    real*8 :: dbxdy, dbxdz, dbydx, dbydz, dbzdx, dbzdy 
    real*8 :: curlbx_scalar,curlby_scalar,curlbz_scalar
    real*8 :: bxav, byav, bzav  
    real*8 :: dexdy, dexdz, deydx, deydz, dezdx,dezdy  

    do k=kb,ke
      do j = jb,je
        do i=2,nx1
          bx1=bx(i+1,j+1,k)  
          bx2=bx(i  ,j+1,k)  
          bx3=bx(i  ,j  ,k)  
          bx4=bx(i+1,j  ,k)  
          by1=by(i+1,j+1,k)  
          by2=by(i  ,j+1,k)  
          by3=by(i  ,j  ,k)  
          by4=by(i+1,j  ,k)  
          bz1=bz(i+1,j+1,k)  
          bz2=bz(i  ,j+1,k)  
          bz3=bz(i  ,j  ,k)  
          bz4=bz(i+1,j  ,k)  

          vixa=(1.-iflag)*(1.5*vix(i,j,k)-0.5*vixo(i,j,k))+iflag*vix(i,j,k)
          viya=(1.-iflag)*(1.5*viy(i,j,k)-0.5*viyo(i,j,k))+iflag*viy(i,j,k)
          viza=(1.-iflag)*(1.5*viz(i,j,k)-0.5*vizo(i,j,k))+iflag*viz(i,j,k)
          
          dena=iflag*0.5*(den(i,j,k)+deno(i,j,k))+(1.-iflag)*den(i,j,k)
          a=1/dena

          dxa=a/(2.*meshX%dxc(i))
          dya=a/(2.*meshY%dxc(j+1)) ! integer index in y direction starts at 0
          dza=a/(2.*meshZ%dxc(k+1)) ! integer index in z direction starts at 0

          dbxdy = bx(i  ,j+1,k  )+bx(i+1,j+1,k  )-bx(i  ,j  ,k  )-bx(i+1,j  ,k  )
          dbxdz = 0.0
          dbydx = by(i+1,j  ,k  )+by(i+1,j+1,k  )-by(i  ,j  ,k  )-by(i  ,j+1,k  )
          dbydz = 0.0
          dbzdx = bz(i+1,j  ,k  )+bz(i+1,j+1,k  )-bz(i  ,j  ,k  )-bz(i  ,j+1,k  )
          dbzdy = bz(i  ,j+1,k  )+bz(i+1,j+1,k  )-bz(i  ,j  ,k  )-bz(i+1,j  ,k  )

          curlbx_scalar=dya*dbzdy-dza*dbydz
          curlby_scalar=dza*dbxdz-dxa*dbzdx
          curlbz_scalar=dxa*dbydx-dya*dbxdy

          bxav=0.25*(bx1+bx2+bx3+bx4)
          byav=0.25*(by1+by2+by3+by4)
          bzav=0.25*(bz1+bz2+bz3+bz4)

          ! 6/25/2006 New eta_par option: tensor eta
          xj = curlbx_scalar
          yj = curlby_scalar
          zj = curlbz_scalar

          if (eta_par.eq.0) then
            tenx=eta(i,j,k)*xj
            teny=eta(i,j,k)*yj
            tenz=eta(i,j,k)*zj
          else
            bxx = bxav
            byy = byav
            bzz = bzav
            btot = sqrt(bxx**2 + byy**2 + bzz**2)
            xtmp1m = sqrt(xj**2 + yj**2 + zj**2)
            xtmp2m = 1.d-12
            curr_tot = max(xtmp2m,xtmp1m)

            if (eta_par.eq.1) then
              tjdotb = eta(i,j,k)*(bxx*xj + byy*yj + bzz*zj)/btot
              tenx = tjdotb*bxx/btot
              teny = tjdotb*byy/btot
              tenz = tjdotb*bzz/btot
            else if (eta_par.eq.2) then
              xtmp1m = sqrt(xj**2 + yj**2 + zj**2)
              xtmp2m = 1.d-12
              curr_tot = max(xtmp2m,xtmp1m)
              tenx = abs(eta(i,j,k)*bxx*xj/(btot*curr_tot))
              tenx = min(resis,tenx)
              tenx = tenx*xj
              teny = abs(eta(i,j,k)*byy*yj/(btot*curr_tot))
              teny = min(resis,teny)
              teny = teny*yj
              tenz = abs(eta(i,j,k)*bzz*zj/(btot*curr_tot))
              tenz = min(resis,tenz)
              tenz = tenz*zj
            endif
          endif

          ex(i,j,k) = (viza*byav-viya*bzav) + (curlby_scalar*bzav-curlbz_scalar*byav) - dpedx(i,j,k) + tenx/a     
          ey(i,j,k) = (vixa*bzav-viza*bxav) + (curlbz_scalar*bxav-curlbx_scalar*bzav) - dpedy(i,j,k) + teny/a     
          ez(i,j,k) = (viya*bxav-vixa*byav) + (curlbx_scalar*byav-curlby_scalar*bxav) - dpedz(i,j,k) + tenz/a                     
        enddo
      enddo
    enddo
  
    ! boundary conditions
    call date_and_time(values=time_begin(:,18))
    call xrealbcc_pack_e_2d(ex,ey,ez,1_8,nx,ny,nz)
    call date_and_time(values=time_end(:,18))
    call add_time(time_begin(1,18),time_end(1,18),time_elapsed(18))

    !VR: boundaries in x & z
    ex(nx2,:,:)=ex(2,:,:)
    ey(nx2,:,:)=ey(2,:,:)
    ez(nx2,:,:)=ez(2,:,:)

    ex(1,:,:)=ex(nx1,:,:)
    ey(1,:,:)=ey(nx1,:,:)
    ez(1,:,:)=ez(nx1,:,:)

    ex(:,:,kb-1)=ex(:,:,kb)
    ey(:,:,kb-1)=ey(:,:,kb)
    ez(:,:,kb-1)=ez(:,:,kb)

    ex(:,:,kb+1)=ex(:,:,kb)
    ey(:,:,kb+1)=ey(:,:,kb)
    ez(:,:,kb+1)=ez(:,:,kb)

    ! calculate curl E
    do k=kb,ke+1
      do j = jb,je+1
        do i=2,nx2
          dexdy= ex(i  ,j  ,k)+ex(i-1,j  ,k)&
                -ex(i  ,j-1,k)-ex(i-1,j-1,k)
          dexdz= 0.0
          deydx= ey(i  ,j  ,k)+ey(i  ,j-1,k)&
                -ey(i-1,j  ,k)-ey(i-1,j-1,k)
          deydz= 0.0
          dezdx= ez(i  ,j  ,k)+ez(i  ,j-1,k)&
                -ez(i-1,j  ,k)-ez(i-1,j-1,k)
          dezdy= ez(i  ,j  ,k)+ez(i-1,j  ,k)&
                -ez(i  ,j-1,k)-ez(i-1,j-1,k)

          curlex(i,j,k)=dezdy/(2.*meshY%dxn(j+1))-deydz/(2.*meshZ%dxn(k+1)) ! integer index in y and z directions start  at 0
          curley(i,j,k)=dexdz/(2.*meshZ%dxn(k+1))-dezdx/(2.*meshX%dxn(i  )) ! integer index in z       direction  starts at 0
          curlez(i,j,k)=deydx/(2.*meshX%dxn(i  ))-dexdy/(2.*meshY%dxn(j+1)) ! integer index in y       direction  starts at 0
        enddo
      enddo
    enddo
  
    return
  end subroutine ecalc_2d


  !---------------------------------------------------------------------
  subroutine bcalc_2d
    integer*8 :: i, j, k, ii
    real*8 :: dts, dts2, dts6
    real*8 :: tempx1(nxmax,jb-1:je+1,kb-1:ke+1)&
                      ,tempy1(nxmax,jb-1:je+1,kb-1:ke+1)&
                      ,tempz1(nxmax,jb-1:je+1,kb-1:ke+1)

    dts=dt/real(n_sub_b)
    dts2=dts/2.
    dts6=dts/6.

    ! subcycle into n_sub_b interations
    do ii = 1, n_sub_b
      bxs=bx
      bys=by
      bzs=bz

      ! R-K first part
      call ecalc_2d( 1 )
      ! B = B(n)+dt*K1/2
      bx=bxs-dts2*curlex
      by=bys-dts2*curley
      bz=bzs-dts2*curlez
      ! temp1 = K1
      tempx1=curlex
      tempy1=curley
      tempz1=curlez

      ! R-K part 2
      call ecalc_2d( 1 )
      ! B = B(n)+dt*K2/2
      bx=bxs-dts2*curlex
      by=bys-dts2*curley
      bz=bzs-dts2*curlez
      ! temp2 = K2
      tempx1=tempx1+2.*curlex
      tempy1=tempy1+2.*curley
      tempz1=tempz1+2.*curlez

      ! R-K  part 3
      call ecalc_2d( 1 )
      ! B = B(n)+dt*K3
      bx=bxs-dts*curlex
      by=bys-dts*curley
      bz=bzs-dts*curlez
      ! temp3 = K3
      tempx1=tempx1+2.*curlex
      tempy1=tempy1+2.*curley
      tempz1=tempz1+2.*curlez

      ! R-K  part 4
      call ecalc_2d( 1 )
      ! B = B(n) + dt*(K1+2K2+2K3+K4)/6
      bx=bxs-dts6*(tempx1+curlex)
      by=bys-dts6*(tempy1+curley)
      bz=bzs-dts6*(tempz1+curlez)

      call xrealbcc_pack_b_2d(bx,by,bz,1_8,nx,ny,nz)
      !VR: exchange X and Z info
      bx(1  ,:,:)=bx(nx1  ,:,:)
      by(1  ,:,:)=by(nx1  ,:,:)
      bz(1  ,:,:)=bz(nx1  ,:,:)
      
      bx(nx2  ,:,:)=bx(2  ,:,:)
      by(nx2  ,:,:)=by(2  ,:,:)
      bz(nx2  ,:,:)=bz(2  ,:,:)

      bx(:,:,kb-1)=bx(:,:,kb)
      by(:,:,kb-1)=by(:,:,kb)
      bz(:,:,kb-1)=bz(:,:,kb)
      
      bx(:,:,kb+1)=bx(:,:,kb)
      by(:,:,kb+1)=by(:,:,kb)
      bz(:,:,kb+1)=bz(:,:,kb)
    end do

    return
  end subroutine bcalc_2d


  !---------------------------------------------------------------------
  subroutine focalc_2d
    real*8 :: bx1,bx2,bx3,bx4,bx5,bx6,bx7,bx8
    real*8 :: by1,by2,by3,by4,by5,by6,by7,by8
    real*8 :: bz1,bz2,bz3,bz4,bz5,bz6,bz7,bz8
    real*8 :: tenx,teny,tenz,xj,yj,zj,bxx,byy,bzz,btot,tjdotb &
                      ,curr_tot
    integer*8 :: i,j,k
    real*8 :: dbxdy,dbydx,dbzdx,dbxdz,dbzdy,dbydz
    real*8 :: curlbx_scalar,curlby_scalar,curlbz_scalar,bxav,byav,bzav

    do k = kb, ke 
      do j = jb,je
        do i = 2, nx1
          dbxdy= bx(i+1,j+1,k) + bx(i,j+1,k) - bx(i+1,j,k) - bx(i,j,k)
          dbxdz= 0.0
          dbydx= by(i+1,j+1,k) + by(i+1,j,k) - by(i,j+1,k) - by(i,j,k)
          dbydz= 0.0
          dbzdx= bz(i+1,j+1,k) + bz(i+1,j,k) - bz(i,j+1,k) - bz(i,j,k)
          dbzdy= bz(i+1,j+1,k) + bz(i,j+1,k) - bz(i+1,j,k) - bz(i,j,k)

          curlbx_scalar=dbzdy/(2.*meshY%dxc(j+1))-dbydz/(2.*meshZ%dxc(k+1))
          curlby_scalar=dbxdz/(2.*meshZ%dxc(k+1))-dbzdx/(2.*meshX%dxc(i  ))
          curlbz_scalar=dbydx/(2.*meshX%dxc(i  ))-dbxdy/(2.*meshY%dxc(j+1))
  
          ! 6/25/2006 New eta_par option: tensor eta
          bx1=bx(i+1,j+1,k)  
          bx2=bx(i  ,j+1,k)  
          bx3=bx(i  ,j  ,k)  
          bx4=bx(i+1,j  ,k)  
          by1=by(i+1,j+1,k)  
          by2=by(i  ,j+1,k)  
          by3=by(i  ,j  ,k)  
          by4=by(i+1,j  ,k)  
          bz1=bz(i+1,j+1,k)  
          bz2=bz(i  ,j+1,k)  
          bz3=bz(i  ,j  ,k)  
          bz4=bz(i+1,j  ,k)  

          bxav=.25*(bx1+bx2+bx3+bx4)
          byav=.25*(by1+by2+by3+by4)
          bzav=.25*(bz1+bz2+bz3+bz4)

          xj = curlbx_scalar
          yj = curlby_scalar
          zj = curlbz_scalar
          bxx = bxav
          byy = byav
          bzz = bzav
          xtmp1m = sqrt(bxx**2 + byy**2 + bzz**2)
          btot = max(1.d-12,xtmp1m)
          xtmp1m = sqrt(xj**2 + yj**2 + zj**2) 
          curr_tot = max(1.d-12,xtmp1m)
          
          if (eta_par.eq.0) then
            tenx=eta(i,j,k)*xj
            teny=eta(i,j,k)*yj
            tenz=eta(i,j,k)*zj
          else
            bxx = bxav
            byy = byav
            bzz = bzav
            xtmp1m = sqrt(bxx**2 + byy**2 + bzz**2)
            btot = max(1.d-12,xtmp1m)
            xtmp1m = sqrt(xj**2 + yj**2 + zj**2)
            curr_tot = max(1.d-12,xtmp1m)
      
            if (eta_par.eq.1) then
              tjdotb = eta(i,j,k)*(bxx*xj + byy*yj + bzz*zj)/btot
              tenx = tjdotb*bxx/btot
              teny = tjdotb*byy/btot
              tenz = tjdotb*bzz/btot
            else if (eta_par.eq.2) then
              tenx = abs(eta(i,j,k)*bxx*xj/(btot*curr_tot))
              tenx = min(resis,tenx)
              tenx = tenx*xj
              teny = abs(eta(i,j,k)*byy*yj/(btot*curr_tot))
              teny = min(resis,teny)
              teny = teny*yj
              tenz = abs(eta(i,j,k)*bzz*zj/(btot*curr_tot))
              tenz = min(resis,tenz)
              tenz = tenz*zj
            endif
            eta_times_b_dot_j(i,j,k) = min(resis,eta(i,j,k)*abs((bxav*xj + byav*yj + bzav*zj))/(curr_tot*btot))
          endif
          
          fox(i,j,k)=-tenx
          foy(i,j,k)=-teny
          foz(i,j,k)=-tenz

        enddo
      enddo
    enddo

    ! boundary conditions
    ! first update internal ghost cells so that all processors
    ! have latest information. Care must be exercised so that
    ! BCs are set wrt proecessor that corresponds to that location. 
    ! To that end, keep z loops set to limits of kb and ke.
    call xrealbcc_pack_e_2d(fox,foy,foz,1_8,nx,ny,nz)

    fox(1  ,:,:)=fox(nx1  ,:,:)
    foy(1  ,:,:)=foy(nx1  ,:,:)
    foz(1  ,:,:)=foz(nx1  ,:,:)
    fox(nx2,:,:)=fox(2,:,:)
    foy(nx2,:,:)=foy(2,:,:)
    foz(nx2,:,:)=foz(2,:,:)

    fox(:,:,kb-1)=fox(:,:,kb)
    foy(:,:,kb-1)=foy(:,:,kb)
    foz(:,:,kb-1)=foz(:,:,kb)
    
    fox(:,:,kb+1)=fox(:,:,kb)
    foy(:,:,kb+1)=foy(:,:,kb)
    foz(:,:,kb+1)=foz(:,:,kb)

    return
  end subroutine focalc_2d


  !---------------------------------------------------------------------
  ! retrun the value for field masking, dependent on z position
  !---------------------------------------------------------------------
  double precision function masking_func(k, flag)
    integer*8, intent(in) :: k
    integer, intent(in) :: flag

    ! if ( (mask .eqv. .true.) .and. (flag == 0 .or. flag == n_sub_b) ) then
    if ( (mask .eqv. .true.) .and. (flag == n_sub_b) ) then
      if ( k <= mask_zs ) then
        masking_func = 1.-(mask_r*(real(k)-mask_zs)/real(mask_zs))**2.
      else if ( k >= nz-mask_zs ) then
        masking_func = 1.-(mask_r*(real(k)-nz+mask_zs)/real(mask_zs))**2.
      else 
        masking_func = 1.
      endif
    else 
      masking_func = 1.
    endif 

    return
  end function masking_func


end module m_field